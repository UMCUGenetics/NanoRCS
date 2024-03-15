from os.path import join as opj
import pandas as pd
import os

wildcard_constraints:
    part2 = '[a-zA-Z0-9]+'

wildcard_constraints:
    part1 = '[a-zA-Z0-9-]+'

out_dir = config['out_dir']
vaf_table_path = config['vaf_table_path']
cancer_type = config['cancer_type']
sample_names,  = glob_wildcards(f'{vaf_table_path}/{{cancer_type_sample_name}}.pickle.gz')
print(sample_names)
#cancer_type = sample_names[0].split('_')[0]
#print(cancer_type)
GEs = config['GEs']
error_rates = config['error_rates']


def get_resource_time(wildcards, attempt):
    ge = float(wildcards.ge)
    hr = round(ge * 1.5,2)
    hr = hr * attempt * 8 + 2
    hr_print = str(hr) + 'h'
    #print(hr_print)
    return hr_print
    

def get_resource_mem(wildcards, attempt):
    mem = 32000* attempt
    return mem


rule all:
    input:
        opj(out_dir, f"combined/{cancer_type}_combined.pickle.gz"),

rule simulate:
    input:
        vcf_table = opj(config["vaf_table_path"], "{cancer_type_sample_name}.pickle.gz"),
    output:
        pickle = opj(out_dir, f"pickle/{cancer_type}/{{cancer_type_sample_name}}_{{ge}}_{{error_rate}}.pickle.gz"),
    params:
        path_pickle = opj(out_dir, f"pickle/{cancer_type}"),
        ge = "{ge}",
        error_rate = "{error_rate}",
        tfmin = config['tumor_fraction_min'],
        tfmax = config['tumor_fraction_max'],
        tumor_fraction_steps = config['tumor_fraction_step'],
        step_type= config['step_type'],
    log:
        opj(out_dir,f"log/{cancer_type}/{{cancer_type_sample_name}}_{{ge}}_{{error_rate}}.log"),
    retries: 3
    conda:
        "envs/simulate.yaml"
    resources:
        runtime=get_resource_time,
        cpus = 2,
        mem_mb = get_resource_mem,
    shell:
        """
        python3 ./scripts/complete_simulation_smk_wrapper.py -i {input.vcf_table} -op {output.pickle} \
        -g {params.ge} -e {params.error_rate} -tfmin {params.tfmin} -tfmax {params.tfmax} \
        -s {params.tumor_fraction_steps} -tft {params.step_type} > {log}
        """

rule concat_pickle:
    input:
        pickle = expand(opj(out_dir,f"pickle/{cancer_type}/{{cancer_type_sample_name}}_{{ge}}_{{error_rate}}.pickle.gz"),ge=GEs,error_rate=error_rates,cancer_type_sample_name=sample_names),
    output:
        combined_pickle=opj(out_dir,f"combined/{cancer_type}_combined.pickle.gz"),
    params:
        pickle_path = opj(out_dir,f"pickle/{cancer_type}")
    log:
        opj(out_dir,f"log/{cancer_type}_combine_pickle.log"),
    resources:
        runtime='1h',
        cpus=2,
        mem_mb=16000,
    shell:
        """
        python scripts/combine_pickle.py -i {params.pickle_path} -o {output.combined_pickle}
        """

