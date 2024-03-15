from os.path import join as opj
import pandas as pd
configfile: "configs/config_ILL.yaml"

wildcard_constraints:
    min_qual = "[0-9]+",
    sample = "[a-zA-Z0-9]+",
    vcf = "[a-zA-Z0-9]+",

rule all:
    input:
        expand(
            opj(config["out_dir"], "overlap/{min_qual}/{control}_{vcf}_ALT.bam"),
            control=config["healthy_controls_bam"].keys(),
            vcf=config["vcf"].keys(),
            min_qual=config["min_qual"],
        ),
        expand(
            opj(config["out_dir"],"overlap/{min_qual}/{sample}_{sample}_ALT.bam"),
            sample=config["samples_bam"].keys(),
            min_qual=config["min_qual"],
            ),
        expand(
            opj(config["out_dir"],"infer_tf/{min_qual}/{sample}_{sample}.txt"),
            sample=config["samples_bam"].keys(),
            min_qual=config["min_qual"],
        ),

        expand(
            opj(config["out_dir"],"infer_tf/{min_qual}/{control}_{vcf}.txt"),
            control=config["healthy_controls_bam"].keys(),
            vcf=config["vcf"].keys(),
            min_qual=config["min_qual"],
        ),
        expand(
            opj(config["out_dir"],"infer_tf/{min_qual}/TF_summary/TF_max_{tf_max}_{control}_{vcf}_tf_summary.csv"),
            control=config["healthy_controls_bam"].keys(),
            vcf=config["vcf"].keys(),
            min_qual=config["min_qual"],
            tf_max=[1.0, 0.1, 0.01, 0.001],
            ),
        expand(
            opj(config["out_dir"],"infer_tf/{min_qual}/TF_summary/TF_max_{tf_max}_{sample}_{sample}_tf_summary.csv"),
            sample=config["samples_bam"].keys(),
            min_qual=config["min_qual"],
            tf_max=[1.0, 0.1, 0.01, 0.001],
            ),
        expand(opj(config["out_dir"], "figure2/{sample}_raw_snv_{min_qual}.pdf"),
            sample=config["samples_bam"].keys(),
            # vcf=config["vcf"].keys(),
            min_qual=config["min_qual"],
        ),


rule overlap_snv_hc:
    input:
        bam = lambda wildcards: config["healthy_controls_bam"][wildcards.control],
        somatic_vcf = lambda wildcards: config["vcf"][wildcards.vcf]
    output:
        ALT = opj(config["out_dir"], "overlap/{min_qual}/{control}_{vcf}_ALT.bam"),
        REF = opj(config["out_dir"], "overlap/{min_qual}/{control}_{vcf}_REF.bam"),
        others= opj(config["out_dir"],"overlap/{min_qual}/{control}_{vcf}_others.bam"),
        pickle= opj(config["out_dir"],"overlap/{min_qual}/{control}_{vcf}_overlap.pickle.gz"),
        all = opj(config["out_dir"], "overlap/{min_qual}/{control}_{vcf}_all.pickle.gz"),
    conda:
        "envs/vis.yaml"
    params:
        min_qual = "{min_qual}"
    resources:
        cpus=2,
        disk_mb=2000,
        mem_mb=1000,
        runtime='1h'
    shell:
        """
        python scripts/01_snv_reads.py --input_bam_path {input.bam} --somatic_vcf_path {input.somatic_vcf} \
        --germline_vcf_path {input.somatic_vcf} --phasing False --min_qual {params.min_qual} \
        --out_pickle {output.pickle} --all_variants {output.all} --out_bam_alt {output.ALT} \
        --out_bam_ref {output.REF} --out_bam_other {output.others}
        """

rule overlap_snv_tumor:
    input:
        bam = lambda wildcards: config["samples_bam"][wildcards.sample],
        somatic_vcf = lambda wildcards: config["vcf"][wildcards.sample]
    output:
        ALT = opj(config["out_dir"], "overlap/{min_qual}/{sample}_{vcf}_ALT.bam"),
        REF = opj(config["out_dir"], "overlap/{min_qual}/{sample}_{vcf}_REF.bam"),
        others= opj(config["out_dir"],"overlap/{min_qual}/{sample}_{vcf}_others.bam"),
        pickle= opj(config["out_dir"],"overlap/{min_qual}/{sample}_{vcf}_overlap.pickle.gz"),
        all = opj(config["out_dir"], "overlap/{min_qual}/{sample}_{vcf}_all.pickle.gz"),
    conda:
        "envs/vis.yaml"
    params:
        min_qual = "{min_qual}"
    resources:
        cpus=2,
        disk_mb=2000,
        mem_mb=1000,
        runtime='1h'
    shell:
        """
        python scripts/01_snv_reads.py --input_bam_path {input.bam} --somatic_vcf_path {input.somatic_vcf} \
        --germline_vcf_path {input.somatic_vcf} --phasing False --min_qual {params.min_qual} \
        --out_pickle {output.pickle} --all_variants {output.all} --out_bam_alt {output.ALT} \
        --out_bam_ref {output.REF} --out_bam_other {output.others}
        """


rule get_vaf:
    input:
        pickle = opj(config["out_dir"],"overlap/{min_qual}/{sample}_{vcf}_overlap.pickle.gz"),
    output:
        bed = opj(config["out_dir"], "infer_tf/{min_qual}/{sample}_{vcf}_ALT.bed"),
        vaf = opj(config["out_dir"], "infer_tf/{min_qual}/{sample}_{vcf}.txt"),
        vaf_raw = opj(config["out_dir"], "infer_tf/{min_qual}/{sample}_{vcf}_raw.txt"),
    conda:
        "envs/vis.yaml"
    params:
        query_map_qual = 60,
        sample= "{sample}",
        vcf_tf = lambda wildcards: config["vcf_tf"][wildcards.vcf]
    resources:
        cpus=2,
        disk_mb=2000,
        mem_mb=1000,
        runtime='1h'
    shell:
        """
        python scripts/vaf_observe.py --query_map_qual {params.query_map_qual} --name {params.sample} \
        --pickle {input.pickle} --tissue_tf {params.vcf_tf} \
        --bed {output.bed} --vaf_observe {output.vaf} --vaf_before_adjust {output.vaf_raw}
        """

rule get_cfdna_tf:
    input:
       txt = opj(config["out_dir"], "infer_tf/{min_qual}/{sample}_{vcf}.txt"),
    output:
        summary = opj(config["out_dir"], "infer_tf/{min_qual}/TF_summary/TF_max_{tf_max}_{sample}_{vcf}_tf_summary.csv"),
        plot = opj(config["out_dir"], "infer_tf/{min_qual}/TF_plot/TF_max_{tf_max}_{sample}_{vcf}_tf_montecarlo.pdf"),
    conda:
        "envs/vis.yaml"
    params:
        error_rates = config['error_rates'],
        name= "{sample}_{vcf}",
        output_folder = opj(config["out_dir"], "infer_tf/{min_qual}/"),
    resources:
        cpus=2,
        disk_mb=4000,
        mem_mb=8000,
        runtime='4h'
    shell:
        """
           python scripts/infer_tf.py -e {params.error_rates} -i {input.txt} -n {params.name} -s {output.summary} -f {output.plot} \
           -o {params.output_folder}
        """


rule get_cfdna_tf_MRD:
    input:
       txt = opj(config["out_dir"], "infer_tf/{min_qual}/{sample}_{vcf}.txt"),
    output:
        summary = opj(config["out_dir"], "infer_tf/{min_qual}/MRD_{tf_max}_{sample}_{vcf}_tf_summary.csv"),
        plot = opj(config["out_dir"], "infer_tf/{min_qual}/MRD_{tf_max}_{sample}_{vcf}_tf_montecarlo.png"),
    conda:
        "envs/vis.yaml"
    params:
        error_rates = config['error_rates'],
        name= "{sample}_{vcf}",
        output_folder = opj(config["out_dir"], "infer_tf/{min_qual}/"),
        tf_max = "{tf_max}",
        tf_sim_type = 'linear',
    resources:
        cpus=2,
        disk_mb=4000,
        mem_mb=8000,
        runtime='4h'
    shell:
        """
           python scripts/infer_tf.py -e {params.error_rates} -i {input.txt} -n {params.name} -s {output.summary} -f {output.plot} \
           -o {params.output_folder} -m {params.tf_max} --type {params.tf_sim_type}
        """

rule plot_raw_snv:
    input:
        snv_pickle = opj(config["out_dir"],"overlap/{min_qual}/{sample}_{sample}_overlap.pickle.gz"),
        vcf_pickle = opj(config["out_dir"], "overlap/{min_qual}/{sample}_{sample}_all.pickle.gz"),
        healthy_pickles= expand(opj(config["out_dir"],"overlap/{{min_qual}}/{control}_{{sample}}_overlap.pickle.gz"),
            control=config['healthy_controls_bam'].keys())

    output:
        figure_raw = opj(config["out_dir"], "figure2/{sample}_raw_snv_{min_qual}.pdf"),
        # figure_bar= opj(config["out_dir"],"figure2/{sample}_snv_summary_{min_qual}.pdf"),
        # figure_found_vaf= opj(config["out_dir"],"figure2/{sample}_{vcf}_tumor_vaf_after_adjust_{min_qual}.pdf"),
        # figure_vaf= opj(config["out_dir"],"figure2/{sample}_{vcf}_tumor_vaf_after_adjust_{min_qual}.pdf"),
    conda:
        "envs/vis.yaml"
    params:
        paths_to_other_figures = opj(config["out_dir"],"figure2/"),
    resources:
        cpus=2,
        disk_mb=2000,
        mem_mb=1000,
        runtime='1h'
    shell:
        """
        python scripts/plot_raw_snv.py -f {output.figure_raw} -p {params.paths_to_other_figures} -v {input.vcf_pickle} \
         {input.snv_pickle} {input.healthy_pickles}
        """
