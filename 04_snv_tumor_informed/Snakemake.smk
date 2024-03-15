from os.path import join as opj
import pandas as pd
#configfile: "configs/config.yaml"

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
            opj(config["out_dir"],"overlap/{min_qual}/others/{sample_repeated}_{sample_repeated}_ALT.bam"),
            sample_repeated=config["samples_bam_repeated"].keys(),
            min_qual=config["min_qual"],
            ),
        expand(
            opj(config["out_dir"],"infer_tf/{min_qual}/others/{sample_repeated}_{sample_repeated}.txt"),
            sample_repeated=config["samples_bam_repeated"].keys(),
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
        ## TODO: add cfDNA TF inference step (integrate infer_tf.py)
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

        expand(
            opj(config["out_dir"],"infer_tf/{min_qual}/others/TF_summary/TF_max_{tf_max}_{sample_repeated}_{sample_repeated}_tf_summary.csv"),
            sample_repeated=config["samples_bam_repeated"].keys(),
            min_qual=config["min_qual"],
            tf_max=[1.0, 0.1, 0.01, 0.02, 0.001],
        ),

        ## TODO: Add these target files after sequencing_summary has been updated. Snakemake has been tested dry-run.
        expand(opj(config["out_dir"],"overlap/{min_qual}/realtime/{control}_{vcf}_SNV_realtime_{time_period}.pickle.gz"),
            control=config["healthy_controls_bam"].keys(),
            vcf=config["vcf"].keys(),
            min_qual=config["min_qual"],
            time_period=[360, 1200],
        ),
        expand(opj(config["out_dir"],"overlap/{min_qual}/realtime/{sample}_{sample}_SNV_realtime_{time_period}.pickle.gz"),
                       sample=config["samples_bam"].keys(),
                       # vcf=config["vcf"].keys(),
                       min_qual=config["min_qual"],
                       time_period = [360, 1200],

        ),

        expand(opj(config["out_dir"],"overlap/{min_qual}/others/realtime/{sample}_{sample}_SNV_realtime_{time_period}.pickle.gz"),
            sample=config["samples_bam_repeated"].keys(),
            min_qual=config["min_qual"],
            time_period=[360, 1200],

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
    params:
        min_qual = "{min_qual}"
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
    params:
        min_qual = "{min_qual}"
    shell:
        """
        python scripts/01_snv_reads.py --input_bam_path {input.bam} --somatic_vcf_path {input.somatic_vcf} \
        --germline_vcf_path {input.somatic_vcf} --phasing False --min_qual {params.min_qual} \
        --out_pickle {output.pickle} --all_variants {output.all} --out_bam_alt {output.ALT} \
        --out_bam_ref {output.REF} --out_bam_other {output.others}
        """


rule overlap_snv_repeated:
    input:
        bam = lambda wildcards: config["samples_bam_repeated"][wildcards.sample_repeated],
        somatic_vcf = lambda wildcards: config["vcf_repeated"][wildcards.sample_repeated]
    output:
        ALT = opj(config["out_dir"], "overlap/{min_qual}/others/{sample_repeated}_{sample_repeated}_ALT.bam"),
        REF = opj(config["out_dir"], "overlap/{min_qual}/others/{sample_repeated}_{sample_repeated}_REF.bam"),
        others= opj(config["out_dir"],"overlap/{min_qual}/others/{sample_repeated}_{sample_repeated}_others.bam"),
        pickle= opj(config["out_dir"],"overlap/{min_qual}/others/{sample_repeated}_{sample_repeated}_overlap.pickle.gz"),
        all = opj(config["out_dir"], "overlap/{min_qual}/others/{sample_repeated}_{sample_repeated}_all.pickle.gz"),
    params:
        min_qual = "{min_qual}"
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
    params:
        query_map_qual = 60,
        sample= "{sample}",
        vcf_tf = lambda wildcards: config["vcf_tf"][wildcards.vcf]
    shell:
        """
        python scripts/vaf_observe.py --query_map_qual {params.query_map_qual} --name {params.sample} \
        --pickle {input.pickle} --tissue_tf {params.vcf_tf} \
        --bed {output.bed} --vaf_observe {output.vaf} --vaf_before_adjust {output.vaf_raw}
        """

rule get_vaf_repeated:
    input:
        pickle = opj(config["out_dir"],"overlap/{min_qual}/others/{sample_repeated}_{sample_repeated}_overlap.pickle.gz"),
    output:
        bed = opj(config["out_dir"], "infer_tf/{min_qual}/others/{sample_repeated}_{sample_repeated}_ALT.bed"),
        vaf = opj(config["out_dir"], "infer_tf/{min_qual}/others/{sample_repeated}_{sample_repeated}.txt"),
        vaf_raw= opj(config["out_dir"],"infer_tf/{min_qual}/{sample_repeated}_{sample_repeated}_raw.txt"),

    params:
        query_map_qual = 60,
        sample_repeated= "{sample_repeated}",
        vcf_tf = lambda wildcards: config["vcf_tf_repeated"][wildcards.sample_repeated]
    shell:
        """
        python scripts/vaf_observe.py --query_map_qual {params.query_map_qual} --name {params.sample_repeated} \
        --pickle {input.pickle} --tissue_tf {params.vcf_tf} \
        --bed {output.bed} --vaf_observe {output.vaf} --vaf_before_adjust {output.vaf_raw}
        """


rule get_cfdna_tf:
    input:
       txt = opj(config["out_dir"], "infer_tf/{min_qual}/{sample}_{vcf}.txt"),
    output:
        summary = opj(config["out_dir"], "infer_tf/{min_qual}/TF_summary/TF_max_{tf_max}_{sample}_{vcf}_tf_summary.csv"),
        plot = opj(config["out_dir"], "infer_tf/{min_qual}/TF_plot/TF_max_{tf_max}_{sample}_{vcf}_tf_montecarlo.png"),
    params:
        error_rates = config['error_rates'],
        name= "{sample}_{vcf}",
        output_folder = opj(config["out_dir"], "infer_tf/{min_qual}/"),
        tf_max = "{tf_max}",
        tf_sim_type = 'linear',
    shell:
        """
           python scripts/infer_tf.py -e {params.error_rates} -i {input.txt} -n {params.name} -s {output.summary} -f {output.plot} \
           -o {params.output_folder} -m {params.tf_max} --type {params.tf_sim_type}
        """



rule get_cfdna_tf_repeated:
    input:
       txt = opj(config["out_dir"], "infer_tf/{min_qual}/others/{sample_repeated}_{sample_repeated}.txt"),
    output:
        summary = opj(config["out_dir"], "infer_tf/{min_qual}/others/TF_summary/TF_max_{tf_max}_{sample_repeated}_{sample_repeated}_tf_summary.csv"),
        plot = opj(config["out_dir"], "infer_tf/{min_qual}/others/TF_plot/TF_max_{tf_max}_{sample_repeated}_{sample_repeated}_tf_montecarlo.png"),
    params:
        error_rates = config['error_rates'],
        name= "{sample_repeated}",
        output_folder= opj(config["out_dir"],"infer_tf/{min_qual}/others/"),
        tf_max = "{tf_max}",
        tf_sim_type = 'linear',
    shell:
        """
           python scripts/infer_tf.py -e {params.error_rates} -i {input.txt} -n {params.name} -s {output.summary} -f {output.plot} \
           -o {params.output_folder}  -m {params.tf_max} --type {params.tf_sim_type}
        """



rule get_timestamp:
    input:
        # vcf_name = "{sample}",
        bam_ref = opj(config["out_dir"],"overlap/{min_qual}/{control}_{vcf}_REF.bam"),
        bam_alt = opj(config["out_dir"],"overlap/{min_qual}/{control}_{vcf}_ALT.bam"),
        summary_file=lambda wildcards: config["summary_file"][wildcards.control]  # Update this path
    params:
        vcf = "{vcf}",
        time_period = "{time_period}",

    output:
        pickle = opj(config["out_dir"], "overlap/{min_qual}/realtime/{control}_{vcf}_SNV_realtime_{time_period}.pickle.gz"),
        fig= opj(config["out_dir"],"overlap/{min_qual}/figures/{control}_{vcf}_{time_period}_SNV_count.pdf")

    shell:
        """
        python scripts/realtime.py --sample {wildcards.control} --sequencing_summary_file {input.summary_file} \
        --vcf_name {params.vcf}  --bam_ref {input.bam_ref} --bam_alt {input.bam_alt} \
        --output_pickle {output.pickle} --output_figure {output.fig} --time_period {params.time_period}
        """

rule get_timestamp_repeated:
    input:
        bam_ref = opj(config["out_dir"],"overlap/{min_qual}/others/{sample_repeated}_{sample_repeated}_REF.bam"),
        bam_alt = opj(config["out_dir"],"overlap/{min_qual}/others/{sample_repeated}_{sample_repeated}_ALT.bam"),
        summary_file= lambda wildcards: config["summary_file"][wildcards.sample_repeated]  # Update this path

    params:
        vcf="{sample_repeated}",
        time_period= "{time_period}",
    output:
        pickle = opj(config["out_dir"], "overlap/{min_qual}/others/realtime/{sample_repeated}_{sample_repeated}_SNV_realtime_{time_period}.pickle.gz"),
        fig= opj(config["out_dir"],"overlap/{min_qual}/figures/{sample_repeated}_{sample_repeated}_{time_period}_SNV_count.pdf")
    shell:
        """
        python scripts/realtime.py --sample {wildcards.sample_repeated} --sequencing_summary_file {input.summary_file} \
        --vcf_name {params.vcf}  --bam_ref {input.bam_ref} --bam_alt {input.bam_alt} \
        --output_pickle {output.pickle} --output_figure {output.fig} --time_period {params.time_period}
        """


rule plot_raw_snv:
    input:
        snv_pickle = opj(config["out_dir"],"overlap/{min_qual}/{sample}_{sample}_overlap.pickle.gz"),
        vcf_pickle = opj(config["out_dir"], "overlap/{min_qual}/{sample}_{sample}_all.pickle.gz"),
        healthy_pickles= expand(opj(config["out_dir"],"overlap/{{min_qual}}/{control}_{{sample}}_overlap.pickle.gz"),
            control=config['healthy_controls_bam'].keys())

    output:
        figure_raw = opj(config["out_dir"], "figure2/{sample}_raw_snv_{min_qual}.pdf"),
    params:
        paths_to_other_figures = opj(config["out_dir"],"figure2/"),
    shell:
        """
        python scripts/plot_raw_snv.py -f {output.figure_raw} -p {params.paths_to_other_figures} -v {input.vcf_pickle} \
         {input.snv_pickle} {input.healthy_pickles}
        """
