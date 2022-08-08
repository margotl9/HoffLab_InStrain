# -----------------------------------------------------
# Snakefile for running InStrain
# -----------------------------------------------------
import pandas as pd
import os


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv("config/samples.tsv", sep="\t")
sample = samples_df["sample"]


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# --------------------------------------------------
# Rules
# --------------------------------------------------
# symlink input paths to new paths
rule symlink_bam:
    input:
        bam=lambda wildcards: samples_df[
            (samples_df["sample"])
            == wildcards.sample
            ].iloc[0]["bam"],
    output:
        bam=results+"00_INPUT/{sample}.bam",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.bam} {output.bam}
        """

# -----------------------------------------------------
# 01 InStrain profile
# -----------------------------------------------------
rule sort_bam:
    input: 
        bam = results+"00_INPUT/{sample}.bam",
    output:
        sort = results+"sorted_bam_files/{sample}_sorted.bam",
    conda:
        "./envs/instrain.yml"
    shell:
        """
        samtools sort {input.bam} -o {output.sort}
        """

rule reference_db_to_fna:
    input:
        config["ref_db"]
    output:
        resources+"ref_contig.fna"
    conda:
        "../envs/instrain.yml"
    shell:
        """
        prodigal -i {input} \
        -o {output}
        """

rule inStrain_profile:
    input:
        fasta=config["ref_db"],
        bam=results+"sorted_bam_files/{sample}_sorted.bam",
    output:
        results+"01_instrain_profile/{sample}_profile.IS/output/{sample}_profile.IS_linkage.tsv",
    params:
        out_dir=results+"01_instrain_profile/{sample}_profile.IS",
        min_breadth=config["min_breadth"],
        extra_args=config["inStrain_profile_extra_args"]

    conda:
        "./envs/instrain.yml"
    # -g {input.fna} \
    shell:
        """
        inStrain profile \
        {input.bam} \
        {input.fasta} \
        -o {params.out_dir} \
        --skip_mm_profiling \
        --skip_genome_wide \
        --min_cov {params.min_breadth} \
        {params.extra_args}
        """



# -----------------------------------------------------
# 02 InStrain compare
# -----------------------------------------------------
rule inStrain_compare:
    input:
        expand(results+"01_instrain_profile/{sample}_profile.IS/output/{sample}_profile.IS_linkage.tsv", sample=sample)
    output:
        results+"02_instrain_compare/compare.IS/output/compare.IS_comparisonsTable.tsv",
    params:
        bam=expand(results+"sorted_bam_files/{sample}_sorted.bam", sample=sample),
        profile=expand(results+"01_instrain_profile/{sample}_profile.IS", sample=sample),
        out_dir=results+"02_instrain_compare/compare.IS",
        extra_args=config["inStrain_compare_extra_args"]

    conda:
        "./envs/instrain.yml"
    shell:
        """
        inStrain compare \
        -i {params.profile} \
        -o {params.out_dir} \
        -bams {params.bam} \
        {params.extra_args}
        """


rule all:
    input:
        results+"02_instrain_compare/compare.IS/output/compare.IS_comparisonsTable.tsv"