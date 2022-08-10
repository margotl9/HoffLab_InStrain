# -----------------------------------------------------
# Snakefile for running InStrain
# -----------------------------------------------------
import pandas as pd
import os


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv("config/samples.tsv", sep="\t")
assembly_sample = samples_df["assembly"] + "_" + samples_df["sample"]


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# --------------------------------------------------
# Reads Rules
# --------------------------------------------------
# symlink input paths to new paths (abundance and sample)
rule symlink_reads:
    input:
        R1=lambda wildcards: samples_df[
            (samples_df["assembly"] + "_" + samples_df["sample"])
            == wildcards.assembly_sample
            ].iloc[0]["R1"],
        R2=lambda wildcards: samples_df[
            (samples_df["assembly"] + "_" + samples_df["sample"])
            == wildcards.assembly_sample
            ].iloc[0]["R2"],
    output:
        R1=results+"00_INPUT/{assembly_sample}_paired_1.fastq.gz",
        R2=results+"00_INPUT/{assembly_sample}_paired_2.fastq.gz",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """

# Align reads to mapping file using bowtie2
rule build_bowtie2db:
    input:
        config["mapping"]
    output:
        results+"00_INPUT/00_bowtie2/map.1.bt2",
    params:
        db=results+"00_INPUT/00_bowtie2/map", 
        extra_args=config["bowtie_build_extra_args"],
    conda:
        "./envs/bowtie.yml"
    threads: config["bowtie_threads"]
    shell:
        """
        # make a bowtie2 db from mapping
        bowtie2-build {input} \
        {params.db} \
        --threads {threads} \
        {params.extra_args}
        """



# Align reads to mapping file using bowtie2
rule bowtie2_align_reads:
    input:
        R1=results+"00_INPUT/{assembly_sample}_paired_1.fastq.gz",
        R2=results+"00_INPUT/{assembly_sample}_paired_2.fastq.gz",
        # R1S=results
        # + "01_READ_PREPROCESSING/03_kneaddata/{assembly_sample}_unmatched_1.fastq",
        # R2S=results
        # + "01_READ_PREPROCESSING/03_kneaddata/{assembly_sample}_unmatched_2.fastq",
        db=results+"00_INPUT/00_bowtie2/map.1.bt2",
    output:
        results
        + "00_INPUT/00_bowtie2/bam_files/{assembly_sample}.bam",
    params:
        db=results+"00_INPUT/00_bowtie2/map",
        sam=results+"00_INPUT/00_bowtie2/{assembly_sample}.sam",
        extra_args=config["bowtie_align_extra_args"],

    log:
        results+"00_INPUT/00_bowtie2/logs/{assembly_sample}_log",
    conda:
        "./envs/bowtie.yml"
    threads: config["bowtie_threads"]
    ##uncomment below and add to shell when using R1S and R2S 
    # -U {input.R1S},{input.R2S} \
    shell:
        """
        # align reads to bowtie2 database
        bowtie2 \
        --threads {threads} \
        -x {params.db} \
        -1 {input.R1} \
        -2 {input.R2} \
        -S {params.sam} \
        > {log} 2>&1 \
        {params.extra_args}
 
        # convert sam to bam
        samtools view -S -b {params.sam} > {output} 
        rm {params.sam}
        """


# --------------------------------------------------
# Bam Rules
# --------------------------------------------------
# symlink input paths to new paths
rule symlink_bam:
    input:
        bam=lambda wildcards: samples_df[
            (samples_df["assembly"] + "_" + samples_df["sample"])
            == wildcards.assembly_sample
            ].iloc[0]["bam"],
    output:
        bam=results+"00_INPUT/{assembly_sample}.bam",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.bam} {output.bam}
        """


#if input file is bam or reads
picked_bam = ""
if config['input_format'] == 'bam':
    picked_bam = results + "00_INPUT/{assembly_sample}.bam"
elif config['input_format'] == 'reads':
    picked_bam = results + "00_INPUT/00_bowtie2/bam_files/{assembly_sample}.bam"

# -----------------------------------------------------
# 01 InStrain profile
# -----------------------------------------------------
rule sort_bam:
    input: 
        bam = picked_bam
    output:
        sort = results+"sorted_bam_files/{assembly_sample}_sorted.bam",
    conda:
        "./envs/instrain.yml"
    shell:
        """
        samtools sort {input.bam} -o {output.sort}
        """

# rule reference_db_to_fna:
#     input:
#         config["gene_fasta"]
#     output:
#         results+"00_INPUT/ref_gene_contig.fna"
#     conda:
#         "./envs/instrain.yml"
#     shell:
#         """
#         prodigal -i {input} \
#         -o {output}
#         """

rule inStrain_profile:
    input:
        mapping_fasta=config["mapping"],
        bam=results+"sorted_bam_files/{assembly_sample}_sorted.bam",
        fna=config["gene_fna"]
    output:
        results+"01_instrain_profile/{assembly_sample}_profile.IS/output/{assembly_sample}_profile.IS_linkage.tsv",
    params:
        out_dir=results+"01_instrain_profile/{assembly_sample}_profile.IS",
        min_breadth=config["min_breadth"],
        extra_args=config["inStrain_profile_extra_args"]
    conda:
        "./envs/instrain.yml"
    shell:
        """
        inStrain profile \
        {input.bam} \
        {input.mapping_fasta} \
        -o {params.out_dir} \
        -g {input.fna} \
        --skip_mm_profiling \
        --min_cov {params.min_breadth} \
        {params.extra_args}
        """



# -----------------------------------------------------
# 02 InStrain compare
# -----------------------------------------------------
rule inStrain_compare:
    input:
        expand(results+"01_instrain_profile/{assembly_sample}_profile.IS/output/{assembly_sample}_profile.IS_linkage.tsv", assembly_sample=assembly_sample)
    output:
        results+"02_instrain_compare/compare.IS/output/compare.IS_comparisonsTable.tsv",
    params:
        bam=expand(results+"sorted_bam_files/{assembly_sample}_sorted.bam", assembly_sample=assembly_sample),
        profile=expand(results+"01_instrain_profile/{assembly_sample}_profile.IS/", assembly_sample=assembly_sample),
        out_dir=results+"02_instrain_compare/compare.IS",
        extra_args=config["inStrain_compare_extra_args"]

    conda:
        "./envs/instrain.yml"
    shell:
        """
        inStrain compare \
        -i {params.profile} \
        -o {params.out_dir} \
        {params.extra_args}
        """


rule all:
    input:
        expand(results+"01_instrain_profile/{assembly_sample}_profile.IS/output/{assembly_sample}_profile.IS_linkage.tsv", assembly_sample=assembly_sample)
