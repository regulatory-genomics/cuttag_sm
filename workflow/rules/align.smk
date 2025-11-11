
DATA_DIR = config["output_base_dir"].rstrip("/")

# Conditional aligner selection
if config.get("aligner", "bowtie2") == "bowtie2":
    # align runs to genome (per-run) using bowtie2
    rule align:
        input:
            r1 = (lambda wc: f"{DATA_DIR}/fastp/{wc.sample}.{wc.run}_R1.fastq.gz") if config["TRIM_ADAPTERS"] else (lambda wc: get_bowtie2_input_by_run(wc)[0]),
            r2 = (lambda wc: f"{DATA_DIR}/fastp/{wc.sample}.{wc.run}_R2.fastq.gz") if config["TRIM_ADAPTERS"] else (lambda wc: get_bowtie2_input_by_run(wc)[1])
        output:
            f"{DATA_DIR}/aligned/{{sample}}.{{run}}.bam"
        log:
            err=f"{DATA_DIR}/logs/bowtie2_{{sample}}.{{run}}.err"
        conda:
            "../envs/align.yml"
        singularity:
            os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "align.sif")
        threads: 8
        shell:
            "bowtie2 --local --very-sensitive-local "
            "--no-unal --no-mixed --threads {threads} "
            "--no-discordant --phred33 "
            "-I 10 -X 700 -x {config[GENOME]} "
            "-1 {input.r1} -2 {input.r2} 2>{log.err} | samtools view -@ {threads} -Sbh - > {output}"
else:
    # align runs to genome (per-run) using bwa-mem2
    rule align:
        input:
            r1 = (lambda wc: f"{DATA_DIR}/fastp/{wc.sample}.{wc.run}_R1.fastq.gz") if config["TRIM_ADAPTERS"] else (lambda wc: get_bowtie2_input_by_run(wc)[0]),
            r2 = (lambda wc: f"{DATA_DIR}/fastp/{wc.sample}.{wc.run}_R2.fastq.gz") if config["TRIM_ADAPTERS"] else (lambda wc: get_bowtie2_input_by_run(wc)[1])
        output:
            f"{DATA_DIR}/aligned/{{sample}}.{{run}}.bam"
        log:
            err=f"{DATA_DIR}/logs/bwa_mem2_{{sample}}.{{run}}.err"
        params:
            rg = (lambda wc: f"@RG\\tID:{wc.sample}.{wc.run}\\tSM:{wc.sample}\\tPL:{config.get('sequencing_platform', 'ILLUMINA')}"),
            ref = (lambda wc: config["FASTA"]),
            bwa_args = config.get("bwa_args", "")
        resources:
            mem_mb=config.get("mem", "64000")
        conda:
            "../envs/bwa.yml"
        singularity:
            os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "bwa.sif")
        threads: 8
        shell:
            "bwa-mem2 mem {params.bwa_args} -t {threads} -R \"{params.rg}\" {params.ref} "
            "{input.r1} {input.r2} 2>{log.err} | samtools view -@ {threads} -Sbh - > {output}"

rule sort:
    input:
        f"{DATA_DIR}/aligned/{{sample}}.{{run}}.bam"
    output: 
        temp(f"{DATA_DIR}/aligned/{{sample}}.{{run}}.sort.bam")
    conda:
        "../envs/sambamba.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "sambamba.sif")
    threads: 4
    log:
        f"{DATA_DIR}/logs/sambamba_sort_{{sample}}.{{run}}.log"
    shell:
        f"sambamba sort --tmpdir={DATA_DIR}/aligned -t {{threads}} -o {{output}} {{input}} > {{log}} 2>&1"

rule merge_runs:
    input:
        get_sorted_bams_for_sample
    output:
        temp(f"{DATA_DIR}/aligned/{{sample}}.merged.bam")
    conda:
        "../envs/align.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "sambamba.sif")
    threads: 4
    log:
        f"{DATA_DIR}/logs/samtools_merge_{{sample}}.log"
    shell:
        "samtools merge -@ {threads} {output} {input} > {log} 2>&1"

rule markdup:
    input:
        rules.merge_runs.output
    output:
        f"{DATA_DIR}/markd/{{sample}}.sorted.markd.bam"
    conda:
        "../envs/sambamba.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "sambamba.sif")
    threads: 4
    log:
        f"{DATA_DIR}/logs/sambamba_markdup_{{sample}}.log"
    shell:
        f"sambamba markdup --tmpdir={DATA_DIR}/markd -t {{threads}} {{input}} {{output}} > {{log}} 2>&1"

rule index:
    input:
        rules.markdup.output
    output:
        f"{DATA_DIR}/markd/{{sample}}.sorted.markd.bam.bai"
    conda:
        "../envs/sambamba.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "sambamba.sif")
    threads: 2
    log:
        f"{DATA_DIR}/logs/samtools_index_{{sample}}.log"
    shell:
        "sambamba index -t {threads} {input} > {log} 2>&1"
