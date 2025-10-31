
DATA_DIR = config["output_base_dir"].rstrip("/")

# align samples to genome
rule bowtie2:
    input:
        r1 = (lambda wc: f"{DATA_DIR}/fastp/{wc.sample}_R1.fastq.gz") if config["TRIM_ADAPTERS"] else (lambda wc: get_bowtie2_input(wc)[0]),
        r2 = (lambda wc: f"{DATA_DIR}/fastp/{wc.sample}_R2.fastq.gz") if config["TRIM_ADAPTERS"] else (lambda wc: get_bowtie2_input(wc)[1])
    output:
        f"{DATA_DIR}/aligned/{{sample}}.bam"
    log:
        err=f"{DATA_DIR}/logs/bowtie2_{{sample}}.err"
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

rule sort:
    input:
        f"{DATA_DIR}/aligned/{{sample}}.bam"
    output: 
        temp(f"{DATA_DIR}/aligned/{{sample}}.sort.bam")
    conda:
        "../envs/sambamba.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "sambamba.sif")
    threads: 4
    log:
        f"{DATA_DIR}/logs/sambamba_sort_{{sample}}.log"
    shell:
        f"sambamba sort --tmpdir={DATA_DIR}/aligned -t {{threads}} -o {{output}} {{input}} > {{log}} 2>&1"

rule markdup:
    input:
        rules.sort.output
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
