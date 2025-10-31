

DATA_DIR = config["output_base_dir"].rstrip("/")

if config["TRIM_ADAPTERS"]:
    # trim adapters from reads before alignment
    rule fastp:
        input:
            r1 = lambda wc: get_bowtie2_input(wc)[0],
            r2 = lambda wc: get_bowtie2_input(wc)[1]
        output:
            r1 = temp(f"{DATA_DIR}/fastp/{{sample}}_R1.fastq.gz"),
            r2 = temp(f"{DATA_DIR}/fastp/{{sample}}_R2.fastq.gz")
        params:
            adapter_fasta_file = config["ADAPTER_FASTA"]
        conda: 
            "../envs/fastp.yml"
        singularity:
            os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "fastp.sif")
        log:
            f"{DATA_DIR}/logs/fastp/{{sample}}.fastp.json"
        threads: 4
        shell:
            "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe --trim_poly_g --adapter_fasta {params.adapter_fasta_file} --thread {threads} -j {log} -h /dev/null"

# fastqc for each read 
rule fastqc:
    input:
        f"{DATA_DIR}/fastp/{{read}}.fastq.gz" if config["TRIM_ADAPTERS"] else get_read_path
    output:
        html=f"{DATA_DIR}/fastqc/{{read}}_fastqc.html",
        zip=f"{DATA_DIR}/fastqc/{{read}}_fastqc.zip"
    conda:
        "../envs/fastqc.yml"
    singularity:
        "docker://staphb/fastqc:0.11.9"
    log:
        f"{DATA_DIR}/logs/fastqc_{{read}}.log"
    threads: 1
    shell:
        f"fastqc -t {{threads}} --outdir {DATA_DIR}/fastqc {{input}} > {{log}} 2>&1"

# detect contaminants
rule fastq_screen:
    input:
        f"{DATA_DIR}/fastp/{{read}}.fastq.gz" if config["TRIM_ADAPTERS"] else get_read_path
    output:
        f"{DATA_DIR}/fastq_screen/{{read}}_screen.txt"
    conda:
        "../envs/fastq_screen.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "fastq_screen.sif")
    log:
        f"{DATA_DIR}/logs/fastq_screen_{{read}}.log"
    threads: 2
    shell:
        f"fastq_screen --aligner bowtie2 --threads {{threads}} --outdir {DATA_DIR}/fastq_screen "
        f"--conf {{config[FASTQ_SCREEN]}} --force {{input}} > {{log}} 2>&1"
