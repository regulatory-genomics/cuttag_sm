

DATA_DIR = config["output_base_dir"].rstrip("/")

if config["TRIM_ADAPTERS"]:
    # trim adapters from reads before alignment
    rule fastp:
        input:
            r1 = lambda wc: get_bowtie2_input_by_run(wc)[0],
            r2 = lambda wc: get_bowtie2_input_by_run(wc)[1]
        output:
            r1 = temp(f"{DATA_DIR}/fastp/{{sample}}.{{run}}_R1.fastq.gz"),
            r2 = temp(f"{DATA_DIR}/fastp/{{sample}}.{{run}}_R2.fastq.gz")
        params:
            adapter_fasta_file = config["ADAPTER_FASTA"]
        conda: 
            "../envs/fastp.yml"
        singularity:
            os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "fastp.sif")
        log:
            f"{DATA_DIR}/logs/fastp/{{sample}}.{{run}}.fastp.json"
        threads: 4
        shell:
            "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe --trim_poly_g --adapter_fasta {params.adapter_fasta_file} --thread {threads} -j {log} -h /dev/null"

# fastqc for each sample/run/mate
rule fastqc:
    input:
        lambda wc: (f"{DATA_DIR}/fastp/{wc.sample}.{wc.run}_{wc.mate}.fastq.gz" if config["TRIM_ADAPTERS"] else (get_bowtie2_input_by_run(wc)[0] if wc.mate == "R1" else get_bowtie2_input_by_run(wc)[1]))
    output:
        html=f"{DATA_DIR}/fastqc/{{sample}}.{{run}}_{{mate}}_fastqc.html",
        zip=f"{DATA_DIR}/fastqc/{{sample}}.{{run}}_{{mate}}_fastqc.zip"
    conda:
        "../envs/fastqc.yml"
    singularity:
        "docker://staphb/fastqc:0.11.9"
    log:
        f"{DATA_DIR}/logs/fastqc_{{sample}}.{{run}}_{{mate}}.log"
    threads: 1
    shell:
        f"fastqc -t {{threads}} --outdir {DATA_DIR}/fastqc {{input}} > {{log}} 2>&1"

# detect contaminants per sample/run/mate
rule fastq_screen:
    input:
        lambda wc: (f"{DATA_DIR}/fastp/{wc.sample}.{wc.run}_{wc.mate}.fastq.gz" if config["TRIM_ADAPTERS"] else (get_bowtie2_input_by_run(wc)[0] if wc.mate == "R1" else get_bowtie2_input_by_run(wc)[1]))
    output:
        f"{DATA_DIR}/fastq_screen/{{sample}}.{{run}}_{{mate}}_screen.txt"
    conda:
        "../envs/fastq_screen.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "fastq_screen.sif")
    log:
        f"{DATA_DIR}/logs/fastq_screen_{{sample}}.{{run}}_{{mate}}.log"
    threads: 2
    shell:
        f"fastq_screen --aligner bowtie2 --threads {{threads}} --outdir {DATA_DIR}/fastq_screen "
        f"--conf {{config[FASTQ_SCREEN]}} --force {{input}} > {{log}} 2>&1"
