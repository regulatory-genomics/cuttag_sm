

DATA_DIR = config["output_base_dir"].rstrip("/")

if config["TRIM_ADAPTERS"]:
    # trim adapters from reads before alignment
    rule fastp:
        input:
            r1 = lambda wc: get_bowtie2_input_by_run(wc)[0],
            r2 = lambda wc: get_bowtie2_input_by_run(wc)[1]
        output:
            r1 = temp(f"{DATA_DIR}/middle_file/Trimmed_fastq/fastp/{{sample}}.{{run}}_R1.fastq.gz"),
            r2 = temp(f"{DATA_DIR}/middle_file/Trimmed_fastq/fastp/{{sample}}.{{run}}_R2.fastq.gz")
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
