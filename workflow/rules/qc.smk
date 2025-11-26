
DATA_DIR = config["output_base_dir"].rstrip("/")
rule plotFinger:
    input:
        f"{DATA_DIR}/markd/{{sample}}.sorted.markd.bam", f"{DATA_DIR}/markd/{{sample}}.sorted.markd.bam.bai"
    output:
        f"{DATA_DIR}/dtools/fingerprint_{{sample}}.tsv"
    conda:
        "../envs/dtools.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "dtools.sif")
    log:
        f"{DATA_DIR}/logs/fingerprint_{{sample}}.log"
    shell:
        "plotFingerprint -b {input[0]} --smartLabels --outRawCounts {output} > {log} 2>&1"

rule frip_plot:
    input:
        expand(f"{DATA_DIR}/plotEnrichment/frip_{{sample}}.tsv", sample = sample_noigg)
    output:
        f"{DATA_DIR}/plotEnrichment/frip.html"
    #conda:
    #    "../envs/plot_report.yml"
    script:
        "../src/frip_plot.py"


rule preseq:
    input:
       rules.markdup.output
    output:
        f"{DATA_DIR}/preseq/estimates_{{sample}}.txt"
    resources:
        defect_mode = defect_mode
    conda:
        "../envs/preseq.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "preseq.sif")
    log:
        f"{DATA_DIR}/logs/preseq_{{sample}}.log"
    shell:
        "preseq c_curve -B {resources.defect_mode} -l 1000000000 -o {output} {input} > {log} 2>&1"

rule preseq_lcextrap:
    input:
        rules.markdup.output
    output:
        f"{DATA_DIR}/preseq/lcextrap_{{sample}}.txt"
    resources:
        defect_mode = defect_mode
    conda:
        "../envs/preseq.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "preseq.sif")
    log:
        f"{DATA_DIR}/logs/preseq_{{sample}}.log"
    shell:
        "preseq lc_extrap -B {resources.defect_mode} -l 1000000000 -e 1000000000 -o {output} {input} > {log} 2>&1 || (echo 'preseq lc_extrap failed; creating empty output' >> {log} 2>&1; : > {output})"

rule multiqc:
    input:
        #get_fastqc_outputs(),
        #get_fastq_screen_outputs(),
        expand(f"{DATA_DIR}/plotEnrichment/frip_{{sample}}.tsv", sample=sample_noigg),
        expand(f"{DATA_DIR}/preseq/lcextrap_{{sample}}.txt", sample=samps)
    output:
        f"{DATA_DIR}/multiqc/multiqc_report.html",
        f"{DATA_DIR}/multiqc/multiqc_data/multiqc_data.json"
    conda:
        "../envs/multiqc.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "multiqc.sif")
    log:
        f"{DATA_DIR}/logs/multiqc.log"
    shell:
        # comment out the "export ..." line if not running pipeline through Singularity
        f"export LC_ALL=C.UTF-8; export LANG=C.UTF-8; "
        f"multiqc "
        #f"{DATA_DIR}/fastqc "
        #f"{DATA_DIR}/fastq_screen "
        f"{DATA_DIR}/plotEnrichment "
        f"{DATA_DIR}/callpeaks "
        f"{DATA_DIR}/preseq "
        f"{DATA_DIR}/logs "
        f"-f -c workflow/src/multiqc_conf.yml -o {DATA_DIR}/multiqc "
        f"--ignore {DATA_DIR}/homer "
        f"--ignore {DATA_DIR}/aligned "
        f"--ignore {DATA_DIR}/markd "
        f"--ignore {DATA_DIR}/tracks "
        f"--ignore {DATA_DIR}/counts "
        f"> {{log}} 2>&1"


rule count_peaks:
    input:
        get_macs2_peak_files(DATA_DIR, config["IGG"])
    output:
        f"{DATA_DIR}/multiqc/peakcount.txt"
    shell:
        """
        if [ -z "{input}" ]; then
            : > {output}
        else
            wc -l {input} | \
            sed '$d' | \
            sed -E 's|.*/macs2_broad_||; s|.*/macs2_narrow_||; s/_peaks.broadPeak//; s/_peaks.narrowPeak//' | \
            awk '{{print $2, $1}}' > {output}
        fi
        """
