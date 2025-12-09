
DATA_DIR = config["output_base_dir"].rstrip("/")
rule plotFinger:
    input:
        rules.markdup.output.bam,
        rules.index_bam.output.bai
    output:
        f"{DATA_DIR}/Report/dtools/fingerprint_{{sample}}.tsv"
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
        expand(f"{DATA_DIR}/Report/plotEnrichment/frip_{{sample}}.tsv", sample = sample_noigg)
    output:
        f"{DATA_DIR}/Report/plotEnrichment/frip.html"
    conda:
        "../envs/plot_report.yml"
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
        f"{DATA_DIR}/Report/preseq/lcextrap_{{sample}}.txt"
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
        expand(f"{DATA_DIR}/Report/plotEnrichment/frip_{{sample}}.tsv", sample=sample_noigg),
        expand(f"{DATA_DIR}/Report/preseq/lcextrap_{{sample}}.txt", sample=samps)
    output:
        f"{DATA_DIR}/Report/multiqc/multiqc_report.html",
        f"{DATA_DIR}/Report/multiqc/multiqc_data/multiqc_data.json"
    conda:
        "../envs/multiqc.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "multiqc.sif")
    log:
        f"{DATA_DIR}/logs/multiqc.log"
    shell:
        """
        export LC_ALL=C.UTF-8; export LANG=C.UTF-8;
        
        multiqc {DATA_DIR}/Report/ {DATA_DIR}/logs/ \
            --ignore {DATA_DIR}/Report/multiqc \
            -o {DATA_DIR}/Report/multiqc \
            -f \
            -c workflow/src/multiqc_conf.yml \
            --cl-config "annotation: {SAMPLE_SHEET}" \
            >> {log} 2>&1
        """


rule count_peaks:
    input:
        get_macs2_peak_files(DATA_DIR, config["IGG"])
    output:
        f"{DATA_DIR}/Report/peak_stat/peakcount.txt"
    log:
        f"{DATA_DIR}/logs/count_peaks.log"
    shell:
        """
        echo "Input files: {input}" > {log}
        echo "Number of input files: $(echo '{input}' | wc -w)" >> {log}
        
        if [ -z "{input}" ]; then
            echo "No input files, creating empty output" >> {log}
            : > {output}
        else
            echo "Processing input files" >> {log}
            wc -l {input} \
            | awk '$2 != "total" {{file=$2; gsub(\".*/macs2_(broad|narrow)_\",\"\",file); gsub(\"_peaks\\\\.(broadPeak|narrowPeak)\",\"\",file); print file, $1}}' \
            > {output}
        fi
        """
