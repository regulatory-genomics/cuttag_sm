
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
        
        # Install local multiqc_cuttag plugin
        # Use --no-build-isolation to use conda's setuptools and avoid proxy issues
        pip install -e workflow/src/multiqc_cuttag --force-reinstall --no-deps --no-build-isolation >> {log} 2>&1 || echo "Warning: Plugin installation failed, continuing anyway..." >> {log} 2>&1
        
        # Run MultiQC
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
        f"{DATA_DIR}/Report/multiqc/peakcount.txt"
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
