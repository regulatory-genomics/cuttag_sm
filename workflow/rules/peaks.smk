
DATA_DIR = config["output_base_dir"].rstrip("/")

rule callpeaks:
    input:
        get_callpeaks
    output:
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/{{sample}}_peaks.bed"
    conda: 
        "../envs/gopeaks.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "gopeaks.sif")
    log:
        f"{DATA_DIR}/logs/{{sample}}_gopeaks.json"
    threads: 4
    params:
        igg = get_igg,
        params = callpeaks_params
    shell:
        f"gopeaks -b {{input[0]}} {{params.igg}} -o {DATA_DIR}/Important_processed/Peaks/callpeaks/{{wildcards.sample}} {{params.params}} > {{log}} 2>&1"

rule callpeaks_macs2_broad:
    input:
        treatment=lambda wc: get_callpeaks(wc)[0]
    output:
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_broad_{{sample}}_peaks.xls",
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_broad_{{sample}}_peaks.broadPeak",
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_broad_{{sample}}_peaks.gappedPeak"
    log:
        f"{DATA_DIR}/logs/{{sample}}_macs2peaks_broad.json"
    params:
        extra=(lambda wc: (
            f"-f BAMPE -g " + (
                "hs" if "hg" in str(config.get("GENES", "")).lower() else (
                    "mm" if "mm" in str(config.get("GENES", "")).lower() else "hs"
                )
            ) + " --nomodel --keep-dup all --broad --broad-cutoff 0.1 -q 0.01"
        ))
    threads: 4
    wrapper:
        "v2.9.1/bio/macs2/callpeak"

rule callpeaks_macs2_narrow:
    input:
        treatment=lambda wc: get_callpeaks(wc)[0]
    output:
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_narrow_{{sample}}_peaks.xls",
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_narrow_{{sample}}_peaks.narrowPeak",
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_narrow_{{sample}}_summits.bed"
    log:
        f"{DATA_DIR}/logs/{{sample}}_macs2peaks_narrow.json"
    params:
        extra=(lambda wc: (
            f"-f BAMPE -g " + (
                "hs" if "hg" in str(config.get("GENES", "")).lower() else (
                    "mm" if "mm" in str(config.get("GENES", "")).lower() else "hs"
                )
            ) + " --nomodel --keep-dup all -q 0.01"
        ))
    threads: 4
    wrapper:
        "v2.9.1/bio/macs2/callpeak"


if os.path.isfile(blacklist_file):
    rule remove_blacklist:
        input:
            f"{DATA_DIR}/Important_processed/Peaks/callpeaks/{{sample}}_peaks.bed"
        output:
            f"{DATA_DIR}/Important_processed/Peaks/callpeaks/{{sample}}_peaks_noBlacklist.bed"
        params:
            blacklist = blacklist_file
        conda:
            "../envs/bedtools.yml"
        singularity:
            "docker://staphb/bedtools:2.30.0"
        threads: 1
        shell:
            "bedtools intersect -v -a {input} -b {params.blacklist} > {output}"

    # merge all peaks to get union peak with at least
    # two reps per condition per peak
    rule make_high_conf_peaks:
        input:
            get_peaks_by_mark_condition_blacklist
        output:
            f"{DATA_DIR}/Important_processed/Peaks/highConf/{{mark_condition}}.highConf.bed"
        conda:
            "../envs/bedtools.yml"
        singularity:
            "docker://staphb/bedtools:2.30.0"
        shell:
            "files=( {input} ); existing=(); "
            "for f in ${{files[@]}}; do [ -s \"$f\" ] && existing+=(\"$f\"); done; "
            "if [ ${{#existing[@]}} -eq 0 ]; then : > {output}; else "
            "cat ${{existing[@]}} | sort -k1,1 -k2,2n | "
            "bedtools merge | "
            "bedtools intersect -a - -b ${{existing[@]}} -c | "
            "awk -v OFS='\t' '$4>=2 {{print}}' > {output}; fi"
else:
    # merge all peaks to get union peak with at least
    # two reps per condition per peak
    rule make_high_conf_peaks:
        input:
            get_peaks_by_mark_condition
        output:
            f"{DATA_DIR}/Important_processed/Peaks/highConf/{{mark_condition}}.highConf.bed"
        conda:
            "../envs/bedtools.yml"
        singularity:
            "docker://staphb/bedtools:2.30.0"
        shell:
            "files=( {input} ); existing=(); "
            "for f in ${{files[@]}}; do [ -s \"$f\" ] && existing+=(\"$f\"); done; "
            "if [ ${{#existing[@]}} -eq 0 ]; then : > {output}; else "
            "cat ${{existing[@]}} | sort -k1,1 -k2,2n | "
            "bedtools merge | "
            "bedtools intersect -a - -b ${{existing[@]}} -c | "
            "awk -v OFS='\t' '$4>=2 {{print}}' > {output}; fi"


# get consensus
rule consensus:
    input:
        expand(f"{DATA_DIR}/Important_processed/Peaks/callpeaks/{{sample}}_peaks_noBlacklist.bed", sample=sample_noigg)
        if os.path.isfile(blacklist_file)
        else expand(
            f"{DATA_DIR}/Important_processed/Peaks/callpeaks/{{sample}}_peaks.bed",
            sample=sample_noigg,
        )
    output:
        consensus_counts = f"{DATA_DIR}/Important_processed/Peaks/counts/{{mark}}_counts.tsv",
        consensus_bed = f"{DATA_DIR}/Important_processed/Peaks/counts/{{mark}}_consensus.bed"
    params:
        blacklist_flag = "-b" if os.path.isfile(blacklist_file) else ""
    conda:
        "../envs/bedtools.yml"
    singularity:
        "docker://staphb/bedtools:2.30.0"
    shell:
        f"OUTPUT_BASE_DIR=\"{DATA_DIR}\" bash workflow/src/consensus_peaks.sh -m {{wildcards.mark}} -n {{config[N_INTERSECTS]}} -o {{output.consensus_counts}} {{params.blacklist_flag}}"

rule frip:
    input:
        peaks=lambda wc: (
            f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_broad_{wc.sample}_peaks.broadPeak"
            if is_broad_mark(wc)
            else f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_narrow_{wc.sample}_peaks.narrowPeak"
        ),
        bam=rules.markdup.output.bam,
        bai=rules.index_bam.output.bai
    output:
        f"{DATA_DIR}/Report/plotEnrichment/frip_{{sample}}.png", f"{DATA_DIR}/Report/plotEnrichment/frip_{{sample}}.tsv"
    conda:
        "../envs/dtools.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "dtools.sif")
    log:
        f"{DATA_DIR}/logs/plotEnrichment_{{sample}}.log"
    shell:
        "bash workflow/src/skip_frip.sh {input.peaks} {input.bam} {output[0]} {output[1]} {log}"

rule genomic_coverage:
    input:
        peaks=lambda wc: (
            f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_broad_{wc.sample}_peaks.broadPeak"
            if is_broad_mark(wc)
            else f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_narrow_{wc.sample}_peaks.narrowPeak"
        ),
        chrom_sizes=config["CSIZES"]
    output:
        f"{DATA_DIR}/Important_processed/Peaks/coverage/{{sample}}_coverage.tsv"
    conda:
        "../envs/bedtools.yml"
    singularity:
        "docker://staphb/bedtools:2.30.0"
    log:
        f"{DATA_DIR}/logs/coverage_{{sample}}.log"
    shell:
        """
        set -euo pipefail
        exec >{log} 2>&1
        # Calculate total genome size
        TOTAL_GENOME_SIZE=$(awk -F'\t' '{{sum += $2}} END {{print sum}}' {input.chrom_sizes})
        
        # Calculate covered bases (merge overlapping intervals first)
        COVERED_BASES=$(sort -k1,1 -k2,2n {input.peaks} | \
                        bedtools merge -i stdin | \
                        awk -F'\t' '{{sum += $3 - $2}} END {{print sum}}')
        
        # Calculate percentage coverage
        COVERAGE_PERCENT=$(echo "scale=6; ($COVERED_BASES / $TOTAL_GENOME_SIZE) * 100" | bc)
        
        # Write output with header
        echo -e "Sample\tCoverage_Percent\tCovered_Bases\tTotal_Genome_Size" > {output}
        echo -e "{wildcards.sample}\t$COVERAGE_PERCENT\t$COVERED_BASES\t$TOTAL_GENOME_SIZE" >> {output}
        """

rule coverage_report:
    input:
        expand(f"{DATA_DIR}/Important_processed/Peaks/coverage/{{sample}}_coverage.tsv", sample=sample_noigg)
    output:
        f"{DATA_DIR}/Report/coverage_report.tsv"
    log:
        f"{DATA_DIR}/logs/coverage_report.log"
    shell:
        """
        set -euo pipefail
        exec >{log} 2>&1
        # Combine all individual coverage files into one report
        echo -e "Sample\tCoverage_Percent\tCovered_Bases\tTotal_Genome_Size" > {output}
        for file in {input}; do
            tail -n +2 "$file" >> {output}
        done
        """
