
DATA_DIR = config["output_base_dir"].rstrip("/")

rule callpeaks:
    input:
        get_callpeaks
    output:
        f"{DATA_DIR}/callpeaks/{{sample}}_peaks.bed"
    conda: 
        "../envs/gopeaks.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "gopeaks.sif")
    log:
        f"{DATA_DIR}/callpeaks/{{sample}}_gopeaks.json"
    params:
        igg = get_igg,
        params = callpeaks_params
    shell:
        f"gopeaks -b {{input[0]}} {{params.igg}} -o {DATA_DIR}/callpeaks/{{wildcards.sample}} {{params.params}} > {{log}} 2>&1"

rule callpeaks_macs2_broad:
    input:
        treatment=lambda wc: get_callpeaks(wc)[0]
    output:
        f"{DATA_DIR}/callpeaks/macs2_broad_{{sample}}_peaks.xls",
        f"{DATA_DIR}/callpeaks/macs2_broad_{{sample}}_peaks.broadPeak",
        f"{DATA_DIR}/callpeaks/macs2_broad_{{sample}}_peaks.gappedPeak"
    log:
        f"{DATA_DIR}/callpeaks/{{sample}}_macs2peaks_broad.json"
    params:
        extra=(lambda wc: (
            f"-f BAMPE -g " + (
                "hs" if "hg" in str(config.get("GENES", "")).lower() else (
                    "mm" if "mm" in str(config.get("GENES", "")).lower() else "hs"
                )
            ) + " --nomodel --keep-dup all --broad"
        ))
    threads: 4
    wrapper:
        "v2.9.1/bio/macs2/callpeak"

rule callpeaks_macs2_narrow:
    input:
        treatment=lambda wc: get_callpeaks(wc)[0]
    output:
        f"{DATA_DIR}/callpeaks/macs2_narrow_{{sample}}_peaks.xls",
        f"{DATA_DIR}/callpeaks/macs2_narrow_{{sample}}_peaks.narrowPeak",
        f"{DATA_DIR}/callpeaks/macs2_narrow_{{sample}}_summits.bed"
    log:
        f"{DATA_DIR}/callpeaks/{{sample}}_macs2peaks_narrow.json"
    params:
        extra=(lambda wc: (
            f"-f BAMPE -g " + (
                "hs" if "hg" in str(config.get("GENES", "")).lower() else (
                    "mm" if "mm" in str(config.get("GENES", "")).lower() else "hs"
                )
            ) + " --nomodel --keep-dup all"
        ))
    threads: 4
    wrapper:
        "v2.9.1/bio/macs2/callpeak"


if os.path.isfile(blacklist_file):
    rule remove_blacklist:
        input:
            f"{DATA_DIR}/callpeaks/{{sample}}_peaks.bed"
        output:
            f"{DATA_DIR}/callpeaks/{{sample}}_peaks_noBlacklist.bed"
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
            f"{DATA_DIR}/highConf/{{mark_condition}}.highConf.bed"
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
            f"{DATA_DIR}/highConf/{{mark_condition}}.highConf.bed"
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
        expand(f"{DATA_DIR}/callpeaks/{{sample}}_peaks_noBlacklist.bed", sample=sample_noigg) if os.path.isfile(blacklist_file) else expand(f"{DATA_DIR}/callpeaks/{{sample}}_peaks.bed", sample=sample_noigg)
    output:
        consensus_counts = f"{DATA_DIR}/counts/{{mark}}_counts.tsv",
        consensus_bed = f"{DATA_DIR}/counts/{{mark}}_consensus.bed"
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
        rules.callpeaks.output, f"{DATA_DIR}/markd/{{sample}}.sorted.markd.bam"
    output:
        f"{DATA_DIR}/plotEnrichment/frip_{{sample}}.png", f"{DATA_DIR}/plotEnrichment/frip_{{sample}}.tsv"
    conda:
        "../envs/dtools.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "dtools.sif")
    log:
        f"{DATA_DIR}/logs/plotEnrichment_{{sample}}.log"
    shell:
        "bash workflow/src/skip_frip.sh {input[0]} {input[1]} {output[0]} {output[1]} {log}"
