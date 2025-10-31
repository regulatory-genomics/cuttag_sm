
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
    container: None
    run:
        pd.options.plotting.backend = "plotly"
        dfs = []
        for i in sorted(input):
            cond_marker = "_".join(i.split("_")[1:3])
            temp_df = pd.read_csv(i, sep = "\t", usecols=["percent"]).rename(columns = {'percent': cond_marker})
            dfs.append(temp_df)
        frip_df = pd.concat(dfs, axis = 1)
        frip_df = frip_df.rename(index={0: 'inside'})
        frip_df.loc["outside"] = 100 - frip_df.loc['inside']
        fig = go.Figure(data=[
            go.Bar(name="inside_peaks", x=frip_df.columns, y=frip_df.loc['inside'], marker_color='rgb(255, 201, 57)'),
            go.Bar(name='outside_peaks', x=frip_df.columns, y=frip_df.loc['outside'], marker_color='rgb(0,39, 118)')
        ])
        fig.update_layout(barmode='stack', 
            title='Fraction of Reads in Peaks by Sample', 
            xaxis_tickfont_size=14, yaxis=dict(title='Fraction of reads in peaks', 
            titlefont_size=16, tickfont_size=14), xaxis=dict(title='Samples'))
        fig.write_html(str(output))


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
        "preseq lc_extrap -B {resources.defect_mode} -l 1000000000 -e 1000000000 -o {output} {input} > {log} 2>&1"

rule multiqc:
    input:
        expand(f"{DATA_DIR}/fastqc/{{read}}_fastqc.zip", read=reads),
        expand(f"{DATA_DIR}/fastq_screen/{{read}}_screen.txt", read=reads),
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
        f"multiqc {DATA_DIR}/ -f -c src/multiqc_conf.yml -o {DATA_DIR}/multiqc --ignore {DATA_DIR}/homer > {{log}} 2>&1"

# export different locales for singularity workaround: https://click.palletsprojects.com/en/8.1.x/unicode-support/


rule custom_report:
    input:
        multiqc_json = f"{DATA_DIR}/multiqc/multiqc_data/multiqc_data.json",
        high_conf_peaks = expand(f"{DATA_DIR}/highConf/{{mark_condition}}.highConf.bed", mark_condition=mark_conditions)
    output:
        f"{DATA_DIR}/custom_report/custom_report.html"
    params:
        rmd = "src/custom_report.Rmd",
        output_dir = f"{DATA_DIR}/custom_report",
        callpeaks_folder = f"{DATA_DIR}/callpeaks",
        high_conf_peaks_folder = f"{DATA_DIR}/highConf",
        blacklist = blacklist_file
    conda:
        "../envs/knit_rmd.yml"
    singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "knit_rmd.sif")
    shell:
        """
        Rscript -e 'rmarkdown::render(input=here::here("{params.rmd}"), output_dir=here::here("{params.output_dir}"), envir = new.env(), params=list(
        multiqc_json=here::here("{input.multiqc_json}"),
        callpeaks_folder=here::here("{params.callpeaks_folder}"),
        high_conf_peaks_folder=here::here("{params.high_conf_peaks_folder}"),
        blacklist_file=here::here("{params.blacklist}")
        ))'
        """