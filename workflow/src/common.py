# common holds pyhton function to be used in the snakefile
import pandas as pd
import os

def get_data_dir():
    return config["output_base_dir"].rstrip("/")

# map samples to fastqs
def get_samples():
    """
    return list of samples from samplesheet.tsv
    """
    return sorted(list(pd.Index(st.index).unique()))

def get_marks():
    """
    return list of marks from samplesheet.tsv
    """
    return sorted(list(st['mark'].astype(str).unique()))

def get_mark_conditions():
    """
    return list of samples by condition
    """
    st['mark_condition']=st['mark'].astype(str)+"_"+st['condition']
    return st['mark_condition'].unique().tolist()

def get_tracks_by_mark_condition(wildcards):
    """
    return list of tracks by mark_condition
    """
    st['mark_condition']=st['mark'].astype(str)+"_"+st['condition']
    samps=st.groupby(["mark_condition"])['sample'].apply(list)[wildcards.mark_condition]
    return [f"{get_data_dir()}/tracks/{s}.bw" for s in samps]
    
def get_peaks_by_mark_condition(wildcards):
    """
    return list of peaks by mark_condition
    """
    st['mark_condition']=st['mark'].astype(str)+"_"+st['condition']
    samps=st.groupby(["mark_condition"])['sample'].apply(list)[wildcards.mark_condition]
    return [f"{get_data_dir()}/callpeaks/{s}_peaks.bed" for s in samps]

def get_peaks_by_mark_condition_blacklist(wildcards):
    """
    return list of peaks by mark_condition
    """
    st['mark_condition']=st['mark'].astype(str)+"_"+st['condition']
    samps=st.groupby(["mark_condition"])['sample'].apply(list)[wildcards.mark_condition]
    return [f"{get_data_dir()}/callpeaks/{s}_peaks_noBlacklist.bed" for s in samps]


def get_bowtie2_input(wildcards):
    """
    returns reads associated with a sample
    """
    try:
        r1=st.loc[wildcards.sample]['R1']
        r2=st.loc[wildcards.sample]['R2']
        return r1,r2
    except Exception:
        row = st[st['sample'] == wildcards.sample].iloc[0]
        return row['R1'], row['R2']

def get_bowtie2_input_by_run(wildcards):
    """
    returns reads associated with a specific sample/run row
    requires 'sample' and 'run' wildcards
    """
    row = st[(st['sample'] == wildcards.sample) & (st['run'] == wildcards.run)]
    if row.empty:
        raise ValueError(f"No samplesheet row for sample={wildcards.sample}, run={wildcards.run}")
    row = row.iloc[0]
    return row['R1'], row['R2']

def get_runs_for_sample(wildcards):
    """
    return list of run identifiers for a biological sample
    """
    runs = st[st['sample'] == wildcards.sample]['run'].astype(str).tolist()
    return runs

def get_sorted_bams_for_sample(wildcards):
    """
    Build paths to per-run sorted BAMs for a sample
    """
    runs = get_runs_for_sample(wildcards)
    return [f"{get_data_dir()}/aligned/{wildcards.sample}.{r}.sort.bam" for r in runs]

def get_reads():
    """
    get list of all reads
    """
    reads = []
    for _, row in st.iterrows():
        sample = row['sample']
        run = str(row['run']) if 'run' in row else '1'
        reads.append(f"{sample}.{run}_R1")
        reads.append(f"{sample}.{run}_R2")
    return reads

def get_fastqc_outputs():
    outputs = []
    for _, row in st.iterrows():
        sample = row['sample']
        run = str(row['run']) if 'run' in row else '1'
        for mate in ['R1','R2']:
            outputs.append(f"{get_data_dir()}/fastqc/{sample}.{run}_{mate}_fastqc.zip")
    return outputs

def get_fastq_screen_outputs():
    outputs = []
    for _, row in st.iterrows():
        sample = row['sample']
        run = str(row['run']) if 'run' in row else '1'
        for mate in ['R1','R2']:
            outputs.append(f"{get_data_dir()}/fastq_screen/{sample}.{run}_{mate}_screen.txt")
    return outputs

def get_igg(wildcards):
    """
    Returns the igg file for the sample unless
    the sample is IgG then no control file is used.
    """ 
    if config['USEIGG']:
        row = st[st['sample'] == wildcards.sample]
        igg = str(row['igg'].iloc[0]) if not row.empty else ''
        iggbam=f'{get_data_dir()}/markd/{igg}.sorted.markd.bam'
        isigg=config['IGG'] in wildcards.sample
        if not isigg:
            return f'-c {iggbam}'
        else:
            return ""
    else:
        return ""

def get_callpeaks(wildcards):
    """
    Returns the callpeaks input files
    """
    bam=f"{get_data_dir()}/markd/{wildcards.sample}.sorted.markd.bam"
    # Only the BAM is needed by gopeaks; index will be created by previous rule
    return [bam]

def callpeaks_params(wildcards):
    """
    Returns callpeaks parameters specified by the user in the samplesheet
    """
    row = st[st['sample'] == wildcards.sample]
    if row.empty or 'gopeaks' not in row.columns:
        return ""
    params = str(row['gopeaks'].iloc[0])
    if params == "-" or params.lower() == "nan":
        return ""
    return params

def defect_mode(wildcards, attempt):
    if attempt == 1:
        return ""
    elif attempt > 1:
        return "-D"
