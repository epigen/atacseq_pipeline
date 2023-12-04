
# create sample annotation file based on MultiQC general stats
rule sample_annotation:
    input:
        multiqc_stats = os.path.join(result_path, "report", "multiqc_report_data", "multiqc_general_stats.txt"),
    output:
        sample_annot = os.path.join(result_path, "counts", "annotation.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    log:
        "logs/rules/sample_annotation.log"
    run:
        annot_df = pd.read_csv(input.multiqc_stats, delimiter='\t', index_col=0).loc[samples_quantify,:]
        annot_df.columns = [col.split("mqc-generalstats-")[1].replace("atac_seq_pipeline-", "").replace('-', '_') for col in annot_df.columns]
        annot_df.index.names = ['sample_name']
        annot_df.to_csv(output.sample_annot)

# generate consensus regions using (py)bedtools
rule get_consensus_regions:
    input:
        summits_bed = expand(os.path.join(result_path,"results","{sample}","peaks","{sample}_summits.bed"), sample=samples_quantify),
    output:
        consensus_regions = os.path.join(result_path,"counts","consensus_regions.bed"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/get_consensus_regions.log"
    script:
        "../scripts/get_consensus_regions.py"

# quantify coverage based on consensus regions support for every sample
rule quantify_support_sample:
    input:
        consensus_regions = os.path.join(result_path,"counts","consensus_regions.bed"),
        peakfile = os.path.join(result_path,"results","{sample}","peaks", "{sample}_summits.bed"),
    output:
        quant_support = os.path.join(result_path,"results","{sample}","peaks", "{sample}_quantification_support.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantify_support_sample_{sample}.log"
    script:
        "../scripts/quantify_support_sample.py"

# quantify coverage based on consensus regions counts for every sample
rule quantify_counts_sample:
    input:
        consensus_regions = os.path.join(result_path,"counts","consensus_regions.bed"),
        bamfile = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"),
    output:
        quant_counts = os.path.join(result_path,"results","{sample}","mapped", "{sample}_quantification_counts.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantify_counts_sample_{sample}.log"
    script:
        "../scripts/quantify_counts_sample.py"
        
# aggregate quantification of counts/support of all samples
rule quantify_aggregate:
    input:
        get_quantifications,
    output:
        os.path.join(result_path,"counts","{kind}.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 4)
    conda:
        "../envs/datamash.yaml",
    log:
        "logs/rules/quantify_aggregate_{kind}.log"
    shell:
        """
        awk 'NR==1 {{print; next}} FNR>1 {{print}}' {input} | datamash transpose -t ',' > {output}
        """