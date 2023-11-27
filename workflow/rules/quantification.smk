
# generate consensus regions using (py)bedtools
rule get_consensus_regions:
    input:
        summits_bed = expand(os.path.join(result_path,"results","{sample}","peaks","{sample}_summits.bed"), sample=samples_quantify),
    output:
        consensus_regions = os.path.join(result_path,"counts","consensus_regions.bed"),
#         annotation_filtered = os.path.join(result_path,"counts","all_annotation.csv"),
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

# aggregate support quantification of all samples
rule quantify_support_aggregate:
    input:
        quant_support = expand(os.path.join(result_path,"results","{sample}","peaks", "{sample}_quantification_support.csv"), sample=samples_quantify),
    output:
        quant_support_agg = os.path.join(result_path,"counts","support.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 4)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantify_support_aggregate.log"
    script:
        "../scripts/aggregate_counts.py"

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

# aggregate count quantification of all samples (8h for >650 samples)
rule quantify_counts_aggregate:
    input:
        quant_counts = expand(os.path.join(result_path,"results","{sample}","mapped", "{sample}_quantification_counts.csv"), sample=samples_quantify),
    output:
        quant_counts_agg = os.path.join(result_path,"counts","counts.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 4)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantify_counts_aggregate.log"
    script:
        "../scripts/aggregate_counts.py"