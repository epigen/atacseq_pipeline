
# generate consensus regions set using (py)bedtools
rule get_consensus_regions:
    input:
        annotation = config["sample_metadata"],
        summits_bed = expand(os.path.join(results_dir,"{sample_name}","peaks","{sample_name}_summits.bed"),sample_name=samples.keys()),
    output:
        consensus_regions = os.path.join(config["project_path"],'all',"consensus_regions.bed"),
        annotation_filtered = os.path.join(config["project_path"],'all',"all_annotation.csv"),
    params:
        # paths
        results_dir = results_dir,
        # bedtools params
        sloop_extension=250,
        # pipeline information
        blacklist_file = config["blacklisted_regions"],
        chrom_file = config["chromosome_sizes"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/get_consensus_regions.log"
    script:
        "../scripts/get_consensus_regions.py"

# quantify coverage based on consensus regions support for every sample
rule quantify_support_sample:
    input:
        consensus_regions = os.path.join(config["project_path"],'all',"consensus_regions.bed"),
        peakfile = os.path.join(results_dir,"{sample}","peaks", "{sample}_summits.bed"),
    output:
        quant_support=os.path.join(results_dir,"{sample}","peaks", "{sample}_quantification_support.csv"),
    params:
        # pipeline information
        chrom_file = config["chromosome_sizes"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/quantify_support_sample_{sample}.log"
    script:
        "../scripts/quantify_support_sample.py"

# aggregate support quantification of all samples
rule quantify_support_aggregate:
    input:
        annotation_filtered = os.path.join(config["project_path"],'all',"all_annotation.csv"),
        quant_support=expand(os.path.join(results_dir,"{sample_name}","peaks", "{sample_name}_quantification_support.csv"), sample_name=samples.keys()),
    output:
        quant_support=os.path.join(config["project_path"],'all',"all_support.csv"),
    params:
        # paths
        results_dir = results_dir,
        # cluster parameters
        partition=partition,
    threads: 4
    resources:
        mem="64G",
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/quantify_support_aggregate.log"
    script:
        "../scripts/quantify_support_aggregate.py"

# quantify coverage based on consensus regions counts for every sample
rule quantify_counts_sample:
    input:
        consensus_regions = os.path.join(config["project_path"],'all',"consensus_regions.bed"),
        bamfile = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam"),
    output:
        quant_counts=os.path.join(results_dir,"{sample}","mapped", "{sample}_quantification_counts.csv"),
    params:
        # paths
        results_dir = results_dir,
        # pipeline information
        chrom_file = config["chromosome_sizes"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/quantify_counts_sample_{sample}.log"
    script:
        "../scripts/quantify_counts_sample.py"

# aggregate count quantification of all samples (8h for >650 samples)
rule quantify_counts_aggregate:
    input:
        annotation_filtered = os.path.join(config["project_path"],'all',"all_annotation.csv"),
        quant_counts=expand(os.path.join(results_dir,"{sample_name}","mapped", "{sample_name}_quantification_counts.csv"), sample_name=samples.keys()),
    output:
        quant_counts=os.path.join(config["project_path"],'all',"all_counts.csv"),
    params:
        # paths
        results_dir = results_dir,
        # cluster parameters
        partition=partition,
    threads: 4
    resources:
        mem="64G",
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/quantify_counts_aggregate.log"
    script:
        "../scripts/quantify_counts_aggregate.py"