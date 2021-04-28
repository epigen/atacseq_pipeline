
# generate consensus regions set using (py)bedtools
rule get_consensus_regions:
    input:
        annotation = config["atacseq.annotation_metadata"],
    output:
        consensus_regions = os.path.join(config["atacseq.project_path"],"consensus_regions.bed"),
    params:
        # paths
        results_dir = results_dir,
        # bedtools params
        sloop_extension=250,
        # pipeline information
        blacklist_file = config["atacseq.blacklisted_regions"],
        chrom_file = config["atacseq.chromosome_sizes"],
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
        
# # quantify coverage based on consensus regions binary
# rule quantify_binary:
#     input:
#         annotation = config["atacseq.annotation_metadata"],
#         consensus_regions = os.path.join(config["atacseq.project_path"],"consensus_regions.bed"),
#     output:
#         quant_binary=os.path.join(config["atacseq.project_path"],"quantification_binary.csv"),
#     params:
#         # paths
#         results_dir = results_dir,
#         # pipeline information
#         chrom_file = config["atacseq.chromosome_sizes"],
#         # cluster parameters
#         partition=partition,
#     threads: threads
#     resources:
#         mem=mem,
#     conda:
#         "../envs/atacseq_analysis.yaml",
#     log:
#         "logs/rules/quantify_binary.log"
#     script:
#         "../scripts/quantify_binary_all.py"

# quantify coverage based on consensus regions binary for every sample
rule quantify_binary_sample:
    input:
        consensus_regions = os.path.join(config["atacseq.project_path"],"consensus_regions.bed"),
        peakfile = os.path.join(results_dir,"{sample}","peaks", "{sample}_summits.bed"),
    output:
        quant_binary=os.path.join(results_dir,"{sample}","peaks", "{sample}_quantification_binary.csv"),
    params:
        # pipeline information
        chrom_file = config["atacseq.chromosome_sizes"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/quantify_binary_sample_{sample}.log"
    script:
        "../scripts/quantify_binary_sample.py"

# aggregate binary quantification of all samples
rule quantify_binary_aggregate:
    input:
        annotation = config["atacseq.annotation_metadata"],
        quant_binary=expand(os.path.join(results_dir,"{sample_name}","peaks", "{sample_name}_quantification_binary.csv"), sample_name=samples.keys()),
    output:
        quant_binary=os.path.join(config["atacseq.project_path"],"quantification_binary.csv"),
    params:
        # paths
        results_dir = results_dir,
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem="64G",
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/quantify_binary_aggregate.log"
    script:
        "../scripts/quantify_binary_aggregate.py"

# # quantify coverage based on consensus regions counts
# rule quantify_counts:
#     input:
#         annotation = config["atacseq.annotation_metadata"],
#         consensus_regions = os.path.join(config["atacseq.project_path"],"consensus_regions.bed"),
#     output:
#         quant_counts=os.path.join(config["atacseq.project_path"],"quantification_counts.csv"),
#     params:
#         # paths
#         results_dir = results_dir,
#         # pipeline information
#         chrom_file = config["atacseq.chromosome_sizes"],
#         # cluster parameters
#         partition=partition,
#     threads: threads
#     resources:
#         mem=mem,
#     conda:
#         "../envs/atacseq_analysis.yaml",
#     log:
#         "logs/rules/quantify_counts.log"
#     script:
#         "../scripts/quantify_counts.py"

# quantify coverage based on consensus regions counts for every sample
rule quantify_counts_sample:
    input:
        consensus_regions = os.path.join(config["atacseq.project_path"],"consensus_regions.bed"),
        bamfile = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam"),
    output:
        quant_counts=os.path.join(results_dir,"{sample}","mapped", "{sample}_quantification_counts.csv"),
    params:
        # paths
        results_dir = results_dir,
        # pipeline information
        chrom_file = config["atacseq.chromosome_sizes"],
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
        annotation = config["atacseq.annotation_metadata"],
        quant_counts=expand(os.path.join(results_dir,"{sample_name}","mapped", "{sample_name}_quantification_counts.csv"), sample_name=samples.keys()),
    output:
        quant_counts=os.path.join(config["atacseq.project_path"],"quantification_counts.csv"),
    params:
        # paths
        results_dir = results_dir,
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem="64G",
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/quantify_counts_aggregate.log"
    script:
        "../scripts/quantify_counts_aggregate.py"