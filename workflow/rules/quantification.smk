
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

# generate promoter regions using (py)bedtools
rule get_promoter_regions:
    input:
        config["gencode_gtf"],
    output:
        promoter_regions = os.path.join(result_path,"counts","promoter_regions.bed"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb = config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/get_promoter_regions.log"
    script:
        "../scripts/get_promoter_regions.py"
        
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
        quant_support = os.path.join(result_path,"results","{sample}","peaks", "{sample}_quantification_support_counts.csv"),
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
        regions = os.path.join(result_path,"counts","{kind}_regions.bed"),
        bamfile = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"),
    output:
        quant_counts = os.path.join(result_path,"results","{sample}","mapped", "{sample}_quantification_{kind}_counts.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/quantify_sample_{sample}_{kind}.log"
    script:
        "../scripts/quantify_counts_sample.py"
        
# aggregate quantification of counts/support of all samples
rule quantify_aggregate:
    input:
        get_quantifications,
    output:
        os.path.join(result_path,"counts","{kind}_counts.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: 2*config.get("threads", 2)
    conda:
        "../envs/datamash.yaml",
    log:
        "logs/rules/quantify_aggregate_{kind}.log"
    shell:
        """
        awk 'NR==1 {{print; next}} FNR>1 {{print}}' {input} | datamash transpose -t ',' > {output}
        """
        
# aggregate HOMER motif enrichment results for all QC'd samples into one CSV
rule homer_aggregate:
    input:
        expand(os.path.join(result_path,"results","{sample}","homer","knownResults.txt"), sample=samples_quantify),
    output:
        os.path.join(result_path,"counts","HOMER_knownMotifs.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 2)
    log:
        "logs/rules/homer_aggregate.log"
    shell:
        """
        first_file=true
        
        > {output}
        
        for file in {input}; do
            sample=$(awk -F'/' '{{print $(NF-2)}}' <<< "$file")

            if $first_file; then
                awk -v prefix="$sample" 'BEGIN {{FS="\\t"; OFS=","}} {{
                    if (NR==1) {{
                        gsub(" ", "_", $0);
                        printf "sample" OFS; 
                        for (i=1; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS);
                    }} else {{
                        printf prefix OFS "\\"" $1 "\\"";
                        for (i=2; i<=NF; i++) printf "%s%s", OFS, $i;
                        print "";
                    }}
                }}' "$file" >> {output}
                
                first_file=false
            else
                awk -v prefix="$sample" 'BEGIN {{FS="\\t"; OFS=","}} FNR > 1 {{printf prefix OFS "\\"" $1 "\\"";
                    for (i=2; i<=NF; i++) printf "%s%s", OFS, $i;
                    print "";
                    }}' "$file" >> {output}
            fi
        done
        """

# map consensus regions to closest TSS per gene
rule map_consensus_tss:
    input:
        region_annotation = os.path.join(result_path,'counts',"region_annotation.csv"),
        consensus_counts = os.path.join(result_path,"counts","consensus_counts.csv"),
    output:
        tss_counts = os.path.join(result_path,"counts","TSS_counts.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/map_consensus_tss.log"
    script:
        "../scripts/map_consensus_tss.py"