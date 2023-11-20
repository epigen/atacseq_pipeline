# alignment with botwtie2 & samtools
rule bowtie2_align:
    input:
        get_raw_bams,
    output:
        output_bam = os.path.join(results_dir,"{sample}","mapped", "{sample}.bam"),
        output_bai =  os.path.join(results_dir,"{sample}","mapped", "{sample}.bam.bai"),
        filtered_bam = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam"),
        filtered_bai = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam.bai"),
    params:
        # paths
        results_dir = results_dir,
        sample_dir = os.path.join(results_dir,"{sample}"),
        bam_dir = os.path.join(results_dir,"{sample}","mapped"),
        # sample information
        sample_name= lambda w: samples["{}".format(w.sample)]["sample_name"],
        read_type= lambda w: samples["{}".format(w.sample)]["read_type"],
        raw_bams= lambda w: samples["{}".format(w.sample)]["raw_bams"],
        genome= lambda w: samples["{}".format(w.sample)]["genome"],
        genome_size= lambda w: samples["{}".format(w.sample)]["genome_size"],
        #alignment parameters
        interleaved_in = lambda w: "--interleaved_in" if samples["{}".format(w.sample)]["read_type"] == "paired"  else " ",
        interleaved = lambda w: "--interleaved" if samples["{}".format(w.sample)]["read_type"] == "paired" else " ",
        filtering = lambda w: "-q 30 -F 2316 -f 2 -L {}".format(config["whitelisted_regions"]) if samples["{}".format(w.sample)]["read_type"] == "paired" else "-q 30 -F 2316 -L {}".format(config["whitelisted_regions"]),
        add_mate_tags = lambda w: "--addMateTags" if samples["{}".format(w.sample)]["read_type"] == "paired" else " ",
        # pipeline information
        adapter_sequence="-a " + config["adapter_sequence"] if config["adapter_sequence"] !="" else " ",
        adapter_fasta="--adapter_fasta " + config["adapter_fasta"] if config["adapter_fasta"] !="" else " ",
        bowtie2_index=config["bowtie2_index"],
        chrM=config["mitochondria_name"],
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/atacseq_pipeline.yaml",
    log:
        "logs/rules/bowtie2_align_{sample}.log"
    shell:
        """
        mkdir -p {params.results_dir};
        mkdir -p {params.sample_dir};
        mkdir -p {params.bam_dir};
        
        RG="--rg-id {params.sample_name} --rg SM:{params.sample_name} --rg PL:illumina --rg CN:CeMM_BSF"

        for i in {params.raw_bams}; do samtools fastq $i 2>> "{params.bam_dir}/{params.sample_name}.samtools.log" ; done | \
            fastp {params.adapter_sequence} {params.adapter_fasta} --stdin {params.interleaved_in} --stdout --html "{params.bam_dir}/{params.sample_name}.fastp.html" --json "{params.bam_dir}/{params.sample_name}.fastp.json" 2> "{params.bam_dir}/{params.sample_name}.fastp.log" | \
            bowtie2 $RG --very-sensitive --no-discordant -p {threads} --maxins 2000 -x {params.bowtie2_index} --met-file "{params.bam_dir}/{params.sample_name}.bowtie2.met" {params.interleaved} - 2> "{params.bam_dir}/{params.sample_name}.txt" | \
            samblaster {params.add_mate_tags} 2> "{params.bam_dir}/{params.sample_name}.samblaster.log" | \
            samtools sort -o "{params.bam_dir}/{params.sample_name}.bam" - 2>> "{params.bam_dir}/{params.sample_name}.samtools.log";
        
        samtools index "{params.bam_dir}/{params.sample_name}.bam" 2>> "{params.bam_dir}/{params.sample_name}.samtools.log";
        samtools idxstats "{params.bam_dir}/{params.sample_name}.bam" | awk '{{ sum += $3 + $4; if($1 == "{params.chrM}") {{ mito_count = $3; }}}}END{{ print "mitochondrial_fraction\t"mito_count/sum }}' > "{params.sample_dir}/{params.sample_name}.stats.tsv";
        samtools flagstat "{params.bam_dir}/{params.sample_name}.bam" > "{params.bam_dir}/{params.sample_name}.samtools_flagstat.log";

        samtools view {params.filtering} -o "{params.bam_dir}/{params.sample_name}.filtered.bam" "{params.bam_dir}/{params.sample_name}.bam";
        samtools index "{params.bam_dir}/{params.sample_name}.filtered.bam";
        """
