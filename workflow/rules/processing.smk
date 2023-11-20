# alignment with botwtie2 & samtools
rule align:
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
        "../envs/bowtie2.yaml",
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

rule coverage:
    input:
        bam = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam"),
        bai = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam.bai"),
    output:
        bigWig = os.path.join(config["project_path"], "atacseq_hub","{sample}.bigWig"),
        tss_hist = os.path.join(results_dir,"{sample}","{sample}.tss_histogram.csv"),
    params:
        # paths
        hub_dir = os.path.join(config["project_path"], "atacseq_hub"),
        sample_dir = os.path.join(results_dir,"{sample}"),
        # sample information
        sample_name= lambda w: samples["{}".format(w.sample)]["sample_name"],
        genome_size= lambda w: samples["{}".format(w.sample)]["genome_size"],
        # misc task parameters
        tss_slop = tss_slop,
        noise_lower = noise_lower,
        noise_upper = noise_upper,
        double_slop = double_slop,
        # pipeline information
        unique_tss = config["unique_tss"],
        chromosome_sizes = config["chromosome_sizes"],
        whitelist = config["whitelisted_regions"],
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/misc_tasks_{sample}.log"
    shell:
        """
        mkdir -p {params.hub_dir};

        bamCoverage --bam {input.bam} \
            -p max --binSize 10  --normalizeUsing RPGC \
            --effectiveGenomeSize {params.genome_size} --extendReads 175 \
            -o "{params.hub_dir}/{params.sample_name}.bigWig" > "{params.hub_dir}/{params.sample_name}.bigWig.log" 2>&1;

        echo "base,count" > {output.tss_hist};
        bedtools slop -b {params.tss_slop} -i {params.unique_tss} -g {params.chromosome_sizes} | \
            bedtools coverage -a - -b {input.bam} -d -sorted | \
            awk '{{if($6 == "+"){{ counts[$7] += $8;}} else counts[{params.double_slop} - $7 + 1] += $8;}} END {{ for(pos in counts) {{ if(pos < {params.noise_lower} || pos > {params.noise_upper}) {{ noise += counts[pos] }} }}; average_noise = noise /(2 * {params.noise_lower}); for(pos in counts) {{print pos-2000-1","(counts[pos]/average_noise) }} }}' | \
            sort -t "," -k1,1n >> {output.tss_hist} ;
        """
        
# peak calling with macs2 & samtools
rule peak_calling:
    input:
        bam = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam"),
        bai = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam.bai"),
        #homer_script=os.path.join("resources","homer","configureHomer.pl"), # needed so that homer installation rule is executed beforehand
    output:
        peak_calls = os.path.join(results_dir,"{sample}","peaks","{sample}_peaks.narrowPeak"),
        macs2_xls = os.path.join(results_dir,"{sample}","peaks","{sample}_peaks.xls"),
        summits_bed = os.path.join(results_dir,"{sample}","peaks","{sample}_summits.bed"),
        homer_log = os.path.join(results_dir,"{sample}","homer","{sample}.homer.log"),
    params:
        # paths
        results_dir = results_dir,
        sample_dir = os.path.join(results_dir,"{sample}"),
        peaks_dir = os.path.join(results_dir,"{sample}","peaks"),
        homer_dir = os.path.join(results_dir,"{sample}","homer"),
        homer_bin = os.path.join(os.getcwd(),"resources","homer","bin"),
        # sample information
        sample_name= lambda w: samples["{}".format(w.sample)]["sample_name"],
        genome= lambda w: samples["{}".format(w.sample)]["genome"],
        genome_size= lambda w: samples["{}".format(w.sample)]["genome_size"],
        # peak calling parameters
        formating = lambda w: '--format BAMPE' if samples["{}".format(w.sample)]["read_type"] == "paired" else '--format BAM',
        # pipeline information
        regulatory_regions = config["regulatory_regions"],
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml",
    log:
        "logs/rules/peak_calling_{sample}.log"
    shell:
        """
        mkdir -p {params.results_dir};
        mkdir -p {params.sample_dir};
        mkdir -p {params.peaks_dir};
        mkdir -p {params.homer_dir};
        
        export PATH="{params.homer_bin}:$PATH";
        
        macs2 callpeak -t {input.bam} {params.formating} \
            --nomodel --keep-dup auto --extsize 147 -g {params.genome_size} \
            -n {params.sample_name} \
            --outdir {params.peaks_dir} > "{params.peaks_dir}/{params.sample_name}.macs2.log" 2>&1;
        
        resources/homer/bin/annotatePeaks.pl {params.peaks_dir}/{params.sample_name}_peaks.narrowPeak {params.genome} \
            > {params.peaks_dir}/{params.sample_name}_peaks.narrowPeak.annotated.tsv \
            2> {params.peaks_dir}/{params.sample_name}_peaks.narrowPeak.annotated.tsv.log;
        
        resources/homer/bin/findMotifsGenome.pl "{params.peaks_dir}/{params.sample_name}_summits.bed" {params.genome} {params.homer_dir} -size 200 -mask \
            > "{params.homer_dir}/{params.sample_name}.homer.log" 2>&1

        cat {params.peaks_dir}/{params.sample_name}_peaks.narrowPeak | wc -l | \
            awk '{{print "peaks\t" $1}}' >> "{params.sample_dir}/{params.sample_name}.stats.tsv"
        
        TOTAL_READS=`samtools idxstats {input.bam} | awk '{{sum += $3}}END{{print sum}}'`;
        samtools view -c -L "{params.peaks_dir}/{params.sample_name}_peaks.narrowPeak" {input.bam} | \
            awk -v total=$TOTAL_READS '{{print "frip\t" $1/total}}' >> "{params.sample_dir}/{params.sample_name}.stats.tsv";

        samtools view -c -L {params.regulatory_regions} {input.bam} | \
            awk -v total=$TOTAL_READS '{{print "regulatory_fraction\t" $1/total}}' >> "{params.sample_dir}/{params.sample_name}.stats.tsv";
        """ 