# alignment with botwtie2 & samtools
rule align:
    input:
        get_raw_bams,
    output:
        bam = os.path.join(result_path,"results","{sample}","mapped", "{sample}.bam"),
        output_bai =  os.path.join(result_path,"results","{sample}","mapped", "{sample}.bam.bai"),
        filtered_bam = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"),
        filtered_bai = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam.bai"),
        bowtie_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.txt'),
        bowtie_met = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.bowtie2.met'),
        fastp_html = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.fastp.html'),
        fastp_json = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.fastp.json'),
        fastp_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.fastp.log'),
        samblaster_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.samblaster.log'),
        samtools_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.samtools.log'),
        samtools_flagstat_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.samtools_flagstat.log'),
        stats = os.path.join(result_path, 'results', "{sample}", '{sample}.align.stats.tsv'),
    params:
        interleaved_in = lambda w: "--interleaved_in" if samples["{}".format(w.sample)]["read_type"] == "paired"  else " ",
        interleaved = lambda w: "--interleaved" if samples["{}".format(w.sample)]["read_type"] == "paired" else " ",
        filtering = lambda w: "-q 30 -F {flag} -f 2 -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]) if samples["{}".format(w.sample)]["read_type"] == "paired" else "-q 30 -F {flag} -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]),
        add_mate_tags = lambda w: "--addMateTags" if samples["{}".format(w.sample)]["read_type"] == "paired" else " ",
        adapter_sequence = "-a " + config["adapter_sequence"] if config["adapter_sequence"] !="" else " ",
        adapter_fasta = "--adapter_fasta " + config["adapter_fasta"] if config["adapter_fasta"] !="" else " ",
        sequencing_platform = config["sequencing_platform"],
        sequencing_center = config["sequencing_center"],
        mitochondria_name = config["mitochondria_name"],
        bowtie2_index = config["bowtie2_index"],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 4*config.get("threads", 2)
    conda:
        "../envs/bowtie2.yaml",
    log:
        "logs/rules/align_{sample}.log"
    shell:
        """
        result_path=$(dirname {output.stats})
        find $result_path -type f -name '*.bam.tmp.*' -exec rm {{}} +;
        
        RG="--rg-id {wildcards.sample} --rg SM:{wildcards.sample} --rg PL:{params.sequencing_platform} --rg CN:{params.sequencing_center}"

        for i in {input}; do samtools fastq $i 2>> "{output.samtools_log}" ; done | \
            fastp {params.adapter_sequence} {params.adapter_fasta} --stdin {params.interleaved_in} --stdout --html "{output.fastp_html}" --json "{output.fastp_json}" 2> "{output.fastp_log}" | \
            bowtie2 $RG --very-sensitive --no-discordant -p {threads} --maxins 2000 -x {params.bowtie2_index} --met-file "{output.bowtie_met}" {params.interleaved} - 2> "{output.bowtie_log}" | \
            samblaster {params.add_mate_tags} 2> "{output.samblaster_log}" | \
            samtools sort -o "{output.bam}" - 2>> "{output.samtools_log}";
        
        samtools index "{output.bam}" 2>> "{output.samtools_log}";
        samtools idxstats "{output.bam}" | awk '{{ sum += $3 + $4; if($1 == "{params.mitochondria_name}") {{ mito_count = $3; }}}}END{{ print "mitochondrial_fraction\t"mito_count/sum }}' > "{output.stats}";
        samtools flagstat "{output.bam}" > "{output.samtools_flagstat_log}";

        samtools view {params.filtering} -o "{output.filtered_bam}" "{output.bam}";
        samtools index "{output.filtered_bam}";
        """

rule tss_coverage:
    input:
        bam = os.path.join(result_path,"results","{sample}","mapped","{sample}.filtered.bam"),
        bai = os.path.join(result_path,"results","{sample}","mapped","{sample}.filtered.bam.bai"),
    output:
        tss_hist = os.path.join(result_path,"results","{sample}","{sample}.tss_histogram.csv"),
    params:
        noise_upper = ( config["tss_slop"] * 2 ) - config["noise_lower"],
        double_slop = ( config["tss_slop"] * 2 ),
        genome_size = config["genome_size"],
        tss_slop = config["tss_slop"],
        unique_tss = config["unique_tss"],
        chromosome_sizes = config["chromosome_sizes"],
        noise_lower = config["noise_lower"],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 4*config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/coverage_{sample}.log"
    shell:
        """
        echo "base,count" > {output.tss_hist};
        bedtools slop -b {params.tss_slop} -i {params.unique_tss} -g {params.chromosome_sizes} | \
            bedtools coverage -a - -b {input.bam} -d -sorted | \
            awk '{{if($6 == "+"){{ counts[$7] += $8;}} else counts[{params.double_slop} - $7 + 1] += $8;}} END {{ for(pos in counts) {{ if(pos < {params.noise_lower} || pos > {params.noise_upper}) {{ noise += counts[pos] }} }}; average_noise = noise /(2 * {params.noise_lower}); for(pos in counts) {{print pos-2000-1","(counts[pos]/average_noise) }} }}' | \
            sort -t "," -k1,1n >> {output.tss_hist} ;
        """
        
# peak calling with MACS2 & samtools and annotation with HOMER
rule peak_calling:
    input:
        bam = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"),
        bai = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam.bai"),
        homer_script = os.path.join(HOMER_path,"configureHomer.pl"),
    output:
        peak_calls = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak"),
        peak_annot = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak.annotated.tsv"),
        peak_annot_log = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak.annotated.tsv.log"),
        macs2_xls = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.xls"),
        summits_bed = os.path.join(result_path,"results","{sample}","peaks","{sample}_summits.bed"),
        homer_knownResults = os.path.join(result_path,"results","{sample}","homer","knownResults.txt"),
        homer_log = os.path.join(result_path,"results","{sample}","homer","{sample}.homer.log"),
        macs2_log = os.path.join(result_path, 'results', "{sample}", 'peaks', '{sample}.macs2.log'),
        stats = os.path.join(result_path, 'results', "{sample}", '{sample}.peak.stats.tsv'),
    params:
        peaks_dir = os.path.join(result_path,"results","{sample}","peaks"),
        homer_dir = os.path.join(result_path,"results","{sample}","homer"),
        homer_bin = os.path.join(HOMER_path,"bin"),
        formating = lambda w: '--format BAMPE' if samples["{}".format(w.sample)]["read_type"] == "paired" else '--format BAM',
        genome_size = config["genome_size"],
        genome = config["genome"],
        regulatory_regions = config["regulatory_regions"],
        keep_dup = config['macs2_keep_dup'],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 4*config.get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml",
    log:
        "logs/rules/peak_calling_{sample}.log"
    shell:
        """
        export PATH="{params.homer_bin}:$PATH";
        
        macs2 callpeak -t {input.bam} {params.formating} \
            --nomodel --keep-dup {params.keep_dup} --extsize 147 -g {params.genome_size} \
            -n {wildcards.sample} \
            --outdir {params.peaks_dir} > "{output.macs2_log}" 2>&1;
        
        {params.homer_bin}/annotatePeaks.pl {output.peak_calls} {params.genome} > {output.peak_annot} 2> {output.peak_annot_log};
        
        {params.homer_bin}/findMotifsGenome.pl "{output.summits_bed}" {params.genome} {params.homer_dir} -size 200 -mask > "{output.homer_log}" 2>&1

        cat {output.peak_calls} | wc -l | awk '{{print "peaks\t" $1}}' >> "{output.stats}"
        
        TOTAL_READS=`samtools idxstats {input.bam} | awk '{{sum += $3}}END{{print sum}}'`;
        
        samtools view -c -L {output.peak_calls} {input.bam} | awk -v total=$TOTAL_READS '{{print "frip\t" $1/total}}' >> "{output.stats}";

        samtools view -c -L {params.regulatory_regions} {input.bam} | awk -v total=$TOTAL_READS '{{print "regulatory_fraction\t" $1/total}}' >> "{output.stats}";
        
        if [ ! -f {output.homer_knownResults} ]; then
            touch {output.homer_knownResults}
        fi
        """
        
rule aggregate_stats:
    input:
        align_stats = os.path.join(result_path, 'results', "{sample}", '{sample}.align.stats.tsv'),
        peak_stats = os.path.join(result_path, 'results', "{sample}", '{sample}.peak.stats.tsv'),
    output:
        os.path.join(result_path, 'results', "{sample}", '{sample}.stats.tsv'),
    resources:
        mem_mb=config.get("mem", "1000"),
    threads: config.get("threads", 2)
    log:
        "logs/rules/aggregate_stats_{sample}.log"
    shell:
        """
        cat {input.align_stats} {input.peak_stats} > {output}
        """
