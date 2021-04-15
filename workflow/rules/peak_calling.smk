# peak calling with macs2 & samtools
rule peak_calling:
    input:
        bam = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam"),
        bai = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam.bai"),
        homer_script=os.path.join("resources","homer","configureHomer.pl"), # needed so that homer installation rule is executed beforehand
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
        regulatory_regions = config["atacseq.regulatory_regions"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_pipeline.yaml",
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