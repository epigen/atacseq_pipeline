rule misc_tasks:
    input:
        bam = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam"),
        bai = os.path.join(results_dir,"{sample}","mapped", "{sample}.filtered.bam.bai"), # TODO: not used
    output:
        bigWig = os.path.join(config["atacseq.project_path"], "atacseq_hub","{sample}.bigWig"),
        tss_hist = os.path.join(results_dir,"{sample}","{sample}.tss_histogram.csv"),
    params:
        # paths
        hub_dir = os.path.join(config["atacseq.project_path"], "atacseq_hub"),
        sample_dir = os.path.join(results_dir,"{sample}"),
        #slopped_tss = os.path.join(results_dir,"{sample}","slopped_tss.bed"), # TODO: not needed? ask BE
        # sample information
        sample_name= lambda w: samples["{}".format(w.sample)]["sample_name"],
        genome_size= lambda w: samples["{}".format(w.sample)]["genome_size"],
        # misc task parameters
        tss_slop = tss_slop,
        noise_lower = noise_lower,
        noise_upper = noise_upper,
        double_slop = double_slop,
        # pipeline information
        unique_tss = config["atacseq.unique_tss"],
        chromosome_sizes = config["atacseq.chromosome_sizes"],
        whitelist = config["atacseq.whitelisted_regions"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "envs/atacseq_pipeline.yaml",
    log:
        "results/logs/rules/misc_tasks_{sample}.log"
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