rule multiqc:
    input:
        expand(os.path.join(results_dir,"{sample_name}","mapped", "{sample_name}.filtered.bam"), sample_name=samples.keys()),
        expand(os.path.join(results_dir,"{sample_name}","peaks","{sample_name}_peaks.narrowPeak"), sample_name=samples.keys()),
        expand(os.path.join(config["atacseq.project_path"], "atacseq_hub","{sample_name}.bigWig"),sample_name=samples.keys()),
    output:
        multiqc_report=report(directory(os.path.join(config["atacseq.project_path"],"atacseq_report")), caption="../report/multiqc.rst", htmlindex="multiqc_report.html", category="QC reports"),  # for inclusion into snakemake report
    params:
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/multiqc.yaml",
    log:
        "logs/rules/multiqc.log"
    shell:
        """
        multiqc -fc {config[atacseq.project_config]} {config[atacseq.project_path]}
        """