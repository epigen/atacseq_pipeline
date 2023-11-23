
rule multiqc:
    input:
        expand(os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"), sample=samples.keys()),
        expand(os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak"), sample=samples.keys()),
        expand(os.path.join(result_path, "hub","{sample}.bigWig"),sample=samples.keys()),
    output:
        multiqc_report = report(os.path.join(result_path,"report","multiqc_report.html"),
                                caption="../report/multiqc.rst",
                                category="{}_{}".format(config["project_name"], module_name),
                                subcategory="QC",
                                labels={
                                    "name": "MultiQC report",
                                    "type": "HTML",
                                }),
#         multiqc_report_dir = report(directory(os.path.join(result_path,"report")), # ERROR with SLURM dependencies when input or output is a directory
#                               caption="../report/multiqc.rst",
#                               htmlindex="multiqc_report.html",
#                               category="{}_{}".format(config["project_name"], module_name),
#                               subcategory="QC",
#                                labels={
#                                   "name": "MultiQC report",
#                                   "type": "HTML",
#                               }),
    params:
        #project_config = config["project_config"], # not easy to resolve -> think about modules... look up docs, maybe parameters can be passed directly into multiqc call?
        result_path = result_path,
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/multiqc.yaml",
    log:
        "logs/rules/multiqc.log"
    shell:
        """
        multiqc --force --outdir {params.result_path}/report --disable-atacseq-report {params.result_path}
        """
        
        # --disable-atacseq-report
#         --cl-config 'qualimap_config: { general_stats_coverage: [20,40,200] }' https://multiqc.info/docs/getting_started/config/#command-line-config
#         multiqc -fc {params.project_config} {params.result_path}