
rule symlink_stats:
    input:
        stats_tsv = os.path.join(result_path, 'results', "{sample}", '{sample}.stats.tsv'),
        tss_csv = os.path.join(result_path, 'results', "{sample}", '{sample}.tss_histogram.csv'),
        mapped_txt = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.txt'),
        fastp_json = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.fastp.json'),
        samblaster_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.samblaster.log'),
        flagstat_log = os.path.join(result_path, 'results', "{sample}", 'mapped', '{sample}.samtools_flagstat.log'),
        macs2_log = os.path.join(result_path, 'results', "{sample}", 'peaks', '{sample}.macs2.log'),
        peaks_xls = os.path.join(result_path, 'results', "{sample}", 'peaks', '{sample}_peaks.xls'),
    output:
        stats_tsv = os.path.join(result_path, 'report', '{sample}.stats.tsv'),
        tss_csv = os.path.join(result_path, 'report', '{sample}_TSS.csv'),
        mapped_txt = os.path.join(result_path, 'report', '{sample}.txt'),
        fastp_json = os.path.join(result_path, 'report', '{sample}.fastp.json'),
        samblaster_log = os.path.join(result_path, 'report', '{sample}.samblaster.log'),
        flagstat_log = os.path.join(result_path, 'report', '{sample}.samtools_flagstat.log'),
        macs2_log = os.path.join(result_path, 'report', '{sample}.macs2.log'),
        peaks_xls = os.path.join(result_path, 'report', '{sample}_peaks.xls'),
    resources:
        mem_mb=config.get("mem", "1000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs", "rules", "symlink_stats_{sample}.log")
    shell:
        """
        ln -sfn $(realpath --relative-to=$(dirname {output.stats_tsv}) {input.stats_tsv}) {output.stats_tsv}
        ln -sfn $(realpath --relative-to=$(dirname {output.tss_csv}) {input.tss_csv}) {output.tss_csv}
        ln -sfn $(realpath --relative-to=$(dirname {output.mapped_txt}) {input.mapped_txt}) {output.mapped_txt}
        ln -sfn $(realpath --relative-to=$(dirname {output.fastp_json}) {input.fastp_json}) {output.fastp_json}
        ln -sfn $(realpath --relative-to=$(dirname {output.samblaster_log}) {input.samblaster_log}) {output.samblaster_log}
        ln -sfn $(realpath --relative-to=$(dirname {output.flagstat_log}) {input.flagstat_log}) {output.flagstat_log}
        ln -sfn $(realpath --relative-to=$(dirname {output.macs2_log}) {input.macs2_log}) {output.macs2_log}
        ln -sfn $(realpath --relative-to=$(dirname {output.peaks_xls}) {input.peaks_xls}) {output.peaks_xls}
        """
            
rule multiqc:
    input:
        expand(os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"), sample=samples.keys()),
        expand(os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak"), sample=samples.keys()),
        expand(os.path.join(result_path, 'report', '{sample}_peaks.xls'), sample=samples.keys()), # representing symlinked stats
        sample_annotation = config["annotation"],
    output:
        multiqc_report = report(os.path.join(result_path,"report","multiqc_report.html"),
                                caption="../report/multiqc.rst",
                                category="{}_{}".format(config["project_name"], module_name),
                                subcategory="QC",
                                labels={
                                    "name": "MultiQC report",
                                    "type": "HTML",
                                }),
        multiqc_stats = os.path.join(result_path, "report", "multiqc_report_data", "multiqc_general_stats.txt"),
    params:
        result_path = result_path,
        multiqc_configs = "{{'title': '{name}', 'intro_text': 'Quality Control Metrics of the ATAC-seq pipeline.', 'fastp': {{'s_name_filenames': true}}, 'annotation': '{annot}', 'genome': '{genome}', 'exploratory_columns': {exploratory_columns}, 'skip_versions_section': true}}".format(name = config["project_name"], annot = config["annotation"], genome = config["genome"], exploratory_columns = config["annot_columns"]),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/multiqc.yaml",
    log:
        "logs/rules/multiqc.log"
    shell:
        """
        multiqc {params.result_path}/report --force --verbose --outdir {params.result_path}/report --filename multiqc_report.html --cl-config "{params.multiqc_configs}"
        """

# visualize sample annotation (including QC metrics)
rule plot_sample_annotation:
    input:
        sample_annotation = config["annotation"],
        sample_annotation_w_QC = os.path.join(result_path, "counts", "sample_annotation.csv"),
    output:
        sample_annotation_plot = os.path.join(result_path,"report","sample_annotation.png"),
        sample_annotation_html = report(os.path.join(result_path,"report","sample_annotation.html"),
                       caption="../report/sample_annotation.rst",
                       category="{}_{}".format(config["project_name"], module_name),
                       subcategory="QC",
                       labels={
                           "name": "Sample annotation",
                           "type": "HTML",
                           }),
    log:
        "logs/rules/plot_sample_annotation.log",
    resources:
        mem_mb="4000",
    conda:
        "../envs/ggplot.yaml"
    script:
        "../scripts/plot_sample_annotation.R"