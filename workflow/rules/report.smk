
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
    params:
        # cluster parameters
        partition=config.get("partition"),
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

# rule ucsc_hub:
#     input:
#         bigwig_files = expand(os.path.join(result_path, "hub", "{sample}.bigWig"), sample=samples.keys()),
#     output:
#         bigwig_symlinks = expand(os.path.join(result_path, "hub", config["genome"], "{sample}.bigWig"), sample=samples.keys()),
#         genomes_file = os.path.join(result_path, "hub", "genomes.txt"),
#         hub_file = os.path.join(result_path, "hub", "hub.txt"),
#         trackdb_file = os.path.join(result_path, "hub", config["genome"], "trackDb.txt"),
#     params:
#         # cluster parameters
#         partition=config.get("partition"),
#     resources:
#         mem_mb=config.get("mem", "1000"),
#     threads: config.get("threads", 1)
#     log:
#         "logs/rules/ucsc_hub.log"
#     run:
#         # create bigwig symlinks
#         for i in range(len(input.bigwig_files)):
#             os.symlink(os.path.join('../',os.path.basename(input.bigwig_files[i])), output.bigwig_symlinks[i])

#         # create genomes.txt
#         with open(output.genomes_file, 'w') as gf:
#             genomes_text = f'genome {config["genome"]}\ntrackDb {config["genome"]}/trackDb.txt\n'
#             gf.write(genomes_text)

#         # create hub file
#         with open(output.hub_file, 'w') as hf:
#             hub_text = [f'hub {config["project_name"]}',
#                         f'shortLabel {config["project_name"]}',
#                         f'longLabel {config["project_name"]}',
#                         'genomesFile genomes.txt',
#                         f'email {config["email"]}\n',]
#             hf.write('\n'.join(hub_text))

#         # create trackdb file
#         with open(output.trackdb_file, 'w') as tf:
#             colors = ['166,206,227', '31,120,180', '51,160,44', '251,154,153', '227,26,28',
#                               '253,191,111', '255,127,0', '202,178,214', '106,61,154', '177,89,40']
            
#             track_db = ['track {}'.format(config["project_name"]),
#                         'type bigWig', 'compositeTrack on', 'autoScale on', 'maxHeightPixels 32:32:8',
#                         'shortLabel {}'.format(config["project_name"][:8]),
#                         'longLabel {}'.format(config["project_name"]),
#                         'visibility full',
#                         '', '']
#             for sample_name in samples.keys():
#                 track_color = '255,40,0'
                
#                 if config["annot_columns"][0]!="":
#                     color_hash = hash(samples[sample_name][config["annot_columns"][0]])
#                     track_color = colors[color_hash % len(colors)]
                
#                 track = ['track {}'.format(sample_name),
#                          'shortLabel {}'.format(sample_name),
#                          'longLabel {}'.format(sample_name),
#                          'bigDataUrl {}.bigWig'.format(sample_name),
#                          'parent {} on'.format(config["project_name"]),
#                          'type bigWig', 'windowingFunction mean',
#                          'color {}'.format(track_color),
#                          '', '']
                
#                 track_db += track

#             tf.write('\n'.join(track_db))
            
rule multiqc:
    input:
        expand(os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"), sample=samples.keys()),
        expand(os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak"), sample=samples.keys()),
#         expand(os.path.join(result_path, "hub","{sample}.bigWig"),sample=samples.keys()),
        expand(os.path.join(result_path, 'report', '{sample}_peaks.xls'), sample=samples.keys()), # representing symlinked stats
        trackdb_file = os.path.join(result_path, "hub", config["genome"], "trackDb.txt"), # representing UCSC hub
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
        multiqc {params.result_path}/report --force --verbose --outdir {params.result_path}/report --filename multiqc_report.html --cl-config "{params.multiqc_configs}"
        """
