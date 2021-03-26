#!/usr/bin/env python
""" MultiQC BSF plugin functions
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""


from __future__ import print_function
from pkg_resources import get_distribution
import logging
import os, csv

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.atacseq_report_version = get_distribution("atacseq_report").version


# Add default config options for the things that are used in atacseq_report
def atacseq_report_execution_start():
    """
    Code to execute after the config files and
    command line flags have been parsed self.
    this setuptools hook is the earliest that will be able
    to use custom command line flags.
    """
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_atacseq_report', True):
        return None

    log.info("Running atacseq_report MultiQC Plugin v{}, use --disable-atacseq-report to disable".format(config.atacseq_report_version))

    # Add to the search patterns used by atacseq module
    if 'atacseq' not in config.sp:
        config.update_dict(config.sp, {'atacseq': {'fn': '*.stats.tsv', 'contents': 'frip'}})
        log.info("updated config.sp for atacseq")
    if 'atacseq/tss' not in config.sp:
        config.update_dict(config.sp, {'atacseq/tss': {'fn': '*TSS.csv', 'contents': 'count'}})

    # Create symlink for the web server
    if hasattr(config, 'base_url') and hasattr(config, 'project_uuid') and hasattr(config, 'public_html_folder'):
        project_url = os.path.join(config.base_url, config.project_uuid)
        os.chdir(config.public_html_folder)
        if not os.path.exists(os.path.join(config.public_html_folder, config.project_uuid)):
            # The symlink has to be relative so that the web server can locate the project folder
            relative_path = os.path.relpath(config.project_path)
            os.symlink(relative_path, config.project_uuid)
        log.info('## You can access the project report from: ##\n{}\n'.format(os.path.join(project_url,
                                                                                          'atacseq_report',
                                                                                          'multiqc_report.html')))
    else:
        log.error('Please provide base_url, project_uuid and public_html_folder in the configuration file')
        exit(1)


    # Setup ATACseq report folder and UCSC track hub
    if hasattr(config, 'sample_annotation'):
        with open(config.sample_annotation, 'r') as sas:
            sas_reader = csv.DictReader(sas)
            samples_dict = {}
            for row in sas_reader:
                if 'sample_name' in row and row['sample_name'] not in samples_dict:
                    samples_dict[row['sample_name']] = row
            log.info('There were {} samples in the sample annotation sheet'.format(len(samples_dict)))
            report_dir = os.path.join(config.project_path, 'atacseq_report')
            if not os.path.exists(report_dir):
                os.mkdir(report_dir)
            config.output_dir = report_dir
            config.analysis_dir = [report_dir]
            os.chdir(report_dir)
            # Create symbolic links to relevant pipeline output files for use in report generation
            for sample_name in samples_dict:
                source_path = os.path.join('../',
                                           'atacseq_results',
                                           sample_name,
                                           '{}.stats.tsv'.format(sample_name))
                if not os.path.islink('{}.stats.tsv'.format(sample_name)):
                    os.symlink(source_path, '{}.stats.tsv'.format(sample_name), )
                source_path = os.path.join('../',
                                        'atacseq_results',
                                        sample_name,
                                        '{}.tss_histogram.csv'.format(sample_name))
                if not os.path.islink('{}_TSS.csv'.format(sample_name)):
                    os.symlink(source_path, '{}_TSS.csv'.format(sample_name))
                source_path = os.path.join('../',
                                           'atacseq_results',
                                           sample_name, 'mapped',
                                           '{}.txt'.format(sample_name))
                if not os.path.islink('{}.txt'.format(sample_name)):
                    os.symlink(source_path, '{}.txt'.format(sample_name))
                source_path = os.path.join('../',
                                        'atacseq_results',
                                        sample_name, 'mapped',
                                        '{}.fastp.json'.format(sample_name))
                if not os.path.islink('{}.fastp.json'.format(sample_name)):
                    os.symlink(source_path, '{}.fastp.json'.format(sample_name))
                source_path = os.path.join('../',
                                          'atacseq_results',
                                          sample_name, 'mapped',
                                          '{}.samblaster.log'.format(sample_name))
                if not os.path.islink('{}.samblaster.log'.format(sample_name)):
                    os.symlink(source_path, '{}.samblaster.log'.format(sample_name))
                source_path = os.path.join('../',
                                           'atacseq_results',
                                           sample_name, 'mapped',
                                           '{}.samtools_flagstat.log'.format(sample_name))
                if not os.path.islink('{}.samtools_flagstat.log'.format(sample_name)):
                    os.symlink(source_path, '{}.samtools_flagstat.log'.format(sample_name))
                source_path = os.path.join('../',
                                           'atacseq_results',
                                           sample_name, 'peaks',
                                           '{}.macs2.log'.format(sample_name))
                if not os.path.islink('{}.macs2.log'.format(sample_name)):
                    os.symlink(source_path, '{}.macs2.log'.format(sample_name))
                source_path = os.path.join('../',
                                           'atacseq_results',
                                           sample_name, 'peaks',
                                           '{}_peaks.xls'.format(sample_name))
                if not os.path.islink('{}_peaks.xls'.format(sample_name)):
                    os.symlink(source_path, '{}_peaks.xls'.format(sample_name))
            # Create UCSC track hub
            if hasattr(config, 'trackhub_dir'):
                hub_dir = os.path.join(config.project_path, config.trackhub_dir) # os.path.join(config.metadata['output_dir'], 'atacseq_hub')
                if not os.path.exists(hub_dir):
                    log.error('Please make sure that trackhub_dir exists')
                track_dir = os.path.join(hub_dir, config.genome)
                if not os.path.exists(track_dir):
                    os.mkdir(track_dir)
                os.chdir(track_dir)
                # Create the bigWig links for the sample coverage tracks
                for sample_name in samples_dict:
                    bigWig_path = os.path.join('../',
                                               '{}.bigWig'.format(sample_name))
                    if not os.path.islink('{}.bigWig'.format(sample_name)):
                        os.symlink(bigWig_path, '{}.bigWig'.format(sample_name))
                genomes_file_path = os.path.join(hub_dir, 'genomes.txt')
                with open(genomes_file_path, 'w') as genomes_file:
                    genomes_text = 'genome {}\ntrackDb {}/trackDb.txt\n'.format(config.genome,
                                                                              config.genome)
                    genomes_file.write(genomes_text)
                hub_file_path = os.path.join(hub_dir, 'hub.txt')
                with open(hub_file_path, 'w') as hub_file:
                    hub_text = ['hub {}'.format(config.trackhub_name),
                                'shortLabel {}'.format(config.trackhub_name),
                                'longLabel {}'.format(config.trackhub_name),
                                'genomesFile genomes.txt',
                                'email {}\n'.format(config.email)]
                    hub_file.write('\n'.join(hub_text))

                trackdb_file_path = os.path.join(hub_dir,
                                                 config.genome,
                                                 'trackDb.txt')
                with open(trackdb_file_path, 'w') as trackdb_file:
                    colors = ['166,206,227', '31,120,180', '51,160,44', '251,154,153', '227,26,28',
                              '253,191,111', '255,127,0', '202,178,214', '106,61,154', '177,89,40']
                    if hasattr(config, 'trackhub_color_by'):
                        color_groups = []
                        for sample_name in samples_dict:
                            if samples_dict[sample_name][config.trackhub_color_by] not in color_groups:
                                color_groups.append(samples_dict[sample_name][config.trackhub_color_by])

                    track_db = ['track {}'.format(config.trackhub_name),
                                'type bigWig', 'compositeTrack on', 'autoScale on', 'maxHeightPixels 32:32:8',
                                'shortLabel {}'.format(config.trackhub_name[:8]),
                                'longLabel {}'.format(config.trackhub_name),
                                'visibility {}'.format(config.trackhub_visibility),
                                '', '']
                    for sample_name in samples_dict:
                        short_label = sample_name
                        if hasattr(config,'trackhub_short_label_column'):
                            short_label = samples_dict[sample_name][config.trackhub_short_label_column]
                        track_color = '255,40,0'
                        if hasattr(config, 'trackhub_color_by'):
                            color_hash = hash(samples_dict[sample_name][config.trackhub_color_by])
                            track_color = colors[color_hash % len(colors)]
                        track = ['track {}'.format(sample_name),
                                 'shortLabel {}'.format(short_label),
                                 'longLabel {}'.format(sample_name),
                                 'bigDataUrl {}.bigWig'.format(sample_name),
                                 'parent {} on'.format(config.trackhub_name),
                                 'type bigWig', 'windowingFunction mean',
                                 'color {}'.format(track_color),
                                 '', '']
                        track_db += track
                    trackdb_file.write('\n'.join(track_db))
            else:
                log.warning('Trackhubs configuration is missing!')
        # Finally, switch back to the report directory for scanning the stats files
        os.chdir(report_dir)
    else:
        log.error('Please provide the location of the ATACseq sample annotation sheet in the configuration file')
        exit(1)
