#!/usr/bin/env python
"""
MultiQC plugin for reporting results of pipelines used at BSF
"""

from setuptools import setup, find_packages

version = '0.2'

setup(
    name = 'atacseq_report',
    version = 0.2,
    author = 'Bekir Erguener',
    author_email = 'berguener@cemm.at',
    description = "MultiQC plugin for reporting results of BSF's ATAC-seq pipeline",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/berguner/atacseq_pipeline',
    download_url = 'https://github.com/berguner/atacseq_pipeline',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    install_requires =[
        'multiqc',
        'click',
        'tables'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'atacseq = atacseq_report.modules.atacseq:MultiqcModule',
        ],
        'multiqc.cli_options.v1': [
            'disable_atacseq_report = atacseq_report.cli:disable_atacseq_report'
        ],
        'multiqc.hooks.v1': [
            'execution_start = atacseq_report.atacseq_report:atacseq_report_execution_start'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ]
)
