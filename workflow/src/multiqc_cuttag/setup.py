#!/usr/bin/env python
"""
MultiQC plugin for reporting results of the CUTTAG-seq pipeline.
"""

from setuptools import setup, find_packages

version = '0.3'

setup(
    name = 'cuttag_report',
    version = 0.1,
    author = 'Gilbert Han',
    author_email = 'hanlitian@westlake.edu.cn',
    description = "MultiQC plugin for reporting results of the CUTTAG-seq pipeline.",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/regulatory-genomics/cuttag_sm/tree/master/workflow/src/multiqc_cuttag',
    download_url = 'https://github.com/regulatory-genomics/cuttag_sm/tree/master/workflow/src/multiqc_cuttag',
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
            'cuttag = cuttag_report.modules.cuttag:MultiqcModule',
        ],
        'multiqc.cli_options.v1': [
            'disable_cuttag_report = cuttag_report.cli:disable_cuttag_report'
        ],
        'multiqc.hooks.v1': [
            'execution_start = cuttag_report.cuttag_report:cuttag_report_execution_start'
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
