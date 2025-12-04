#!/usr/bin/env python
"""
MultiQC command line options - we tie into the MultiQC
core here and add some new command line parameters.
See the Click documentation for more command line flag types:
http://click.pocoo.org/5/
"""

import click

# Sets config.kwargs['disable_plugin'] to True if specified (will be False otherwise)
disable_cuttag_report = click.option('--disable-cuttag-report', 'disable_cuttag_report',
    is_flag = True,
    help = "Disable the CUTTAG-seq pipeline report"
)