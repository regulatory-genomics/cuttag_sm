#!/usr/bin/env python
""" MultiQC BSF plugin functions
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""


from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import config

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.cuttag_report_version = get_distribution("cuttag_report").version


# Add default config options for the things that are used in cuttag_report
def cuttag_report_execution_start():
    """
    Code to execute after the config files and
    command line flags have been parsed self.
    this setuptools hook is the earliest that will be able
    to use custom command line flags.
    """
    # Halt execution if we've disabled the plugin
    if config.kwargs.get("disable_cuttag_report", True):
        return None

    log.info(
        "Running cuttag_report MultiQC Plugin v{}, use --disable-cuttag-report to disable".format(
            config.cuttag_report_version
        )
    )

    # ------------------------------------------------------------------
    # Search patterns (equivalent to 'sp:' section in multiqc_conf.yml)
    # ------------------------------------------------------------------
    if "cuttag" not in config.sp:
        config.update_dict(
            config.sp,
            {"cuttag": {"fn": "*.stats.tsv", "contents": "frip"}},
        )
        log.info("Updated config.sp for 'cuttag' stats files")

    # Peaks counts, fragment length, differential peaks
    config.update_dict(
        config.sp,
        {
            "peaks": {"fn": "*.counts.txt"},
            "fraglen": {"fn": "*.fraglen.tsv"},
            "diffpeaks": {"fn": "*diffexp-summary.txt"},
        },
    )

    # ------------------------------------------------------------------
    # Custom data sections (equivalent to 'custom_data:' in multiqc_conf.yml)
    # ------------------------------------------------------------------
    if not hasattr(config, "custom_data") or config.custom_data is None:
        config.custom_data = {}

    config.update_dict(
        config.custom_data,
        {
            "peaks": {
                "file_format": "txt",
                "section_name": "Peak Counts",
                "description": "Callpeaks peak counts",
                "plot_type": "bargraph",
                "pconfig": {
                    "id": "peakcounts bargraph",
                    "title": "Sample Peak Counts",
                    "ylab": "Number of Peaks",
                    "xlab": "Sample",
                },
            },
            "fraglen": {
                "file_format": "txt",
                "section_name": "Fragment length distribution",
                "description": "Sample Bam fragment length distribution",
                "plot_type": "linegraph",
                "pconfig": {
                    "id": "fraglen dist",
                    "title": "sample fragment length distribution",
                    "xlab": "fragment length",
                    "ylab": "counts",
                },
            },
            "diffpeaks": {
                "file_format": "txt",
                "section_name": "Differential peaks",
                "description": "Differential peaks for each condition",
                "plot_type": "bargraph",
                "pconfig": {
                    "id": "diff peaks bargraph",
                    "title": "Condition Differential Peaks",
                    "ylab": "Number of differential peaks",
                    "xlab": "condition",
                },
            },
        },
    )

