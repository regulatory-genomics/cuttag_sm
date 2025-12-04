#!/usr/bin/env python

"""
MultiQC module to parse CUT&Tag pipeline stats
"""

from __future__ import print_function

from collections import OrderedDict
import logging
import csv

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger("multiqc")


class MultiqcModule(BaseMultiqcModule):
    """
    CUT&Tag module class
    """

    def __init__(self):
        # Halt execution if we've disabled the plugin
        if config.kwargs.get("disable_cuttag_report", True):
            return None

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="The CUT&Tag Pipeline",
            anchor="cuttag",
            href="https://github.com/regulatory-genomics/cuttag_sm",
            info="Processes, quantifies and annotates CUT&Tag data.",
        )

        log.info("Initialized cuttag module")

        # Parse CUT&Tag stats for each sample from '*.stats.tsv'
        self.cuttag_data = {}
        for f in self.find_log_files(sp_key="cuttag"):
            self.cuttag_data[f["s_name"]] = self.parse_cuttag_stats(f["f"])
        log.info(
            "Found stats file for {} CUT&Tag samples".format(len(self.cuttag_data))
        )

        # Halt module if no stats were found
        if len(self.cuttag_data) == 0:
            raise UserWarning

        # Remove ignored samples if there is any
        self.cuttag_data = self.ignore_samples(self.cuttag_data)

        # Try to load optional sample annotation sheet (safe to skip if missing)
        self.sample_sas_dict = {}
        annotation_path = getattr(config, "annotation", None)
        if annotation_path is not None:
            try:
                with open(annotation_path, "r") as fh:
                    sample_sas = csv.DictReader(fh)
                    for row in sample_sas:
                        if "sample_name" in row:
                            self.sample_sas_dict[row["sample_name"]] = row
            except Exception as exc:
                log.warning(
                    "Could not load annotation file '{}': {}".format(
                        annotation_path, str(exc)
                    )
                )

        # Add stats to general table
        self.add_cuttag_to_general_stats()

    def parse_cuttag_stats(self, fhandle):
        data = {}
        for line in fhandle.splitlines():
            if not line:
                continue
            s = line.split("\t")
            if len(s) >= 2:
                data[s[0]] = s[1]
        return data

    def add_cuttag_to_general_stats(self):
        data = {}
        for sample_name in self.cuttag_data:
            data[sample_name] = {}

            if hasattr(config, "exploratory_columns"):
                for column in config.exploratory_columns:
                    if column in self.sample_sas_dict.get(sample_name, {}):
                        data[sample_name][column] = self.sample_sas_dict[sample_name][
                            column
                        ]

            # Copy and type-cast metrics
            if (
                "NSC" in self.cuttag_data[sample_name]
                and self.cuttag_data[sample_name]["NSC"] != "nan"
            ):
                try:
                    value = float(self.cuttag_data[sample_name]["NSC"])
                except ValueError as err:  # pragma: no cover - defensive
                    log.debug(err)
                    value = None
                data[sample_name]["NSC"] = value

            if (
                "RSC" in self.cuttag_data[sample_name]
                and self.cuttag_data[sample_name]["RSC"] != "nan"
            ):
                try:
                    value = float(self.cuttag_data[sample_name]["RSC"])
                except Exception:  # pragma: no cover - defensive
                    value = None
                data[sample_name]["RSC"] = value

            if "peaks" in self.cuttag_data[sample_name]:
                try:
                    value = int(self.cuttag_data[sample_name]["peaks"])
                except Exception:
                    value = None
                data[sample_name]["peaks"] = value

            if "filtered_peaks" in self.cuttag_data[sample_name]:
                try:
                    value = int(self.cuttag_data[sample_name]["filtered_peaks"])
                except Exception:
                    value = None
                data[sample_name]["filtered_peaks"] = value

            if "frip" in self.cuttag_data[sample_name]:
                try:
                    value = float(self.cuttag_data[sample_name]["frip"])
                except Exception:
                    value = None
                data[sample_name]["frip"] = value

            if "regulatory_fraction" in self.cuttag_data[sample_name]:
                try:
                    value = float(self.cuttag_data[sample_name]["regulatory_fraction"])
                except Exception:
                    value = None
                data[sample_name]["regulatory_fraction"] = value

            if "tss_max" in self.cuttag_data[sample_name]:
                data[sample_name]["tss_max"] = self.cuttag_data[sample_name]["tss_max"]

            # The following fields are kept minimal and only use values
            # present in the main '*.stats.tsv' files.

        headers = OrderedDict()
        if hasattr(config, "exploratory_columns"):
            for column in config.exploratory_columns:
                log.info("Adding exploratory column {}".format(column))
                headers[column] = {
                    "description": column,
                    "title": column,
                    "scale": False,
                }
        else:
            log.warning("No exploratory columns were specified in the config")

        # Get list of columns to show/hide from config (if provided)
        visible_columns = getattr(config, "cuttag_general_stats_columns", None)
        if visible_columns is not None:
            log.info(
                "Using config-specified columns for CUT&Tag General Statistics: {}".format(
                    visible_columns
                )
            )

        # Helper function to conditionally add columns
        def add_header_if_visible(key, header_config):
            """Add header only if visible_columns is None (all visible) or key is in visible_columns"""
            if visible_columns is None or key in visible_columns:
                headers[key] = header_config
            else:
                # Still add but mark as hidden
                header_config["hidden"] = True
                headers[key] = header_config

        add_header_if_visible(
            "peaks",
            {
                "description": "Number of detected peaks",
                "title": "Peaks",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
        )

        add_header_if_visible(
            "filtered_peaks",
            {
                "description": "Number of peaks remaining after filtering",
                "title": "Filtered\nPeaks",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
        )

        add_header_if_visible(
            "NSC",
            {
                "description": "Normalized Strand Cross-correlation Coefficient",
                "title": "NSC",
                "scale": "Reds",
                "min": 0.0,
                "max": 2.0,
                "format": "{:.2f}",
            },
        )

        add_header_if_visible(
            "NSC_PCT",
            {
                "description": "NSC Percentile Among All CUT&Tag samples",
                "title": "NSC_PCT",
                "scale": "Reds",
                "suffix": "%",
                "max": 100,
                "format": "{:,.0f}",
            },
        )

        add_header_if_visible(
            "RSC",
            {
                "description": "Relative Strand Cross-correlation Coefficient",
                "title": "RSC",
                "scale": "Reds",
                "min": 0.0,
                "max": 2.0,
                "format": "{:.2f}",
            },
        )

        add_header_if_visible(
            "RSC_PCT",
            {
                "description": "RSC Percentile Among All CUT&Tag samples",
                "title": "RSC_PCT",
                "scale": "Reds",
                "suffix": "%",
                "max": 100,
                "format": "{:,.0f}",
            },
        )

        add_header_if_visible(
            "frip",
            {
                "description": "Fraction of Reads in Peaks",
                "title": "FRiP",
                "scale": "Reds-rev",
                "min": 0.0,
                "max": 1.0,
                "format": "{:.2f}",
            },
        )

        add_header_if_visible(
            "regulatory_fraction",
            {
                "description": "Fraction of Reads in Regulatory Regions",
                "title": "Regulatory",
                "scale": "Reds-rev",
                "min": 0.0,
                "max": 1.0,
                "format": "{:.2f}",
            },
        )

        add_header_if_visible(
            "tss_max",
            {
                "description": "The peak value of TSS enrichment",
                "title": "TSS",
                "scale": "Reds-rev",
                "format": "{:.1f}",
            },
        )

        self.general_stats_addcols(data, headers)


