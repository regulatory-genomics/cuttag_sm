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
        if config.kwargs.get("disable_cuttag_report", False):
            return None

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="The CUT&Tag Pipeline",
            anchor="cuttag",
            href="https://github.com/regulatory-genomics/cuttag_sm",
            info="Processes, quantifies and annotates CUT&Tag data.",
        )

        log.info("Initialized cuttag module")

        # Try to load optional sample annotation sheet FIRST (to enable sample name normalization)
        self.sample_sas_dict = {}
        self.sample_name_map = {}  # Maps log file names to canonical sample names
        annotation_path = getattr(config, "annotation", None)
        if annotation_path is not None:
            try:
                with open(annotation_path, "r") as fh:
                    sample_sas = csv.DictReader(fh)
                    for row in sample_sas:
                        sample_name = row.get("sample_name", "").strip() or row.get("sample", "").strip()
                        if sample_name:
                            self.sample_sas_dict[sample_name] = row
                            # Build mapping for common log file name patterns
                            run = row.get("run", "").strip() or "1"
                            # Map various patterns to canonical sample name
                            self.sample_name_map[f"{sample_name}_{run}"] = sample_name
                            self.sample_name_map[f"{sample_name}{run}_{sample_name}{run}"] = sample_name
                            self.sample_name_map[sample_name] = sample_name
                log.info("Loaded sample annotation for {} samples".format(len(self.sample_sas_dict)))
            except Exception as exc:
                log.warning(
                    "Could not load annotation file '{}': {}".format(
                        annotation_path, str(exc)
                    )
                )

        # Parse CUT&Tag stats for each sample from '*.stats.tsv'
        self.cuttag_data = {}
        for f in self.find_log_files(sp_key="cuttag"):
            # Normalize the sample name using the annotation mapping
            original_name = f["s_name"]
            normalized_name = self._normalize_sample_name(original_name)
            self.cuttag_data[normalized_name] = self.parse_cuttag_stats(f["f"])
            if original_name != normalized_name:
                log.debug(f"Normalized sample name: '{original_name}' -> '{normalized_name}'")
        log.info(
            "Found stats file for {} CUT&Tag samples".format(len(self.cuttag_data))
        )

        # Halt module if no stats were found
        if len(self.cuttag_data) == 0:
            raise UserWarning

        # Remove ignored samples if there is any
        self.cuttag_data = self.ignore_samples(self.cuttag_data)

        # Add stats to general table
        self.add_cuttag_to_general_stats()

    def _normalize_sample_name(self, log_name):
        """
        Normalize sample names from log files to canonical sample names.
        Handles patterns like:
        - sample1_sample2 (fastp concatenated) -> sample
        - bowtie2_sample.1.err -> sample
        - sambamba_markdup_sample -> sample
        """
        import re
        
        # Try direct mapping first
        if log_name in self.sample_name_map:
            return self.sample_name_map[log_name]
        
        # Extract sample name from common log file patterns
        patterns = [
            (r"^bowtie2_(.+?)\.\d+\.err?$", 1),
            (r"^bowtie2_(.+?)\.err?$", 1),
            (r"^bwa_mem2_(.+?)\.\d+\.err?$", 1),
            (r"^bwa_mem2_(.+?)\.err?$", 1),
            (r"^sambamba_markdup_(.+?)(?:\.log)?$", 1),
            (r"^sambamba_sort_(.+?)\.\d+(?:\.log)?$", 1),
            (r"^samtools_merge_(.+?)(?:\.log)?$", 1),
            (r"^samtools_index_(.+?)(?:\.log)?$", 1),
        ]
        
        for pattern, group_idx in patterns:
            match = re.match(pattern, log_name)
            if match:
                extracted = match.group(group_idx)
                # Try to find canonical name for extracted sample
                if extracted in self.sample_name_map:
                    return self.sample_name_map[extracted]
                # Check if it's already a canonical sample name
                if extracted in self.sample_sas_dict:
                    return extracted
                return extracted
        
        # Handle fastp concatenated names like "sample1_sample2"
        if "_" in log_name and not any(log_name.startswith(prefix) for prefix in ["bowtie2_", "bwa_mem2_", "sambamba_", "samtools_"]):
            parts = log_name.split("_")
            if len(parts) >= 2:
                # Try to find common prefix
                first_part = "_".join(parts[:len(parts)//2])
                second_part = "_".join(parts[len(parts)//2:])
                first_clean = re.sub(r'\d+$', '', first_part)
                second_clean = re.sub(r'\d+$', '', second_part)
                if first_clean == second_clean and first_clean in self.sample_sas_dict:
                    return first_clean
                # Try checking if any part matches a known sample
                for i in range(1, len(parts) + 1):
                    candidate = "_".join(parts[:i])
                    candidate_clean = re.sub(r'\d+$', '', candidate)
                    if candidate_clean in self.sample_sas_dict:
                        return candidate_clean
        
        # Return original if no normalization found
        return log_name

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

        self.general_stats_addcols(data, headers)


