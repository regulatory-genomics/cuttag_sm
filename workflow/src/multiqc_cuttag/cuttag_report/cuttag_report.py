#!/usr/bin/env python
""" MultiQC CUT&Tag plugin functions
Custom Python functions for CUT&Tag pipeline MultiQC integration.
"""

from __future__ import print_function
from pkg_resources import get_distribution
from collections import OrderedDict
import logging
import os
import re
import csv

from multiqc.utils import config

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.cuttag_report_version = get_distribution("cuttag_report").version


def cuttag_report_execution_start():
    """
    Code to execute after config files and command line flags have been parsed.
    This is the earliest hook that can use custom command line flags.
    """
    if config.kwargs.get("disable_cuttag_report", False):
        return None

    log.info(
        "Running cuttag_report MultiQC Plugin v{}, use --disable-cuttag-report to disable".format(
            config.cuttag_report_version
        )
    )

    # Load annotation path from config
    _load_annotation_config()
    
    # Set up search patterns
    _setup_search_patterns()
    
    # Set up custom data sections
    _setup_custom_data()
    
    # Set up sample name renaming
    _setup_sample_name_renaming()


def _load_annotation_config():
    """Load annotation path from config files."""
    if hasattr(config, "annotation") and config.annotation:
        return
    
    log.debug("Attempting to load annotation from config...")
    
    # Check nondefault_config dict
    if hasattr(config, "nondefault_config") and isinstance(config.nondefault_config, dict):
        if 'annotation' in config.nondefault_config:
            config.annotation = config.nondefault_config['annotation']
            log.info(f"Loaded annotation from config: {config.annotation}")
            return
    
    # Check explicit config files
    config_files = getattr(config, "explicit_user_config_files", [])
    for config_file in config_files:
        try:
            import yaml
            with open(config_file, 'r') as f:
                cfg_data = yaml.safe_load(f)
                if cfg_data and 'annotation' in cfg_data:
                    config.annotation = cfg_data['annotation']
                    log.info(f"Loaded annotation from {config_file}: {config.annotation}")
                    return
        except Exception as e:
            log.debug(f"Could not load annotation from {config_file}: {e}")


def _setup_search_patterns():
    """Set up file search patterns for MultiQC modules."""
    if "cuttag" not in config.sp:
        config.update_dict(
            config.sp,
            {"cuttag": {"fn": "*.stats.tsv", "contents": "frip"}},
        )
        log.info("Updated config.sp for 'cuttag' stats files")

    config.update_dict(
        config.sp,
        {
            "peaks": {"fn": "*.counts.txt"},
            "fraglen": {"fn": "*.fraglen.tsv"},
            "diffpeaks": {"fn": "*diffexp-summary.txt"},
        },
    )


def _setup_custom_data():
    """Set up custom data sections for MultiQC report."""
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


def _setup_sample_name_renaming():
    """Set up sample name renaming using MultiQC's general_stats_rename config."""
    sample_mapping = _load_sample_mapping()
    if not sample_mapping:
        return
    
    if not hasattr(config, "general_stats_rename"):
        config.general_stats_rename = {}
    
    # Add all mappings to general_stats_rename
    config.general_stats_rename.update(sample_mapping)
    log.info(f"Set up {len(sample_mapping)} sample name mappings")


def _load_sample_mapping():
    """
    Load sample sheet and create a mapping from log-derived names to canonical sample names.
    Returns a dict mapping log names -> canonical sample names.
    """
    annotation_path = getattr(config, "annotation", None)
    
    if not annotation_path:
        log.debug("No annotation file configured for sample name normalization")
        return {}
    
    if not os.path.exists(annotation_path):
        log.warning(f"Annotation file not found: {annotation_path}")
        return {}
    
    sample_mapping = {}
    
    try:
        with open(annotation_path, "r") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sample_name = row.get("sample_name", "").strip() or row.get("sample", "").strip()
                if not sample_name:
                    continue
                
                run = row.get("run", "").strip() or "1"
                
                # Map various log file name patterns to canonical sample name
                patterns = [
                    f"bowtie2_{sample_name}.{run}",
                    f"bowtie2_{sample_name}.{run}.err",
                    f"bowtie2_{sample_name}",
                    f"bwa_mem2_{sample_name}.{run}",
                    f"bwa_mem2_{sample_name}.{run}.err",
                    f"bwa_mem2_{sample_name}",
                    f"sambamba_markdup_{sample_name}",
                    f"sambamba_markdup_{sample_name}.log",
                    f"sambamba_sort_{sample_name}.{run}",
                    f"sambamba_sort_{sample_name}.{run}.log",
                    f"sambamba_sort_{sample_name}",
                    f"samtools_merge_{sample_name}",
                    f"samtools_merge_{sample_name}.log",
                    f"samtools_index_{sample_name}",
                    f"samtools_index_{sample_name}.log",
                    f"{sample_name}.{run}_fastp",
                    f"{sample_name}.{run}.fastp",
                    f"{sample_name}.{run}",
                ]
                
                for pattern in patterns:
                    sample_mapping[pattern] = sample_name
        
        log.info(f"Loaded sample name mapping for {len(set(sample_mapping.values()))} samples")
    except Exception as exc:
        log.warning(f"Could not load annotation file: {exc}")
    
    return sample_mapping


def cuttag_report_after_modules():
    """
    Hook that runs after all modules are initialized.
    Injects custom metrics and normalizes sample names.
    """
    if config.kwargs.get("disable_cuttag_report", False):
        return None
    
    try:
        from multiqc.utils import report
        
        # Inject custom metrics
        _inject_deeptools_enrichment(report)
        _add_sambamba_calculated_metrics(report)
        
        # Normalize sample names
        _normalize_sample_names(report)
    
    except Exception as exc:
        import traceback
        log.error(f"Error in after_modules hook: {exc}")
        log.debug(traceback.format_exc())


def _inject_deeptools_enrichment(report):
    """
    Parse raw DeepTools plotEnrichment output (frip_*.tsv) and add to general stats.
    """
    # Search for frip_*.tsv files
    frip_files = _find_frip_files()
    if not frip_files:
        log.debug("No frip_*.tsv files found")
        return
    
    log.info(f"Found {len(frip_files)} FRiP files")
    
    # Parse FRiP files
    data = {}
    for frip_file in frip_files:
        try:
            with open(frip_file, "r") as fh:
                lines = [l.strip() for l in fh.readlines() if l.strip()]
            
            if len(lines) < 2:
                continue
            
            header = lines[0].split("\t")
            try:
                percent_idx = header.index("percent")
            except ValueError:
                log.warning(f"No 'percent' column in {frip_file}")
                continue
            
            # Extract sample name from filename (frip_SAMPLE.tsv)
            sample = os.path.basename(frip_file).replace("frip_", "").replace(".tsv", "")
            
            # Parse data line
            for line in lines[1:]:
                parts = line.split("\t")
                if len(parts) > percent_idx:
                    try:
                        val = float(parts[percent_idx])
                        data[sample] = {"deeptools_enrichment_pct": val}
                        log.info(f"Added FRiP for sample '{sample}': {val}%")
                    except ValueError:
                        continue
        except Exception as exc:
            log.warning(f"Could not parse {frip_file}: {exc}")
    
    if not data:
        return
    
    # Add to general stats
    _ensure_report_structures(report)
    
    module_name = "deeptools_enrichment"
    report.general_stats_data[module_name] = data
    report.general_stats_headers[module_name] = OrderedDict()
    report.general_stats_headers[module_name]["deeptools_enrichment_pct"] = {
        "title": "FRiP %",
        "description": "DeepTools percent in features",
        "format": "{:.2f}",
        "scale": "Blues",
    }


def _find_frip_files():
    """Search for frip_*.tsv files in analysis directories."""
    search_dirs = []
    
    # Add analysis directories
    if hasattr(config, "analysis_dir") and config.analysis_dir:
        analysis_dirs = config.analysis_dir if isinstance(config.analysis_dir, list) else [config.analysis_dir]
        for adir in analysis_dirs:
            if adir:
                search_dirs.extend([
                    adir,
                    os.path.join(adir, "plotEnrichment"),
                    os.path.join(adir, "Report", "plotEnrichment"),
                ])
    
    # Add output directory
    if hasattr(config, "output_dir") and config.output_dir:
        search_dirs.extend([
            config.output_dir,
            os.path.join(config.output_dir, "plotEnrichment"),
        ])
    
    # Add current working directory
    search_dirs.extend([".", "plotEnrichment", "Report/plotEnrichment"])
    
    # Search for files
    frip_files = []
    for search_dir in search_dirs:
        if not os.path.isdir(search_dir):
            continue
        try:
            for fname in os.listdir(search_dir):
                if fname.startswith("frip_") and fname.endswith(".tsv"):
                    frip_files.append(os.path.join(search_dir, fname))
        except OSError:
            continue
    
    return frip_files


def _add_sambamba_calculated_metrics(report):
    """
    Add calculated metrics from sambamba markdup data:
    - nrf_sm: Non-Redundant Fraction = 100 - duplicate_rate
    - duplicate_read: Number of duplicate reads
    - unique_read: (sorted_end_pairs * 2 - duplicate_reads) / 2
    """
    # Find Sambamba module
    sambamba_data = _find_module_data(report, ["sambamba", "markdup"])
    if not sambamba_data:
        log.debug("Sambamba module not found")
        return
    
    module_name, module_data = sambamba_data
    log.info(f"Found Sambamba module: {module_name}")
    
    # Process each sample
    for sample_name, sample_data in module_data.items():
        sample_dict = _extract_sample_dict(sample_data)
        if not sample_dict:
            continue
        
        # Calculate nrf_sm = 100 - duplicate_rate
        if "duplicate_rate" in sample_dict:
            try:
                duplicate_rate = float(sample_dict["duplicate_rate"])
                nrf_sm = 100.0 - duplicate_rate
                sample_dict["nrf_sm"] = nrf_sm
                log.info(f"Added nrf_sm for {sample_name}: {nrf_sm:.2f}%")
            except (ValueError, TypeError) as e:
                log.debug(f"Could not calculate nrf_sm for {sample_name}: {e}")
        
        # Add duplicate_reads to general stats
        if "duplicate_reads" in sample_dict and "duplicate_read" not in sample_dict:
            sample_dict["duplicate_read"] = sample_dict["duplicate_reads"]
            log.info(f"Added duplicate_read for {sample_name}: {sample_dict['duplicate_reads']}")
        
        # Calculate unique_read
        if "sorted_end_pairs" in sample_dict and "duplicate_reads" in sample_dict:
            try:
                sorted_end_pairs = int(sample_dict["sorted_end_pairs"])
                duplicate_reads = int(sample_dict["duplicate_reads"])
                unique_read = (sorted_end_pairs * 2 - duplicate_reads) / 2
                sample_dict["unique_read"] = unique_read
                log.info(f"Added unique_read for {sample_name}: {unique_read:.0f}")
            except (ValueError, TypeError) as e:
                log.debug(f"Could not calculate unique_read for {sample_name}: {e}")
    
    # Add headers for new metrics
    _add_sambamba_headers(report, module_name)
    log.info("Successfully added calculated sambamba metrics")


def _find_module_data(report, search_terms):
    """Find module data in general_stats_data by searching for terms in module name."""
    if not hasattr(report, "general_stats_data") or not report.general_stats_data:
        return None
    
    # Handle both list and dict formats
    if isinstance(report.general_stats_data, list):
        for module_dict in report.general_stats_data:
            for module_name, module_data in module_dict.items():
                if any(term in module_name.lower() for term in search_terms):
                    return (module_name, module_data)
    else:
        for module_name, module_data in report.general_stats_data.items():
            if any(term in module_name.lower() for term in search_terms):
                return (module_name, module_data)
    
    return None


def _extract_sample_dict(sample_data):
    """Extract dict from sample data (handles InputRow objects and lists)."""
    if isinstance(sample_data, dict):
        return sample_data
    elif isinstance(sample_data, list) and len(sample_data) > 0:
        # Check for InputRow object with data attribute
        if hasattr(sample_data[0], 'data'):
            return sample_data[0].data
        elif isinstance(sample_data[0], dict):
            return sample_data[0]
    return None


def _add_sambamba_headers(report, module_name):
    """Add headers for calculated sambamba metrics."""
    _ensure_report_structures(report)
    
    if module_name not in report.general_stats_headers:
        report.general_stats_headers[module_name] = OrderedDict()
    
    headers = {
        "nrf_sm": {
            "title": "NRF",
            "description": "Non-Redundant Fraction (100 - duplicate rate)",
            "format": "{:.2f}",
            "suffix": "%",
            "scale": "RdYlGn",
            "min": 0,
            "max": 100,
        },
        "duplicate_read": {
            "title": "Dup Reads",
            "description": "Number of duplicate reads",
            "format": "{:,.0f}",
            "scale": "Reds",
        },
        "unique_read": {
            "title": "Unique Reads",
            "description": "Number of unique reads: (sorted_end_pairs * 2 - duplicate_reads) / 2",
            "format": "{:,.0f}",
            "scale": "Greens",
        },
    }
    
    for key, header in headers.items():
        if key not in report.general_stats_headers[module_name]:
            report.general_stats_headers[module_name][key] = header


def _normalize_sample_names(report):
    """Normalize sample names in general stats data."""
    sample_mapping = _load_sample_mapping()
    if not sample_mapping:
        log.debug("No sample mapping for normalization")
        return
    
    if not hasattr(report, "general_stats_data") or not report.general_stats_data:
        return
    
    canonical_samples = set(sample_mapping.values())
    log.info(f"Normalizing sample names (canonical samples: {sorted(canonical_samples)})")
    
    # Apply normalization
    renamed_count = 0
    
    if isinstance(report.general_stats_data, list):
        for module_dict in report.general_stats_data:
            for module_name, module_data in list(module_dict.items()):
                renamed_count += _rename_samples_in_module(module_data, sample_mapping, canonical_samples, module_name)
    else:
        for module_name, module_data in list(report.general_stats_data.items()):
            renamed_count += _rename_samples_in_module(module_data, sample_mapping, canonical_samples, module_name)
    
    if renamed_count > 0:
        log.info(f"Normalized {renamed_count} sample names")
    else:
        log.debug("No sample names were normalized")


def _rename_samples_in_module(module_data, sample_mapping, canonical_samples, module_name):
    """Rename samples in a single module's data."""
    renamed_count = 0
    
    for old_name in list(module_data.keys()):
        new_name = _find_canonical_name(old_name, sample_mapping, canonical_samples)
        
        if new_name and new_name != old_name:
            if new_name in module_data:
                # Merge if target exists
                _merge_sample_data(module_data, old_name, new_name)
                log.info(f"Merged '{old_name}' into '{new_name}' in {module_name}")
            else:
                # Simple rename
                module_data[new_name] = module_data.pop(old_name)
                renamed_count += 1
                log.debug(f"Renamed '{old_name}' -> '{new_name}' in {module_name}")
    
    return renamed_count


def _find_canonical_name(sample_name, sample_mapping, canonical_samples):
    """Find canonical name for a sample."""
    # Exact match
    if sample_name in sample_mapping:
        return sample_mapping[sample_name]
    
    # Try pattern extraction
    extracted = _extract_sample_from_pattern(sample_name)
    if extracted and extracted in canonical_samples:
        return extracted
    
    # Try prefix matching
    for canon in canonical_samples:
        if canon.startswith(sample_name) or sample_name.startswith(canon):
            return canon
    
    return None


def _extract_sample_from_pattern(name):
    """Extract sample name from log file patterns."""
    patterns = [
        r"^bowtie2_(.+?)\.\d+\.err?$",
        r"^bowtie2_(.+?)\.err?$",
        r"^bwa_mem2_(.+?)\.\d+\.err?$",
        r"^bwa_mem2_(.+?)\.err?$",
        r"^sambamba_markdup_(.+?)(?:\.log)?$",
        r"^sambamba_sort_(.+?)\.\d+(?:\.log)?$",
        r"^samtools_merge_(.+?)(?:\.log)?$",
        r"^samtools_index_(.+?)(?:\.log)?$",
    ]
    
    for pattern in patterns:
        match = re.match(pattern, name)
        if match:
            return match.group(1)
    
    # Handle fastp concatenated names
    if "_" in name and not any(name.startswith(p) for p in ["bowtie2_", "bwa_mem2_", "sambamba_", "samtools_"]):
        parts = name.split("_")
        if len(parts) >= 2:
            # Try to find common base by removing trailing digits
            base = re.sub(r'\d+$', '', "_".join(parts[:len(parts)//2]))
            if base:
                return base
    
    return None


def _merge_sample_data(module_data, old_name, new_name):
    """Merge data from old_name into new_name."""
    old_data = module_data[old_name]
    new_data = module_data[new_name]
    
    if isinstance(new_data, dict) and isinstance(old_data, dict):
        new_data.update(old_data)
    elif isinstance(new_data, (int, float)) and isinstance(old_data, (int, float)):
        module_data[new_name] = (new_data + old_data) / 2
    
    del module_data[old_name]


def _ensure_report_structures(report):
    """Ensure report has required data structures."""
    if not hasattr(report, "general_stats_data") or report.general_stats_data is None:
        report.general_stats_data = {}
    if not hasattr(report, "general_stats_headers") or report.general_stats_headers is None:
        report.general_stats_headers = {}


def cuttag_report_before_report_generation():
    """
    Hook that runs just before report generation.
    Final chance to normalize sample names before output.
    """
    if config.kwargs.get("disable_cuttag_report", False):
        return None
    
    try:
        from multiqc.utils import report
        
        sample_mapping = _load_sample_mapping()
        if not sample_mapping:
            return
        
        # Apply final normalization
        _normalize_sample_names(report)
        
        log.info("Completed final sample name normalization")
    
    except Exception as exc:
        import traceback
        log.error(f"Error in before_report_generation hook: {exc}")
        log.debug(traceback.format_exc())
