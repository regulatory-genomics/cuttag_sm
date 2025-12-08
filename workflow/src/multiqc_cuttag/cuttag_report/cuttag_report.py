#!/usr/bin/env python
""" MultiQC BSF plugin functions
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""


from __future__ import print_function
from pkg_resources import get_distribution
import logging
import os
import re
import csv

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
    if config.kwargs.get("disable_cuttag_report", False):
        return None

    log.info(
        "Running cuttag_report MultiQC Plugin v{}, use --disable-cuttag-report to disable".format(
            config.cuttag_report_version
        )
    )
    
    # Load annotation path from config files if not already set
    # MultiQC loads YAML config but doesn't always set custom keys as attributes
    if not hasattr(config, "annotation") or config.annotation is None:
        log.debug("Attempting to load annotation from nondefault_config...")
        
        # Check nondefault_config dict which contains loaded config values
        if hasattr(config, "nondefault_config") and isinstance(config.nondefault_config, dict):
            if 'annotation' in config.nondefault_config:
                config.annotation = config.nondefault_config['annotation']
                log.info(f"Loaded annotation from nondefault_config: {config.annotation}")
            else:
                log.debug(f"nondefault_config keys: {list(config.nondefault_config.keys())}")
        
        # Also check explicit_user_config_files and try to load from them
        if not hasattr(config, "annotation") or config.annotation is None:
            config_files = getattr(config, "explicit_user_config_files", [])
            log.debug(f"Found {len(config_files)} explicit config files: {config_files}")
            
            for config_file in config_files:
                try:
                    import yaml
                    log.debug(f"Checking {config_file} for annotation...")
                    with open(config_file, 'r') as f:
                        cfg_data = yaml.safe_load(f)
                        if cfg_data and 'annotation' in cfg_data:
                            config.annotation = cfg_data['annotation']
                            log.info(f"Loaded annotation path from {config_file}: {config.annotation}")
                            break
                except Exception as e:
                    log.debug(f"Could not load annotation from {config_file}: {e}")

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
    
    # Set up sample name normalization via general_stats_rename
    # This must be done before modules run
    _setup_sample_name_renaming()


def _setup_sample_name_renaming():
    """
    Set up sample name renaming using MultiQC's general_stats_rename config.
    This must be called before modules run.
    """
    sample_mapping = _load_sample_mapping()
    if not sample_mapping:
        return
    
    if not hasattr(config, "general_stats_rename"):
        config.general_stats_rename = {}
    
    canonical_samples = set(sample_mapping.values())
    
    def extract_sample_from_log_name(log_name):
        """Extract sample name from various log file name patterns"""
        patterns = [
            (r"^bowtie2_(.+?)\.\d+\.err?$", 1),
            (r"^bowtie2_(.+?)\.err?$", 1),  # Without run number
            (r"^bwa_mem2_(.+?)\.\d+\.err?$", 1),
            (r"^bwa_mem2_(.+?)\.err?$", 1),  # Without run number
            (r"^sambamba_markdup_(.+?)(?:\.log)?$", 1),
            (r"^sambamba_sort_(.+?)\.\d+(?:\.log)?$", 1),
            (r"^samtools_merge_(.+?)(?:\.log)?$", 1),
            (r"^samtools_index_(.+?)(?:\.log)?$", 1),
        ]
        
        for pattern, group_idx in patterns:
            match = re.match(pattern, log_name)
            if match:
                return match.group(group_idx)
        
        # Handle fastp concatenated sample names like "P6-1_test1_P6-1_test2"
        # This happens when fastp merges multiple runs
        # Pattern: {sample}{run}_{sample}{run} -> extract {sample}
        if "_" in log_name and not any(log_name.startswith(prefix) for prefix in ["bowtie2_", "bwa_mem2_", "sambamba_", "samtools_"]):
            # Try to extract common prefix from concatenated pattern
            parts = log_name.split("_")
            if len(parts) >= 2:
                # Check if it looks like concatenated samples: sample1_sample2
                # Try to find the longest common prefix
                first_part = "_".join(parts[:len(parts)//2])
                second_part = "_".join(parts[len(parts)//2:])
                
                # Remove trailing digits from both parts
                first_clean = re.sub(r'\d+$', '', first_part)
                second_clean = re.sub(r'\d+$', '', second_part)
                
                if first_clean == second_clean:
                    return first_clean
                
                # If they don't match, try removing last digit from each part
                # and see if the base matches
                for i in range(len(parts)):
                    candidate = "_".join(parts[:i+1])
                    candidate_clean = re.sub(r'\d+$', '', candidate)
                    # Check if this candidate appears in canonical samples
                    if candidate_clean:
                        return candidate_clean
        
        return None
    
    def find_canonical_name(potential_sample):
        """Find canonical sample name from potential sample name"""
        # Exact match first
        if potential_sample in canonical_samples:
            return potential_sample
        
        # Try prefix matching: if potential_sample is a prefix of a canonical sample
        # This handles cases like "P6-1_test" -> "P6-1_test1"
        best_match = None
        best_length = 0
        for canon_sample in canonical_samples:
            if canon_sample.startswith(potential_sample):
                if len(canon_sample) > best_length:
                    best_match = canon_sample
                    best_length = len(canon_sample)
        
        if best_match:
            return best_match
        
        # Try reverse: if canonical starts with potential (shouldn't happen but be safe)
        for canon_sample in canonical_samples:
            if potential_sample.startswith(canon_sample):
                return canon_sample
        
        return None
    
    # Build rename mapping for all possible log-derived names
    rename_count = 0
    for log_pattern, canon_name in sample_mapping.items():
        config.general_stats_rename[log_pattern] = canon_name
        rename_count += 1
    
    # Also add pattern-based mappings for names we might encounter
    # We'll handle these in after_modules hook since we need to see actual names first
    log.info(f"Set up {rename_count} sample name mappings in general_stats_rename config")


def _load_sample_mapping():
    """
    Load sample sheet and create a mapping from log-derived names to canonical sample names.
    Returns a dict mapping log names -> canonical sample names.
    """
    sample_mapping = {}
    
    # Try multiple ways to get the annotation path
    annotation_path = None
    
    # Method 1: Direct attribute
    if hasattr(config, "annotation"):
        annotation_path = config.annotation
        log.info(f"Found annotation via attribute: {annotation_path}")
    
    # Method 2: Check if config has a dict-like structure
    if annotation_path is None and hasattr(config, "__dict__"):
        annotation_path = config.__dict__.get("annotation")
        if annotation_path:
            log.info(f"Found annotation via __dict__: {annotation_path}")
    
    # Method 3: Check if config has get method
    if annotation_path is None and hasattr(config, "get"):
        annotation_path = config.get("annotation")
        if annotation_path:
            log.info(f"Found annotation via get(): {annotation_path}")
    
    # Debug: log what we're getting from config
    log.info(f"DEBUG: Final annotation_path = {annotation_path}")
    log.info(f"DEBUG: config type = {type(config)}")
    log.info(f"DEBUG: hasattr(config, 'annotation') = {hasattr(config, 'annotation')}")
    
    if annotation_path is None:
        log.warning("No annotation file found for sample name normalization")
        log.warning("Make sure to pass --cl-config \"annotation: '/path/to/samplesheet.csv'\" to MultiQC")
        return sample_mapping
    
    if not os.path.exists(annotation_path):
        log.warning(f"Annotation file not found: {annotation_path}")
        return sample_mapping
    
    try:
        with open(annotation_path, "r") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sample_name = row.get("sample_name", "").strip() or row.get("sample", "").strip()
                if not sample_name:
                    continue
                
                run = row.get("run", "").strip() or "1"
                
                # Map various log file name patterns to canonical sample name
                # Pattern: bowtie2_{sample}.{run}.err -> {sample}
                sample_mapping[f"bowtie2_{sample_name}.{run}"] = sample_name
                sample_mapping[f"bowtie2_{sample_name}.{run}.err"] = sample_name
                
                # Pattern: bwa_mem2_{sample}.{run}.err -> {sample}
                sample_mapping[f"bwa_mem2_{sample_name}.{run}"] = sample_name
                sample_mapping[f"bwa_mem2_{sample_name}.{run}.err"] = sample_name
                
                # Pattern: sambamba_markdup_{sample} -> {sample}
                sample_mapping[f"sambamba_markdup_{sample_name}"] = sample_name
                sample_mapping[f"sambamba_markdup_{sample_name}.log"] = sample_name
                
                # Pattern: sambamba_sort_{sample}.{run} -> {sample}
                sample_mapping[f"sambamba_sort_{sample_name}.{run}"] = sample_name
                sample_mapping[f"sambamba_sort_{sample_name}.{run}.log"] = sample_name
                
                # Pattern: samtools_merge_{sample} -> {sample}
                sample_mapping[f"samtools_merge_{sample_name}"] = sample_name
                sample_mapping[f"samtools_merge_{sample_name}.log"] = sample_name
                
                # Pattern: samtools_index_{sample} -> {sample}
                sample_mapping[f"samtools_index_{sample_name}"] = sample_name
                sample_mapping[f"samtools_index_{sample_name}.log"] = sample_name
                
                # Also try patterns without run number for cases where run is omitted
                sample_mapping[f"bowtie2_{sample_name}"] = sample_name
                sample_mapping[f"bwa_mem2_{sample_name}"] = sample_name
                sample_mapping[f"sambamba_sort_{sample_name}"] = sample_name
        
        # After loading all samples, create mappings for fastp concatenated patterns
        # Group samples by their base name (without run number)
        samples_by_base = {}
        for row_sample_name in set(sample_mapping.values()):
            # Remove trailing digits to get base name
            base_name = re.sub(r'\d+$', '', row_sample_name)
            if base_name not in samples_by_base:
                samples_by_base[base_name] = []
            samples_by_base[base_name].append(row_sample_name)
        
        # Create concatenated patterns for fastp
        # fastp concatenates samples like: sample1_sample2 when merging runs
        for base_name, samples in samples_by_base.items():
            if len(samples) >= 2:
                # Sort samples to ensure consistent ordering
                sorted_samples = sorted(samples)
                # Create concatenated pattern
                concatenated = "_".join(sorted_samples)
                sample_mapping[concatenated] = base_name
                log.debug(f"Added fastp concatenated pattern: {concatenated} -> {base_name}")
                
        log.info(f"Loaded sample name mapping for {len(set(sample_mapping.values()))} samples")
    except Exception as exc:
        log.warning(f"Could not load annotation file for sample name normalization: {exc}")
    
    return sample_mapping


def cuttag_report_after_modules():
    """
    Hook that runs after all modules are initialized.
    Normalize sample names in general stats to merge log-derived names with canonical names.
    """
    # Halt execution if we've disabled the plugin
    if config.kwargs.get("disable_cuttag_report", False):
        return None
    
    try:
        from multiqc.utils import report
        
        # Load sample name mapping
        sample_mapping = _load_sample_mapping()
        if not sample_mapping:
            log.warning("No sample mapping loaded for normalization")
            return
        
        log.info(f"Loaded sample mapping with {len(sample_mapping)} entries")
        
        # Get the general stats data
        if not hasattr(report, "general_stats_data") or not report.general_stats_data:
            log.warning("No general_stats_data found in report")
            return
        
        # Build reverse mapping: log_name -> canonical_name
        # Also build a mapping from extracted sample names to canonical names
        reverse_mapping = {}
        canonical_samples = set(sample_mapping.values())
        
        # Pattern matching functions
        def extract_sample_from_log_name(log_name):
            """Extract sample name from various log file name patterns"""
            patterns = [
                (r"^bowtie2_(.+?)\.\d+\.err?$", 1),
                (r"^bowtie2_(.+?)\.err?$", 1),  # Without run number
                (r"^bwa_mem2_(.+?)\.\d+\.err?$", 1),
                (r"^bwa_mem2_(.+?)\.err?$", 1),  # Without run number
                (r"^sambamba_markdup_(.+?)(?:\.log)?$", 1),
                (r"^sambamba_sort_(.+?)\.\d+(?:\.log)?$", 1),
                (r"^samtools_merge_(.+?)(?:\.log)?$", 1),
                (r"^samtools_index_(.+?)(?:\.log)?$", 1),
            ]
            
            for pattern, group_idx in patterns:
                match = re.match(pattern, log_name)
                if match:
                    return match.group(group_idx)
            
            # Handle fastp concatenated sample names like "P6-1_test1_P6-1_test2"
            # This happens when fastp merges multiple runs
            # Pattern: {sample}{run}_{sample}{run} -> extract {sample}
            if "_" in log_name and not any(log_name.startswith(prefix) for prefix in ["bowtie2_", "bwa_mem2_", "sambamba_", "samtools_"]):
                # Try to extract common prefix from concatenated pattern
                parts = log_name.split("_")
                if len(parts) >= 2:
                    # Check if it looks like concatenated samples: sample1_sample2
                    # Try to find the longest common prefix
                    first_part = "_".join(parts[:len(parts)//2])
                    second_part = "_".join(parts[len(parts)//2:])
                    
                    # Remove trailing digits from both parts
                    first_clean = re.sub(r'\d+$', '', first_part)
                    second_clean = re.sub(r'\d+$', '', second_part)
                    
                    if first_clean == second_clean:
                        return first_clean
                    
                    # If they don't match, try removing last digit from each part
                    # and see if the base matches
                    for i in range(len(parts)):
                        candidate = "_".join(parts[:i+1])
                        candidate_clean = re.sub(r'\d+$', '', candidate)
                        # Check if this candidate appears in canonical samples
                        if candidate_clean:
                            return candidate_clean
            
            return None
        
        def find_canonical_name(potential_sample):
            """Find canonical sample name from potential sample name"""
            # Exact match first
            if potential_sample in canonical_samples:
                return potential_sample
            
            # Try prefix matching: if potential_sample is a prefix of a canonical sample
            for canon_sample in canonical_samples:
                if canon_sample.startswith(potential_sample) or potential_sample.startswith(canon_sample):
                    # Prefer the longer match (more specific)
                    if len(canon_sample) >= len(potential_sample):
                        return canon_sample
            
            # Try substring matching as last resort
            for canon_sample in canonical_samples:
                if potential_sample in canon_sample or canon_sample in potential_sample:
                    return canon_sample
            
            return None
        
        # Log all current sample names for debugging
        all_sample_names = set()
        
        # Handle both list and dict formats of general_stats_data
        if isinstance(report.general_stats_data, list):
            # Newer MultiQC versions use a list of dicts
            for module_dict in report.general_stats_data:
                for module_name, module_data in module_dict.items():
                    all_sample_names.update(module_data.keys())
        else:
            # Older versions use a dict
            for module_name, module_data in report.general_stats_data.items():
                all_sample_names.update(module_data.keys())
        
        log.info(f"Found {len(all_sample_names)} unique sample names in general stats: {sorted(all_sample_names)[:10]}...")
        log.info(f"Canonical samples from mapping: {sorted(canonical_samples)}")
        
        # Rename samples in general stats
        # Note: general_stats_rename should already be set in execution_start hook
        if not hasattr(config, "general_stats_rename"):
            config.general_stats_rename = {}
        
        renamed_count = 0
        rename_map = {}  # Track all renames to avoid conflicts
        
        # Handle both list and dict formats
        if isinstance(report.general_stats_data, list):
            # Newer MultiQC versions use a list of dicts
            for module_dict in report.general_stats_data:
                for module_name, module_data in list(module_dict.items()):
                    for old_name in list(module_data.keys()):
                        new_name = None
                        
                        # Try exact match first
                        if old_name in sample_mapping:
                            new_name = sample_mapping[old_name]
                            log.debug(f"Exact match: '{old_name}' -> '{new_name}'")
                        else:
                            # Try pattern matching
                            extracted = extract_sample_from_log_name(old_name)
                            if extracted:
                                log.debug(f"Extracted '{extracted}' from '{old_name}'")
                                new_name = find_canonical_name(extracted)
                                if new_name:
                                    log.debug(f"Found canonical name '{new_name}' for extracted '{extracted}'")
                        
                        if new_name and new_name != old_name:
                            # Check if we've already renamed this sample in another module
                            if old_name in rename_map:
                                new_name = rename_map[old_name]
                            
                            # Set MultiQC's general_stats_rename config (recommended approach)
                            config.general_stats_rename[old_name] = new_name
                            
                            # Also directly modify the data structure (backup approach)
                            if new_name in module_data and new_name != old_name:
                                # Merge data if target already exists - check if it's a dict
                                if isinstance(module_data[new_name], dict) and isinstance(module_data[old_name], dict):
                                    module_data[new_name].update(module_data[old_name])
                                    del module_data[old_name]
                                    log.info(f"Merged '{old_name}' into existing '{new_name}' in {module_name}")
                                elif isinstance(module_data[new_name], (int, float)) and isinstance(module_data[old_name], (int, float)):
                                    # Average numeric values from multiple runs
                                    module_data[new_name] = (module_data[new_name] + module_data[old_name]) / 2
                                    del module_data[old_name]
                                    log.info(f"Averaged '{old_name}' into existing '{new_name}' in {module_name}")
                                else:
                                    # Can't merge - just delete the duplicate and keep the first one
                                    del module_data[old_name]
                                    log.info(f"Removed duplicate '{old_name}' (keeping '{new_name}') in {module_name}")
                            else:
                                module_data[new_name] = module_data.pop(old_name)
                                rename_map[old_name] = new_name
                                renamed_count += 1
                                log.info(f"Renamed '{old_name}' -> '{new_name}' in {module_name}")
        else:
            # Older MultiQC versions use a dict
            for module_name, module_data in list(report.general_stats_data.items()):
                for old_name in list(module_data.keys()):
                    new_name = None
                    
                    # Try exact match first
                    if old_name in sample_mapping:
                        new_name = sample_mapping[old_name]
                        log.debug(f"Exact match: '{old_name}' -> '{new_name}'")
                    else:
                        # Try pattern matching
                        extracted = extract_sample_from_log_name(old_name)
                        if extracted:
                            log.debug(f"Extracted '{extracted}' from '{old_name}'")
                            new_name = find_canonical_name(extracted)
                            if new_name:
                                log.debug(f"Found canonical name '{new_name}' for extracted '{extracted}'")
                    
                    if new_name and new_name != old_name:
                        # Check if we've already renamed this sample in another module
                        if old_name in rename_map:
                            new_name = rename_map[old_name]
                        
                        # Set MultiQC's general_stats_rename config (recommended approach)
                        config.general_stats_rename[old_name] = new_name
                        
                        # Also directly modify the data structure (backup approach)
                        if new_name in module_data and new_name != old_name:
                            # Merge data if target already exists
                            if isinstance(module_data[new_name], dict) and isinstance(module_data[old_name], dict):
                                module_data[new_name].update(module_data[old_name])
                                del module_data[old_name]
                                log.info(f"Merged '{old_name}' into existing '{new_name}' in {module_name}")
                            elif isinstance(module_data[new_name], (int, float)) and isinstance(module_data[old_name], (int, float)):
                                # Average numeric values from multiple runs
                                module_data[new_name] = (module_data[new_name] + module_data[old_name]) / 2
                                del module_data[old_name]
                                log.info(f"Averaged '{old_name}' into existing '{new_name}' in {module_name}")
                            else:
                                # Can't merge - just delete the duplicate and keep the first one
                                del module_data[old_name]
                                log.info(f"Removed duplicate '{old_name}' (keeping '{new_name}') in {module_name}")
                        else:
                            module_data[new_name] = module_data.pop(old_name)
                            rename_map[old_name] = new_name
                            renamed_count += 1
                            log.info(f"Renamed '{old_name}' -> '{new_name}' in {module_name}")
        
        if renamed_count > 0:
            log.info(f"Normalized {renamed_count} sample names in MultiQC general stats")
            log.info(f"Set general_stats_rename config with {len(config.general_stats_rename)} entries")
        else:
            log.warning("No sample names were normalized. Check sample mapping and log file patterns.")
            log.debug(f"Sample mapping keys (first 10): {list(sample_mapping.keys())[:10]}")
            log.debug(f"Sample names in report (first 10): {sorted(all_sample_names)[:10]}")
    
    except Exception as exc:
        import traceback
        log.error(f"Error during sample name normalization: {exc}")
        log.debug(traceback.format_exc())


def cuttag_report_before_report_generation():
    """
    Hook that runs just before report generation.
    This is our last chance to normalize sample names before they're written to files.
    """
    # Halt execution if we've disabled the plugin
    if config.kwargs.get("disable_cuttag_report", False):
        return None
    
    try:
        from multiqc.utils import report
        
        # Load sample name mapping
        sample_mapping = _load_sample_mapping()
        if not sample_mapping:
            log.debug("No sample mapping for before_report_generation hook")
            return
        
        # Apply sample name normalization one more time to ensure it's in the output
        if hasattr(report, "general_stats_data") and report.general_stats_data:
            canonical_samples = set(sample_mapping.values())
            
            # Pattern matching function (same as before)
            def extract_sample_from_log_name(log_name):
                """Extract sample name from various log file name patterns"""
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
                        return match.group(group_idx)
                
                # Handle fastp concatenated sample names
                if "_" in log_name and not any(log_name.startswith(prefix) for prefix in ["bowtie2_", "bwa_mem2_", "sambamba_", "samtools_"]):
                    parts = log_name.split("_")
                    if len(parts) >= 2:
                        first_part = "_".join(parts[:len(parts)//2])
                        second_part = "_".join(parts[len(parts)//2:])
                        first_clean = re.sub(r'\d+$', '', first_part)
                        second_clean = re.sub(r'\d+$', '', second_part)
                        if first_clean == second_clean:
                            return first_clean
                        for i in range(len(parts)):
                            candidate = "_".join(parts[:i+1])
                            candidate_clean = re.sub(r'\d+$', '', candidate)
                            if candidate_clean:
                                return candidate_clean
                return None
            
            def find_canonical_name(potential_sample):
                """Find canonical sample name from potential sample name"""
                if potential_sample in canonical_samples:
                    return potential_sample
                for canon_sample in canonical_samples:
                    if canon_sample.startswith(potential_sample) or potential_sample.startswith(canon_sample):
                        if len(canon_sample) >= len(potential_sample):
                            return canon_sample
                for canon_sample in canonical_samples:
                    if potential_sample in canon_sample or canon_sample in potential_sample:
                        return canon_sample
                return None
            
            # Process all samples and build a comprehensive rename map
            final_rename_map = {}
            
            # Debug: log all current sample names
            all_current_samples = set()
            if isinstance(report.general_stats_data, list):
                for module_dict in report.general_stats_data:
                    for module_name, module_data in module_dict.items():
                        all_current_samples.update(module_data.keys())
            else:
                for module_name, module_data in report.general_stats_data.items():
                    all_current_samples.update(module_data.keys())
            log.info(f"Current samples in general_stats_data before rename: {sorted(all_current_samples)}")
            
            # Check for data_sources or saved_raw_data
            if hasattr(report, "data_sources"):
                log.info(f"Found data_sources attribute - type: {type(report.data_sources)}")
                if isinstance(report.data_sources, dict):
                    log.info(f"data_sources keys: {list(report.data_sources.keys())[:10]}")
                    # data_sources is a dict of module_name -> list of sample names or dict
                    # We need to rename the sample names in each module's data
                    for module_name in list(report.data_sources.keys()):
                        module_data = report.data_sources[module_name]
                        log.info(f"data_sources[{module_name}] type: {type(module_data)}, value: {module_data if not isinstance(module_data, dict) else f'dict with {len(module_data)} keys'}")
                        if isinstance(module_data, dict):
                            log.info(f"data_sources[{module_name}] sample keys: {list(module_data.keys())}")
                        
                        if isinstance(module_data, list):
                            # Get the list of samples for this module
                            old_samples = module_data[:]
                            new_samples = []
                            for sample in old_samples:
                                # Check if this sample needs renaming
                                renamed = False
                                for old_name, new_name in final_rename_map.items():
                                    if sample == old_name:
                                        if new_name not in new_samples:
                                            new_samples.append(new_name)
                                        renamed = True
                                        log.info(f"Renamed in data_sources[{module_name}]: {old_name} -> {new_name}")
                                        break
                                if not renamed:
                                    new_samples.append(sample)
                            # Update the list
                            report.data_sources[module_name] = new_samples
                            log.info(f"Updated data_sources[{module_name}]: {old_samples} -> {new_samples}")
                        elif isinstance(module_data, dict):
                            # Rename keys in the dict
                            for old_name, new_name in final_rename_map.items():
                                if old_name in module_data:
                                    if new_name in module_data and new_name != old_name:
                                        # Merge if both exist
                                        if isinstance(module_data[new_name], list) and isinstance(module_data[old_name], list):
                                            module_data[new_name].extend(module_data[old_name])
                                        module_data[new_name] = module_data.pop(old_name)
                                    else:
                                        module_data[new_name] = module_data.pop(old_name)
                                    log.info(f"Renamed key in data_sources[{module_name}]: {old_name} -> {new_name}")
            
            if hasattr(report, "saved_raw_data"):
                log.info(f"Found saved_raw_data attribute - type: {type(report.saved_raw_data)}")
                if isinstance(report.saved_raw_data, dict):
                    log.info(f"saved_raw_data keys (first 10): {list(report.saved_raw_data.keys())[:10]}")
                    # Rename samples in saved_raw_data
                    for module_name in list(report.saved_raw_data.keys()):
                        if isinstance(report.saved_raw_data[module_name], dict):
                            for old_name, new_name in final_rename_map.items():
                                if old_name in report.saved_raw_data[module_name]:
                                    report.saved_raw_data[module_name][new_name] = report.saved_raw_data[module_name].pop(old_name)
                                    log.info(f"Renamed in saved_raw_data[{module_name}]: {old_name} -> {new_name}")
            
            if isinstance(report.general_stats_data, list):
                for module_dict in report.general_stats_data:
                    for module_name, module_data in module_dict.items():
                        for sample_name in list(module_data.keys()):
                            if sample_name in sample_mapping:
                                final_rename_map[sample_name] = sample_mapping[sample_name]
                            else:
                                extracted = extract_sample_from_log_name(sample_name)
                                if extracted:
                                    canonical = find_canonical_name(extracted)
                                    if canonical:
                                        final_rename_map[sample_name] = canonical
            else:
                for module_name, module_data in report.general_stats_data.items():
                    for sample_name in list(module_data.keys()):
                        if sample_name in sample_mapping:
                            final_rename_map[sample_name] = sample_mapping[sample_name]
                        else:
                            extracted = extract_sample_from_log_name(sample_name)
                            if extracted:
                                canonical = find_canonical_name(extracted)
                                if canonical:
                                    final_rename_map[sample_name] = canonical
            
            # Apply the rename map to general_stats_rename config
            if not hasattr(config, "general_stats_rename"):
                config.general_stats_rename = {}
            config.general_stats_rename.update(final_rename_map)
            
            # Also directly rename in the data structure to ensure it appears in output files
            if isinstance(report.general_stats_data, list):
                for module_dict in report.general_stats_data:
                    for module_name, module_data in list(module_dict.items()):
                        for old_name, new_name in final_rename_map.items():
                            if old_name in module_data and old_name != new_name:
                                if new_name in module_data:
                                    # Merge if target exists
                                    if isinstance(module_data[new_name], dict) and isinstance(module_data[old_name], dict):
                                        module_data[new_name].update(module_data[old_name])
                                    elif isinstance(module_data[new_name], (int, float)) and isinstance(module_data[old_name], (int, float)):
                                        module_data[new_name] = (module_data[new_name] + module_data[old_name]) / 2
                                    del module_data[old_name]
                                else:
                                    # Simple rename
                                    module_data[new_name] = module_data.pop(old_name)
            else:
                for module_name, module_data in list(report.general_stats_data.items()):
                    for old_name, new_name in final_rename_map.items():
                        if old_name in module_data and old_name != new_name:
                            if new_name in module_data:
                                # Merge if target exists
                                if isinstance(module_data[new_name], dict) and isinstance(module_data[old_name], dict):
                                    module_data[new_name].update(module_data[old_name])
                                elif isinstance(module_data[new_name], (int, float)) and isinstance(module_data[old_name], (int, float)):
                                    module_data[new_name] = (module_data[new_name] + module_data[old_name]) / 2
                                del module_data[old_name]
                            else:
                                # Simple rename
                                module_data[new_name] = module_data.pop(old_name)
            
            log.info(f"Before report generation: Applied {len(final_rename_map)} sample name mappings")
    
    except Exception as exc:
        import traceback
        log.error(f"Error in before_report_generation hook: {exc}")
        log.debug(traceback.format_exc())

