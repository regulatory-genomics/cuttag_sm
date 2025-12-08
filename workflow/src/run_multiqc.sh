#!/bin/bash
# Standalone bash script to run MultiQC for CUT&Tag pipeline
# Usage: ./run_multiqc.sh [DATA_DIR] [SAMPLE_SHEET]

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PLUGIN_DIR="$WORKFLOW_DIR/src/multiqc_cuttag"

# Set default values or use command line arguments
DATA_DIR="${1:-${DATA_DIR:-./results}}"
SAMPLE_SHEET="${2:-${SAMPLE_SHEET:-config/samplesheet.csv}}"
LOG_FILE="${LOG_FILE:-${DATA_DIR}/logs/multiqc.log}"

# Remove trailing slash from DATA_DIR
DATA_DIR="${DATA_DIR%/}"

# Check if DATA_DIR exists
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: DATA_DIR '$DATA_DIR' does not exist" >&2
    exit 1
fi

# Check if SAMPLE_SHEET exists
if [ ! -f "$SAMPLE_SHEET" ]; then
    echo "Warning: SAMPLE_SHEET '$SAMPLE_SHEET' not found. Sample name normalization may not work." >&2
fi

# Check if plugin directory exists and ensure it's accessible to MultiQC
if [ ! -d "$PLUGIN_DIR" ]; then
    echo "Warning: Plugin directory '$PLUGIN_DIR' not found. CUT&Tag MultiQC plugin may not be loaded." >&2
else
    # Check if plugin is installed, if not install it in development mode
    if ! python -c "import cuttag_report" 2>/dev/null; then
        echo "Installing CUT&Tag MultiQC plugin in development mode..."
        cd "$PLUGIN_DIR" || exit 1
        pip install -e . > /dev/null 2>&1
        if [ $? -eq 0 ]; then
            echo "Plugin installed successfully."
        else
            echo "Warning: Failed to install plugin. Adding to PYTHONPATH instead..." >&2
            export PYTHONPATH="$PLUGIN_DIR:$PYTHONPATH"
        fi
        cd - > /dev/null || exit 1
    else
        echo "CUT&Tag MultiQC plugin is already installed."
    fi
fi

# Create log directory if it doesn't exist
mkdir -p "$(dirname "$LOG_FILE")"

# Run MultiQC
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# Convert SAMPLE_SHEET to absolute path
SAMPLE_SHEET_ABS="$(realpath "$SAMPLE_SHEET" 2>/dev/null || echo "$SAMPLE_SHEET")"

echo "Running MultiQC with annotation file: $SAMPLE_SHEET_ABS"

multiqc "$DATA_DIR" \
  -f \
  -c workflow/src/multiqc_conf.yml \
  -o "$DATA_DIR/Report/multiqc" \
  --cl-config "annotation: $SAMPLE_SHEET_ABS" \
  --ignore "$DATA_DIR/downtream_res/Homer" \
  --ignore "$DATA_DIR/middle_file" \
  --ignore "$DATA_DIR/Important_processed/Bam" \
  --ignore "$DATA_DIR/Important_processed/Track" \
  --ignore "$DATA_DIR/Important_processed/Peaks" \
  > "$LOG_FILE" 2>&1

# Check exit status
if [ $? -eq 0 ]; then
    echo "MultiQC completed successfully. Report: $DATA_DIR/Report/multiqc/multiqc_report.html"
else
    echo "Error: MultiQC failed. Check log: $LOG_FILE" >&2
    exit 1
fi

