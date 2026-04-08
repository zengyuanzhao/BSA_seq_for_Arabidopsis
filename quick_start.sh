#!/bin/bash
# quick_start.sh - BSA-seq quick launch script
# Use this if the conda environment (BSA_seq) is already set up.

set -e

echo "=========================================="
echo "  BSA-seq Pipeline - Quick Start"
echo "=========================================="
echo ""

# Activate conda environment
if command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate BSA_seq
elif command -v conda &> /dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate BSA_seq
else
    echo "Error: micromamba or conda not found."
    echo "Please run: bash setup_and_run.sh"
    exit 1
fi

# Check configuration file
if [ ! -f "config.sh" ]; then
    echo "Error: config.sh not found."
    echo "Please run: bash setup_and_run.sh"
    exit 1
fi
source config.sh

echo "[OK] Environment activated"
echo ""

# Check input files
echo "Checking input files..."
for f in "$SHORT_R1" "$SHORT_R2" "$LONG_R1" "$LONG_R2"; do
    if [ ! -f "$f" ]; then
        echo "[MISSING] $f"
        exit 1
    fi
done
echo "[OK] All input files present"

# Check reference genome
if [ ! -f "$REF" ]; then
    echo "[MISSING] Reference genome: $REF"
    exit 1
fi
echo "[OK] Reference genome: $REF"
echo ""

# Show current configuration
echo "Analysis configuration:"
echo "  Short-root pool R1: $SHORT_R1"
echo "  Long-root pool  R1: $LONG_R1"
echo "  Target region:      ${FINE_CHR}:${FINE_START}-${FINE_END}"
echo "  Threads:            $THREADS"
echo ""

# Remove old results if requested
if [ -d "$OUTDIR" ]; then
    read -p "Old results directory $OUTDIR/ detected. Delete it? (y/n): " clean
    if [ "$clean" = "y" ] || [ "$clean" = "Y" ]; then
        rm -rf "$OUTDIR"
        echo "[OK] Old results removed"
    fi
fi

# Launch pipeline
echo "Launching pipeline..."
echo "Log file: run.log"
echo "Monitor: tail -f run.log"
echo ""

nohup bash run_pipeline.sh > run.log 2>&1 &
PID=$!
echo "[OK] Pipeline started, PID: $PID"
echo ""
echo "Useful commands:"
echo "  Monitor log:    tail -f run.log"
echo "  Check process:  ps aux | grep $PID"
echo "  View results:   ls $OUTDIR/result/"
echo ""
echo "Estimated runtime: 2-4 hours"
echo ""

# Tail log in foreground (Ctrl+C exits tail but does NOT stop the background pipeline)
sleep 2
tail -f run.log
