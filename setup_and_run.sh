#!/bin/bash
# =============================================================================
# BSA-seq Analysis Pipeline - One-click Install and Run
# Designed for users without a bioinformatics background
# =============================================================================

set -e

# Colored output helpers
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[OK]${NC} $1"; }
print_warn()    { echo -e "${YELLOW}[WARN]${NC} $1"; }
print_error()   { echo -e "${RED}[ERROR]${NC} $1"; }

# =============================================================================
# Step 1: Detect and install a package manager
# =============================================================================
setup_package_manager() {
    print_info "Checking for a package manager..."

    if command -v micromamba &> /dev/null; then
        print_success "Found micromamba"
        PKG_MGR="micromamba"
        PKG_MGR_BIN="$(command -v micromamba)"
        return 0
    fi
    if command -v mamba &> /dev/null; then
        print_success "Found mamba"
        PKG_MGR="mamba"
        PKG_MGR_BIN="$(command -v mamba)"
        return 0
    fi
    if command -v conda &> /dev/null; then
        print_success "Found conda"
        PKG_MGR="conda"
        PKG_MGR_BIN="$(command -v conda)"
        return 0
    fi

    print_warn "No conda/mamba/micromamba detected. Installing micromamba..."
    mkdir -p "$HOME/.local/bin"

    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | \
            tar -xvj -C "$HOME/.local/bin/" --strip-components=1 bin/micromamba
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        curl -Ls https://micro.mamba.pm/api/micromamba/osx-64/latest | \
            tar -xvj -C "$HOME/.local/bin/" --strip-components=1 bin/micromamba
    else
        print_error "Unsupported OS: $OSTYPE"
        exit 1
    fi

    export PATH="$HOME/.local/bin:$PATH"
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
    PKG_MGR="micromamba"
    PKG_MGR_BIN="$HOME/.local/bin/micromamba"
    print_success "micromamba installed successfully"
}

# =============================================================================
# Step 2: Create and configure the conda environment
# =============================================================================
setup_environment() {
    print_info "Setting up analysis environment..."
    ENV_NAME="BSA_seq"

    if [ "$PKG_MGR" = "micromamba" ]; then
        if [ ! -f "$HOME/.bashrc" ] || ! grep -q "micromamba" "$HOME/.bashrc"; then
            print_info "Initializing micromamba shell..."
            "$PKG_MGR_BIN" shell init -s bash
        fi
        eval "$("$PKG_MGR_BIN" shell hook --shell bash)"

        if ! "$PKG_MGR_BIN" env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
            print_info "Creating environment: $ENV_NAME"
            "$PKG_MGR_BIN" create -n "$ENV_NAME" -y
        fi
        micromamba activate "$ENV_NAME"
        print_info "Installing bioinformatics tools (this may take 10-20 minutes)..."
        micromamba install -c bioconda -c conda-forge -y \
            trimmomatic bwa samtools gatk4 snpeff python pandas tabix bcftools
    else
        if ! "$PKG_MGR_BIN" env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
            print_info "Creating environment: $ENV_NAME"
            "$PKG_MGR_BIN" create -n "$ENV_NAME" -y
        fi
        source "$("$PKG_MGR_BIN" info --base)/etc/profile.d/conda.sh"
        conda activate "$ENV_NAME"
        print_info "Installing bioinformatics tools (this may take 10-20 minutes)..."
        "$PKG_MGR_BIN" install -c bioconda -c conda-forge -y \
            trimmomatic bwa samtools gatk4 snpeff python pandas tabix bcftools
    fi
    print_success "Environment ready"
}

# =============================================================================
# Step 3: Download the snpEff Arabidopsis database
# =============================================================================
setup_snpeff_database() {
    print_info "Checking snpEff Arabidopsis database..."
    _DB_LIST=$(snpEff databases 2>/dev/null || true)
    AVAILABLE=$(echo "$_DB_LIST" | awk '$1 == "Arabidopsis_thaliana" {print $1; exit}' || true)
    if [ -n "$AVAILABLE" ]; then
        print_success "snpEff Arabidopsis_thaliana database already present"
    else
        print_info "Downloading Arabidopsis_thaliana database (~200 MB)..."
        if ! snpEff download Arabidopsis_thaliana; then
            print_error "Failed to download snpEff Arabidopsis_thaliana database."
            exit 1
        fi
        print_success "Database download complete"
    fi
}

# =============================================================================
# Step 4: Data preparation wizard
# =============================================================================
data_setup_wizard() {
    print_info "=========================================="
    print_info " Data Preparation Wizard"
    print_info "=========================================="
    echo ""

    if [ -d "RawData" ]; then
        print_info "Detected RawData directory:"
        ls -1 RawData/
        echo ""
        read -p "Use existing data? (y/n): " use_existing
        if [ "$use_existing" = "y" ] || [ "$use_existing" = "Y" ]; then
            if [ -f "config.sh" ]; then
                source config.sh
                return 0
            fi
            print_warn "config.sh not found, continuing with the setup wizard."
        fi
    fi

    echo "Please prepare your data in the following structure:"
    echo ""
    echo "  RawData/"
    echo "  |-- mut_pool_name/"
    echo "  |   |-- *_1.fq.gz  (R1 file)"
    echo "  |   \`-- *_2.fq.gz  (R2 file)"
    echo "  \`-- wt_pool_name/"
    echo "      |-- *_1.fq.gz  (R1 file)"
    echo "      \`-- *_2.fq.gz  (R2 file)"
    echo ""
    read -p "Press Enter to continue..."

    if [ ! -d "RawData" ]; then
        print_error "RawData directory not found. Please prepare your data first."
        exit 1
    fi

    print_info "Select the mutant pool folder:"
    select mut_dir in $(ls RawData/) "Enter manually"; do
        [ "$mut_dir" = "Enter manually" ] && read -p "Mutant pool folder name: " mut_dir
        break
    done

    print_info "Select the wild-type pool folder:"
    select wt_dir in $(ls RawData/) "Enter manually"; do
        [ "$wt_dir" = "Enter manually" ] && read -p "Wild-type pool folder name: " wt_dir
        break
    done

    mut_r1=$(ls "RawData/$mut_dir/"*_1.fq.gz 2>/dev/null | head -1 || true)
    mut_r2=$(ls "RawData/$mut_dir/"*_2.fq.gz 2>/dev/null | head -1 || true)
    wt_r1=$(ls  "RawData/$wt_dir/"*_1.fq.gz  2>/dev/null | head -1 || true)
    wt_r2=$(ls  "RawData/$wt_dir/"*_2.fq.gz  2>/dev/null | head -1 || true)

    for pair_name in "mut_r1:RawData/$mut_dir/*_1.fq.gz" \
                     "mut_r2:RawData/$mut_dir/*_2.fq.gz" \
                     "wt_r1:RawData/$wt_dir/*_1.fq.gz" \
                     "wt_r2:RawData/$wt_dir/*_2.fq.gz"; do
        var_name="${pair_name%%:*}"
        pattern="${pair_name##*:}"
        val="${!var_name}"
        if [[ -z "$val" ]]; then
            print_error "No file matching '$pattern' was found."
            print_error "Please check your RawData directory structure and file naming."
            exit 1
        fi
    done

    update_config "$mut_r1" "$mut_r2" "$wt_r1" "$wt_r2"
}

# =============================================================================
# Step 5: Write configuration file
# =============================================================================
update_config() {
    local mut_r1=$1 mut_r2=$2 wt_r1=$3 wt_r2=$4
    print_info "Writing configuration file..."

    echo ""
    echo "Please set the fine-mapping interval (BSA target region):"
    read -p "Chromosome (default: 1): "          chr_input;   chr=${chr_input:-1}
    read -p "Start position (default: 10000000): " start_input; start=${start_input:-10000000}
    read -p "End position   (default: 15000000): " end_input;   end=${end_input:-15000000}

    cat > config.sh << EOF
# ===== User Configuration =====
REF="./reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
KNOWN_VCF=""
MUT_R1="$mut_r1"
MUT_R2="$mut_r2"
WT_R1="$wt_r1"
WT_R2="$wt_r2"
FINE_CHR="$chr"
FINE_START=$start
FINE_END=$end
THREADS=8
OUTDIR="./bsa_output"
EOF

    print_success "config.sh written successfully"
}

# =============================================================================
# Step 6: Check reference genome
# =============================================================================
check_reference() {
    print_info "Checking reference genome..."
    mkdir -p reference
    REF_FILE="reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

    if [ ! -f "$REF_FILE" ]; then
        print_warn "Reference genome not found: $REF_FILE"
        print_info "Auto-downloading from Ensembl Plants..."
        wget -q --show-progress \
            ftp://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
            -O reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
        gunzip reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
        print_success "Reference genome downloaded and ready"
    else
        print_success "Reference genome ready"
    fi
}

# =============================================================================
# Step 7: Launch the analysis pipeline
# =============================================================================
run_pipeline() {
    print_info "=========================================="
    print_info " Starting BSA-seq Analysis Pipeline"
    print_info "=========================================="
    echo ""

    if [ -d "bsa_output" ]; then
        read -p "Old results directory bsa_output/ detected. Delete it? (y/n): " clean_old
        if [ "$clean_old" = "y" ] || [ "$clean_old" = "Y" ]; then
            rm -rf bsa_output/
            rm -f run.log
            print_success "Old results removed"
        fi
    fi

    print_info "Launching pipeline in the background..."
    print_info "Log file: run.log"
    print_info "Monitor progress: tail -f run.log"
    echo ""

    nohup bash run_pipeline.sh > run.log 2>&1 &
    PID=$!
    print_success "Pipeline started, PID: $PID"
    echo ""
    echo "View log:     tail -f run.log"
    echo "View results: ls bsa_output/result/"
}

# =============================================================================
# Main
# =============================================================================
main() {
    echo ""
    echo "=========================================="
    echo "  BSA-seq Pipeline - One-click Setup & Run"
    echo "=========================================="
    echo ""

    setup_package_manager
    setup_environment
    setup_snpeff_database
    data_setup_wizard
    check_reference
    run_pipeline

    echo ""
    echo "=========================================="
    print_success "Setup complete!"
    echo "=========================================="
    echo ""
}

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    main "$@"
fi
