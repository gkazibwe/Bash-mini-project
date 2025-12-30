#!/bin/bash

# ==============================================================================
# SCRIPT: comparative_analysis.sh
# DESCRIPTION: Comparative variant analysis pipeline between two microorganism strains.
#              Includes Download, QC, Alignment, Variant Calling, Annotation, and
#              shared/unique missense SNP analysis.
#
# AUTHORS:
#   1. KAZIBWE GEORGE 2025/HD07/25965U <gkazibwe@gmail.com>
#   2. WASSWA CHARLES LWANGA <wasswacharleslwanga4@gmail.com>
#   3. DDUMBA FRANCIS SEMAKUBA <frncsddumbasema@gmail.com>
#
# REQUIRES: conda activate variant_analysis_env
# ==============================================================================


# --- 1. CONFIGURATION VARIABLES ---
# (Change these variables to analyze different strains/organisms)

# Sample Accessions (SRA Run IDs)
STRAIN_1_ID="SRR36143512" # Herpes simplex Strain A
STRAIN_2_ID="SRR23265797" # Herpes simplex Strain B
SAMPLES=($STRAIN_1_ID $STRAIN_2_ID)

# Reference Genome (Human herpesvirus 1 strain 17)
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/859/985/GCF_000859985.2_ViralProj15217/GCF_000859985.2_ViralProj15217_genomic.fna.gz"
REF_NAME="Human_herpesvirus_1_strain_17"
SNPEFF_DB="HSV1_strain17" # Database name for snpEff
SNPEFF_CONFIG="$CONDA_PREFIX/share/snpeff-*/snpEff.config"

# CPU Threads
THREADS=8

# --- 2. DIRECTORY SETUP ---
BASE_DIR=$(pwd)
DIR_REF="$BASE_DIR/1_reference"
DIR_RAW="$BASE_DIR/2_raw_reads"
DIR_QC="$BASE_DIR/3_qc_reports"
DIR_ALIGN="$BASE_DIR/4_alignment"
DIR_VARIANTS="$BASE_DIR/5_variants"
DIR_RESULTS="$BASE_DIR/6_final_results"
CHECKPOINT_FILE="$BASE_DIR/.checkpoint"

# Create Directories
mkdir -p $DIR_REF $DIR_RAW $DIR_QC $DIR_ALIGN $DIR_VARIANTS $DIR_RESULTS

# --- 3. ERROR HANDLING & CHECKPOINT FUNCTIONS ---

# Exit immediately if a command exits with a non-zero status
set -e
set -o pipefail

log() {
    echo -e "\n[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

check_step() {
    local step_name=$1
    if grep -q "$step_name" "$CHECKPOINT_FILE" 2>/dev/null; then
        log "Step '$step_name' already completed. Skipping."
        return 0 # True, step exists
    else
        return 1 # False, step needs running
    fi
}

mark_step() {
    echo "$1" >> "$CHECKPOINT_FILE"
    log "Step '$1' completed successfully."
}

log "Starting Comparative Analysis Pipeline..."
log "Comparing $STRAIN_1_ID vs $STRAIN_2_ID"

# ==============================================================================
# STEP 1: DOWNLOAD REFERENCE GENOME
# ==============================================================================
if ! check_step "DOWNLOAD_REF"; then
    log "Downloading Reference Genome..."
    wget -O "$DIR_REF/${REF_NAME}.fna.gz" "$REF_URL"
    gunzip -f "$DIR_REF/${REF_NAME}.fna.gz"
    mark_step "DOWNLOAD_REF"
fi

# ==============================================================================
# STEP 1b: DOWNLOAD GENE ANNOTATION (GFF)
# ==============================================================================
if ! check_step "DOWNLOAD_GFF"; then
    log "Downloading HSV-1 gene annotation (GFF)..."
    wget -O "$DIR_REF/${REF_NAME}.gff.gz" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/859/985/GCF_000859985.2_ViralProj15217/GCF_000859985.2_ViralProj15217_genomic.gff.gz"
    gunzip -f "$DIR_REF/${REF_NAME}.gff.gz"
    mark_step "DOWNLOAD_GFF"
fi

# ==============================================================================
# STEP 2: DOWNLOAD RAW READS (SRA)
# ==============================================================================
if ! check_step "DOWNLOAD_READS"; then
    log "Downloading Raw Reads from NCBI SRA..."
    for sample in "${SAMPLES[@]}"; do
        log "Downloading $sample..."
        # fastq-dump splits files into forward/reverse (_1.fastq, _2.fastq)
        fasterq-dump "$sample" -O "$DIR_RAW" -e $THREADS
        gzip "$DIR_RAW/${sample}_1.fastq"
        gzip "$DIR_RAW/${sample}_2.fastq"

    done
    mark_step "DOWNLOAD_READS"
fi

# ==============================================================================
# STEP 2b: REPAIR PAIRED-END READS
# ==============================================================================
if ! check_step "REPAIR_READS"; then
    log "Repairing paired-end FASTQ files..."

    for sample in "${SAMPLES[@]}"; do
        repair.sh \
            in1="$DIR_RAW/${sample}_1.fastq.gz" \
            in2="$DIR_RAW/${sample}_2.fastq.gz" \
            out1="$DIR_RAW/${sample}_1.repaired.fastq.gz" \
            out2="$DIR_RAW/${sample}_2.repaired.fastq.gz" \
            outs="$DIR_RAW/${sample}_singletons.fastq.gz"
    done

    mark_step "REPAIR_READS"
fi

# ==============================================================================
# STEP 3: QUALITY CONTROL (FastQC + MultiQC)
# ==============================================================================
if ! check_step "RUN_QC"; then
    log "Running FastQC on raw reads..."
    fastqc -t $THREADS "$DIR_RAW"/*.fastq.gz -o "$DIR_QC"
    
    log "Aggregating reports with MultiQC..."
    multiqc "$DIR_QC" -o "$DIR_QC" -n "Combined_QC_Report"
    mark_step "RUN_QC"
fi

# ==============================================================================
# STEP 4: INDEX REFERENCE & ALIGNMENT
# ==============================================================================
if ! check_step "ALIGNMENT"; then
    log "Indexing Reference Genome..."
    bwa index "$DIR_REF/${REF_NAME}.fna"
    
    for sample in "${SAMPLES[@]}"; do
        log "Aligning $sample to reference..."
        
        R1="$DIR_RAW/${sample}_1.repaired.fastq.gz"
        R2="$DIR_RAW/${sample}_2.repaired.fastq.gz"

        
        # Helper variables for Read Group (RG) needed for variant calling
        # ID=Sample, SM=SampleName, PL=Illumina
        RG_TAG="@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA"
        
        # BWA MEM -> Samtools View (to BAM) -> Samtools Sort -> Write File
        bwa mem -t $THREADS -R "$RG_TAG" "$DIR_REF/${REF_NAME}.fna" "$R1" "$R2" | \
        samtools view -bS - | \
        samtools sort -@ $THREADS -o "$DIR_ALIGN/${sample}.sorted.bam"
        
        log "Indexing BAM file for $sample..."
        samtools index "$DIR_ALIGN/${sample}.sorted.bam"
    done
    mark_step "ALIGNMENT"
fi

# ==============================================================================
# STEP 5: VARIANT CALLING (Joint Calling)
# ==============================================================================
if ! check_step "VARIANT_CALLING"; then
    log "Calling variants (Joint Calling)..."
    
    # mpileup: generates genotype likelihoods
    # call: performs actual SNP calling
    # -m: multiallelic-caller (standard)
    # -v: output variants only
    
    bcftools mpileup -f "$DIR_REF/${REF_NAME}.fna" \
        "$DIR_ALIGN/$STRAIN_1_ID.sorted.bam" \
        "$DIR_ALIGN/$STRAIN_2_ID.sorted.bam" | \
    bcftools call -mv -Ob -o "$DIR_VARIANTS/joint_variants.bcf"
    
    mark_step "VARIANT_CALLING"
fi

# ==============================================================================
# STEP 5b: BUILD CUSTOM SNPEFF DATABASE (HSV-1)
# ==============================================================================
if ! check_step "BUILD_SNPEFF_DB"; then
    log "Building custom SnpEff database for HSV-1..."

    # 1. Locate the default config in Conda safely
    DEFAULT_CONFIG=$(find "$CONDA_PREFIX" -name snpEff.config | head -n 1)
    
    if [ -z "$DEFAULT_CONFIG" ]; then
        log "Error: Could not find snpEff.config in $CONDA_PREFIX"
        exit 1
    fi

    # 2. Copy config to current directory (Local Config)
    cp "$DEFAULT_CONFIG" ./snpEff.config
    
    # 3. Configure 'data.dir' to point to a local folder './data'
    sed -i 's|^data.dir = .*|data.dir = ./data|' ./snpEff.config

    # 4. Add the Genome Entry (only if not already there to prevent duplicates)
    if ! grep -q "$SNPEFF_DB.genome" ./snpEff.config; then
        echo -e "\n# Custom Genome: HSV-1" >> ./snpEff.config
        echo -e "$SNPEFF_DB.genome : Human herpes simplex virus 1 (strain 17)" >> ./snpEff.config
    fi

    # 5. Prepare Data Directories
    mkdir -p "./data/$SNPEFF_DB"

    # Copy Reference and GFF to the expected location
    cp "$DIR_REF/${REF_NAME}.fna" "./data/$SNPEFF_DB/sequences.fa"
    cp "$DIR_REF/${REF_NAME}.gff" "./data/$SNPEFF_DB/genes.gff"

    # 6. Build Database using the LOCAL config
    # ADDED FLAGS: -noCheckCds -noCheckProtein (Fixes the "File not found" error)
    snpEff build -c ./snpEff.config -gff3 -noCheckCds -noCheckProtein -v "$SNPEFF_DB"

    mark_step "BUILD_SNPEFF_DB"
fi

# ==============================================================================
# STEP 6: VARIANT ANNOTATION (SnpEff)
# ==============================================================================
if ! check_step "ANNOTATION"; then
    log "Annotating Variants..."
    
    # Note the addition of '-c ./snpEff.config'
    bcftools view "$DIR_VARIANTS/joint_variants.bcf" | \
    snpEff eff -c ./snpEff.config -v "$SNPEFF_DB" - > "$DIR_VARIANTS/annotated_variants.vcf"
    
    # Compress and index the VCF for bcftools operations
    bgzip -f "$DIR_VARIANTS/annotated_variants.vcf"
    tabix -p vcf "$DIR_VARIANTS/annotated_variants.vcf.gz"
    
    mark_step "ANNOTATION"
fi

# ==============================================================================
# STEP 7: IDENTIFY UNIQUE/SHARED MISSENSE VARIANTS
# ==============================================================================
if ! check_step "ANALYSIS"; then
    log "Filtering for Missense Variants and Comparing Strains..."
    
    # 1. Locate SnpSift.jar dynamically (Searching both uppercase and lowercase variants)
    # We search in $CONDA_PREFIX/share to find the jar within the snpeff installation
    SNPSIFT_JAR=$(find "$CONDA_PREFIX/share" -iname "SnpSift.jar" | head -n 1)

    if [ -z "$SNPSIFT_JAR" ]; then
        log "Error: SnpSift.jar not found. Attempting to locate via snpSift wrapper..."
        # Some conda versions use a wrapper script; we try to find the jar relative to it
        SNPSIFT_PATH=$(which snpSift 2>/dev/null || which SnpSift 2>/dev/null)
        if [ -n "$SNPSIFT_PATH" ]; then
             # Extract directory and look for the jar
             SNPSIFT_JAR=$(dirname $(readlink -f "$SNPSIFT_PATH"))/SnpSift.jar
        fi
    fi

    if [ ! -f "$SNPSIFT_JAR" ]; then
        log "Critical Error: Could not find SnpSift.jar. Please ensure 'snpeff' is installed."
        exit 1
    fi
    
    INPUT_VCF="$DIR_VARIANTS/annotated_variants.vcf.gz"
    
    # 2. Extract only Missense variants using SnpSift filter
    log "Filtering missense variants..."
    if [ -n "$SNPSIFT_JAR" ]; then
        zcat "$INPUT_VCF" | java -jar "$SNPSIFT_JAR" filter \
            "ANN[*].EFFECT has 'missense_variant'" \
            > "$DIR_RESULTS/all_missense.vcf"
    else
        zcat "$INPUT_VCF" | $SNPSIFT_CMD filter \
            "ANN[*].EFFECT has 'missense_variant'" \
            > "$DIR_RESULTS/all_missense.vcf"
    fi
    
    bgzip -f "$DIR_RESULTS/all_missense.vcf"
    tabix -p vcf "$DIR_RESULTS/all_missense.vcf.gz"
    
    # 3. Split the joint VCF into individual VCFs
    log "Splitting VCF by sample..."
    bcftools view -s "$STRAIN_1_ID" "$DIR_RESULTS/all_missense.vcf.gz" -Oz -o "$DIR_RESULTS/${STRAIN_1_ID}_missense.vcf.gz"
    bcftools view -s "$STRAIN_2_ID" "$DIR_RESULTS/all_missense.vcf.gz" -Oz -o "$DIR_RESULTS/${STRAIN_2_ID}_missense.vcf.gz"
    
    tabix -p vcf "$DIR_RESULTS/${STRAIN_1_ID}_missense.vcf.gz"
    tabix -p vcf "$DIR_RESULTS/${STRAIN_2_ID}_missense.vcf.gz"
    
    # 4. Find Intersection and Complements (Unique/Shared)
    # Using -n~11 to find variants in at least one file (creates separate outputs)
    log "Calculating intersections..."
    bcftools isec -p "$DIR_RESULTS/comparison" \
        "$DIR_RESULTS/${STRAIN_1_ID}_missense.vcf.gz" \
        "$DIR_RESULTS/${STRAIN_2_ID}_missense.vcf.gz"
    
    # bcftools isec output explanation:
    # 0000.vcf = unique to first file (STRAIN_1)
    # 0001.vcf = unique to second file (STRAIN_2)
    # 0002.vcf = shared between both files
    # 0003.vcf = present in both (alternative representation)
    
    # Check which files were created
    log "Checking generated comparison files..."
    ls -lh "$DIR_RESULTS/comparison/"
    
    # Safely rename outputs (check if files exist first)
    if [ -f "$DIR_RESULTS/comparison/0000.vcf" ]; then
        mv "$DIR_RESULTS/comparison/0000.vcf" "$DIR_RESULTS/${STRAIN_1_ID}_UNIQUE_missense.vcf"
        log "Created: ${STRAIN_1_ID}_UNIQUE_missense.vcf"
    else
        touch "$DIR_RESULTS/${STRAIN_1_ID}_UNIQUE_missense.vcf"
        log "Warning: No unique variants found for $STRAIN_1_ID"
    fi
    
    if [ -f "$DIR_RESULTS/comparison/0001.vcf" ]; then
        mv "$DIR_RESULTS/comparison/0001.vcf" "$DIR_RESULTS/${STRAIN_2_ID}_UNIQUE_missense.vcf"
        log "Created: ${STRAIN_2_ID}_UNIQUE_missense.vcf"
    else
        touch "$DIR_RESULTS/${STRAIN_2_ID}_UNIQUE_missense.vcf"
        log "Warning: No unique variants found for $STRAIN_2_ID"
    fi
    
    if [ -f "$DIR_RESULTS/comparison/0002.vcf" ]; then
        mv "$DIR_RESULTS/comparison/0002.vcf" "$DIR_RESULTS/SHARED_missense.vcf"
        log "Created: SHARED_missense.vcf"
    elif [ -f "$DIR_RESULTS/comparison/0003.vcf" ]; then
        mv "$DIR_RESULTS/comparison/0003.vcf" "$DIR_RESULTS/SHARED_missense.vcf"
        log "Created: SHARED_missense.vcf (from 0003)"
    else
        touch "$DIR_RESULTS/SHARED_missense.vcf"
        log "Warning: No shared variants found between strains"
    fi
    
    # Clean up any remaining files
    rm -f "$DIR_RESULTS/comparison/"*.vcf
    
    # Generate statistics
    COUNT_1=$(grep -v "^#" "$DIR_RESULTS/${STRAIN_1_ID}_UNIQUE_missense.vcf" 2>/dev/null | wc -l || echo "0")
    COUNT_2=$(grep -v "^#" "$DIR_RESULTS/${STRAIN_2_ID}_UNIQUE_missense.vcf" 2>/dev/null | wc -l || echo "0")
    COUNT_SHARED=$(grep -v "^#" "$DIR_RESULTS/SHARED_missense.vcf" 2>/dev/null | wc -l || echo "0")
    
    # Get total missense count
    COUNT_TOTAL=$(zcat "$DIR_RESULTS/all_missense.vcf.gz" 2>/dev/null | grep -v "^#" | wc -l || echo "0")
    
    # Create summary report
    echo "================================================" > "$DIR_RESULTS/Summary_Report.txt"
    echo "    Comparative Missense Variant Analysis" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "================================================" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Analysis Date: $(date +'%Y-%m-%d %H:%M:%S')" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Strains Compared:" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "  - Strain 1: $STRAIN_1_ID" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "  - Strain 2: $STRAIN_2_ID" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "------------------------------------------------" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Results:" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "------------------------------------------------" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Total Missense Variants Detected: $COUNT_TOTAL" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Unique to $STRAIN_1_ID: $COUNT_1" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Unique to $STRAIN_2_ID: $COUNT_2" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Shared between both strains: $COUNT_SHARED" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "================================================" >> "$DIR_RESULTS/Summary_Report.txt"
    
    # Display results to console
    log "=========================================="
    log "Analysis Complete!"
    log "=========================================="
    log "Total Missense Variants: $COUNT_TOTAL"
    log "Unique to $STRAIN_1_ID: $COUNT_1"
    log "Unique to $STRAIN_2_ID: $COUNT_2"
    log "Shared: $COUNT_SHARED"
    log "=========================================="
    log "Full report saved to: $DIR_RESULTS/Summary_Report.txt"
    
    mark_step "ANALYSIS"
fi

log "Pipeline finished successfully."
