#!/bin/bash
# ==============================================================================
# SCRIPT: comparative_analysis.sh
# DESCRIPTION: Comparative variant analysis pipeline between two HSV-1 strains
#	       Includes Download, QC, Alignment, Variant Calling, Annotation, and
#              shared/unique missense SNP analysis.
#              Hybrid version: runs both on HPC and local computers
# AUTHORS:
#   1. KAZIBWE GEORGE 2025/HD07/25965U <gkazibwe@gmail.com>
#   2. WASSWA CHARLES LWANGA 2025/HD07/26027U <wasswacharleslwanga4@gmail.com>
#   3. DDUMBA FRANCIS SEMAKUBA 2025/HD07/25984U <frncsddumbasema@gmail.com>
# ==============================================================================

# --- 0. ENVIRONMENT DETECTION ---
HPC_ENV=false
if [ -n "$SLURM_JOB_ID" ] || [ -n "$PBS_JOBID" ]; then
    HPC_ENV=true
fi

# --- 1. CONFIGURATION VARIABLES ---
# Sample Accessions (SRA Run IDs)
STRAIN_1_ID="SRR36143512" # Herpes simplex first Strain
STRAIN_2_ID="SRR23265797" # Herpes simplex second Strain
SAMPLES=($STRAIN_1_ID $STRAIN_2_ID)

# Reference Genome (Human herpesvirus 1 strain 17)
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/859/985/GCF_000859985.2_ViralProj15217/GCF_000859985.2_ViralProj15217_genomic.fna.gz"
REF_NAME="Human_herpesvirus_1_strain_17"
SNPEFF_DB="HSV1_strain17" # Database name for snpEff

THREADS=8
JAVA_MEM="16g"  # Java memory for snpEff/SnpSift

# --- 2. DIRECTORIES ---
BASE_DIR=$(pwd)
if [ "$HPC_ENV" = true ]; then
    echo "HPC environment detected. Using scratch storage..."
    BASE_DIR="${SLURM_TMPDIR:-/tmp}/hsv1_pipeline"
fi
mkdir -p "$BASE_DIR"
DIR_REF="$BASE_DIR/1_reference"
DIR_RAW="$BASE_DIR/2_raw_reads"
DIR_QC="$BASE_DIR/3_qc_reports"
DIR_ALIGN="$BASE_DIR/4_alignment"
DIR_VARIANTS="$BASE_DIR/5_variants"
DIR_RESULTS="$BASE_DIR/6_final_results"
CHECKPOINT_FILE="$BASE_DIR/.checkpoint"

# Create Directories
mkdir -p $DIR_REF $DIR_RAW $DIR_QC $DIR_ALIGN $DIR_VARIANTS $DIR_RESULTS
touch "$CHECKPOINT_FILE"

# --- 3. ERROR HANDLING & LOGGING ---

# Exit immediately if a command exits with a non-zero status
set -e
set -o pipefail

log() { echo -e "\n[$(date +'%Y-%m-%d %H:%M:%S')] $1"; }
check_step() { grep -q "$1" "$CHECKPOINT_FILE" 2>/dev/null && return 0 || return 1; }
mark_step() { echo "$1" >> "$CHECKPOINT_FILE"; log "Step '$1' completed."; }

log "Starting Comparative Analysis Pipeline..."
log "Comparing $STRAIN_1_ID vs $STRAIN_2_ID"

# ==============================================================================
# STEP 1: DOWNLOAD REFERENCE GENOME
# ==============================================================================
if ! check_step "DOWNLOAD_REF"; then
    log "Downloading Reference Genome..."
    wget -O "$DIR_REF/${REF_NAME}.fna.gz" "$REF_URL"
    gunzip -f "$DIR_REF/${REF_NAME}.fna"
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
# STEP 2: DOWNLOAD RAW READS
# ==============================================================================
if ! check_step "DOWNLOAD_READS"; then
    log "Downloading Raw Reads from NCBI SRA..."
    for sample in "${SAMPLES[@]}"; do
        log "Downloading $sample..."
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
# STEP 3: QUALITY CONTROL
# ==============================================================================
if ! check_step "RUN_QC"; then
    log "Running FastQC..."
    fastqc -t $THREADS "$DIR_RAW"/*.fastq.gz -o "$DIR_QC"
    log "Aggregating reports with MultiQC..."
    multiqc "$DIR_QC" -o "$DIR_QC" -n "Combined_QC_Report"
    mark_step "RUN_QC"
fi

# ==============================================================================
# STEP 4: ALIGNMENT
# ==============================================================================
if ! check_step "ALIGNMENT"; then
    log "Indexing Reference Genome..."
    bwa index "$DIR_REF/${REF_NAME}.fna"

    for sample in "${SAMPLES[@]}"; do
        log "Aligning $sample..."
        R1="$DIR_RAW/${sample}_1.repaired.fastq.gz"
        R2="$DIR_RAW/${sample}_2.repaired.fastq.gz"
        RG_TAG="@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA"

        bwa mem -t $THREADS -R "$RG_TAG" "$DIR_REF/${REF_NAME}.fna" "$R1" "$R2" | \
        samtools view -bS - | \
        samtools sort -@ $THREADS -o "$DIR_ALIGN/${sample}.sorted.bam"

        samtools index "$DIR_ALIGN/${sample}.sorted.bam"
    done
    mark_step "ALIGNMENT"
fi

# ==============================================================================
# STEP 5: VARIANT CALLING
# ==============================================================================
if ! check_step "VARIANT_CALLING"; then
    log "Calling variants..."
    bcftools mpileup -f "$DIR_REF/${REF_NAME}.fna" \
        "$DIR_ALIGN/$STRAIN_1_ID.sorted.bam" \
        "$DIR_ALIGN/$STRAIN_2_ID.sorted.bam" | \
    bcftools call -mv -Ob -o "$DIR_VARIANTS/joint_variants.bcf"
    mark_step "VARIANT_CALLING"
fi

# ==============================================================================
# STEP 5b: BUILD SNPEFF DATABASE
# ==============================================================================
if ! check_step "BUILD_SNPEFF_DB"; then
    log "Building SnpEff database..."
    DEFAULT_CONFIG=$(find "$CONDA_PREFIX" -name snpEff.config | head -n 1)
    if [ -z "$DEFAULT_CONFIG" ]; then
        log "Error: snpEff.config not found"
        exit 1
    fi

    cp "$DEFAULT_CONFIG" ./snpEff.config
    sed -i 's|^data.dir = .*|data.dir = ./data|' ./snpEff.config
    if ! grep -q "$SNPEFF_DB.genome" ./snpEff.config; then
        echo -e "\n# Custom Genome: HSV-1" >> ./snpEff.config
        echo -e "$SNPEFF_DB.genome : Human herpes simplex virus 1 (strain 17)" >> ./snpEff.config
    fi

    mkdir -p "./data/$SNPEFF_DB"
    cp "$DIR_REF/${REF_NAME}.fna" "./data/$SNPEFF_DB/sequences.fa"
    cp "$DIR_REF/${REF_NAME}.gff" "./data/$SNPEFF_DB/genes.gff"

    snpEff build -Xmx$JAVA_MEM -c ./snpEff.config -gff3 -noCheckCds -noCheckProtein -v "$SNPEFF_DB"
    mark_step "BUILD_SNPEFF_DB"
fi

# ==============================================================================
# STEP 6: VARIANT ANNOTATION
# ==============================================================================
if ! check_step "ANNOTATION"; then
    log "Annotating variants..."
    bcftools view "$DIR_VARIANTS/joint_variants.bcf" | \
    snpEff -Xmx$JAVA_MEM eff -c ./snpEff.config -v "$SNPEFF_DB" - > "$DIR_VARIANTS/annotated_variants.vcf"

    bgzip -f "$DIR_VARIANTS/annotated_variants.vcf"
    tabix -p vcf "$DIR_VARIANTS/annotated_variants.vcf.gz"
    mark_step "ANNOTATION"
fi

# ==============================================================================
# STEP 7: MISSENSE VARIANTS ANALYSIS
# ==============================================================================
if ! check_step "ANALYSIS"; then
    log "Filtering missense variants..."
    SNPSIFT_JAR=$(find "$CONDA_PREFIX/share" -iname "SnpSift.jar" | head -n1)
    INPUT_VCF="$DIR_VARIANTS/annotated_variants.vcf.gz"

    zcat "$INPUT_VCF" | java -Xmx$JAVA_MEM -jar "$SNPSIFT_JAR" filter \
        "ANN[*].EFFECT has 'missense_variant'" > "$DIR_RESULTS/all_missense.vcf"
    bgzip -f "$DIR_RESULTS/all_missense.vcf"
    tabix -p vcf "$DIR_RESULTS/all_missense.vcf.gz"

    # Split by sample
    for sample in "${SAMPLES[@]}"; do
        bcftools view -s "$sample" "$DIR_RESULTS/all_missense.vcf.gz" -Oz -o "$DIR_RESULTS/${sample}_missense.vcf.gz"
        tabix -p vcf "$DIR_RESULTS/${sample}_missense.vcf.gz"
    done

    # Intersections and unique variants
    mkdir -p "$DIR_RESULTS/comparison"
    bcftools isec -p "$DIR_RESULTS/comparison" \
        "$DIR_RESULTS/${STRAIN_1_ID}_missense.vcf.gz" \
        "$DIR_RESULTS/${STRAIN_2_ID}_missense.vcf.gz"

    # Rename outputs
    mv "$DIR_RESULTS/comparison/0000.vcf" "$DIR_RESULTS/${STRAIN_1_ID}_UNIQUE_missense.vcf" 2>/dev/null || touch "$DIR_RESULTS/${STRAIN_1_ID}_UNIQUE_missense.vcf"
    mv "$DIR_RESULTS/comparison/0001.vcf" "$DIR_RESULTS/${STRAIN_2_ID}_UNIQUE_missense.vcf" 2>/dev/null || touch "$DIR_RESULTS/${STRAIN_2_ID}_UNIQUE_missense.vcf"
    if [ -f "$DIR_RESULTS/comparison/0002.vcf" ]; then
        mv "$DIR_RESULTS/comparison/0002.vcf" "$DIR_RESULTS/SHARED_missense.vcf"
    elif [ -f "$DIR_RESULTS/comparison/0003.vcf" ]; then
        mv "$DIR_RESULTS/comparison/0003.vcf" "$DIR_RESULTS/SHARED_missense.vcf"
    else
        touch "$DIR_RESULTS/SHARED_missense.vcf"
    fi
    rm -f "$DIR_RESULTS/comparison/"*.vcf

    # Generate summary
    COUNT_1=$(grep -v "^#" "$DIR_RESULTS/${STRAIN_1_ID}_UNIQUE_missense.vcf" 2>/dev/null | wc -l || echo "0")
    COUNT_2=$(grep -v "^#" "$DIR_RESULTS/${STRAIN_2_ID}_UNIQUE_missense.vcf" 2>/dev/null | wc -l || echo "0")
    COUNT_SHARED=$(grep -v "^#" "$DIR_RESULTS/SHARED_missense.vcf" 2>/dev/null | wc -l || echo "0")
    COUNT_TOTAL=$(zcat "$DIR_RESULTS/all_missense.vcf.gz" 2>/dev/null | grep -v "^#" | wc -l || echo "0")

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

    log "=========================================="
    log "Analysis Complete!"
    log "Total Missense Variants: $COUNT_TOTAL"
    log "Unique to $STRAIN_1_ID: $COUNT_1"
    log "Unique to $STRAIN_2_ID: $COUNT_2"
    log "Shared: $COUNT_SHARED"
    log "Full report saved to: $DIR_RESULTS/Summary_Report.txt"
    log "=========================================="

    mark_step "ANALYSIS"
fi

# ==============================================================================
# STEP 8: COPY BACK RESULTS IF HPC
# ==============================================================================
if [ "$HPC_ENV" = true ]; then
    log "Copying results back to home directory..."
    cp -r "$BASE_DIR" "$PWD/"
fi

log "Pipeline finished successfully."
