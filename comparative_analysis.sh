#!/bin/bash

# ==============================================================================
# SCRIPT: comparative_analysis.sh
# DESCRIPTION: Comparative variant analysis pipeline between two microorganism strains.
#              Includes Download, QC, Alignment, Variant Calling, Annotation, and
#              shared/unique missense SNP analysis.
#
# AUTHORS:
#   1. KAZIBWE GEORGE <gkazibwe@gmail.com>
#   2. WASSWA CHARLES LWANGA <wasswacharleslwanga4@gmail.com>
#   3. DDUMBA FRANCIS SEMAKUBA <frncsddumbasema@gmail.com>
#
# REQUIRES: conda activate variant_analysis_env
# ==============================================================================


# --- 1. CONFIGURATION VARIABLES ---
# (Change these variables to analyze different strains/organisms)

# Sample Accessions (SRA Run IDs)
STRAIN_1_ID="SRR1770413" # E. coli Strain A (Example)
STRAIN_2_ID="SRR1770414" # E. coli Strain B (Example)
SAMPLES=($STRAIN_1_ID $STRAIN_2_ID)

# Reference Genome (E. coli K-12 MG1655)
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
REF_NAME="ecoli_k12_ref"
SNPEFF_DB="Escherichia_coli_str_k_12_substr_mg1655" # Database name for snpEff

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
# STEP 2: DOWNLOAD RAW READS (SRA)
# ==============================================================================
if ! check_step "DOWNLOAD_READS"; then
    log "Downloading Raw Reads from NCBI SRA..."
    for sample in "${SAMPLES[@]}"; do
        log "Downloading $sample..."
        # fastq-dump splits files into forward/reverse (_1.fastq, _2.fastq)
        fastq-dump --split-files --gzip --outdir "$DIR_RAW" "$sample"
    done
    mark_step "DOWNLOAD_READS"
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
        
        R1="$DIR_RAW/${sample}_1.fastq.gz"
        R2="$DIR_RAW/${sample}_2.fastq.gz"
        
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
# STEP 6: VARIANT ANNOTATION (SnpEff)
# ==============================================================================
if ! check_step "ANNOTATION"; then
    log "Downloading SnpEff database for $SNPEFF_DB..."
    # Only download if not already present in snpEff config
    snpEff download -v "$SNPEFF_DB"
    
    log "Annotating Variants..."
    # Output as VCF for readability and downstream filtering
    bcftools view "$DIR_VARIANTS/joint_variants.bcf" | \
    snpEff eff -v "$SNPEFF_DB" - > "$DIR_VARIANTS/annotated_variants.vcf"
    
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
    
    INPUT_VCF="$DIR_VARIANTS/annotated_variants.vcf.gz"
    
    # 1. Extract only Missense variants using SnpSift filter
    # "ANN" is the annotation field added by SnpEff
    cat "$INPUT_VCF" | java -jar $(which SnpSift.jar) filter "ANN[*].EFFECT has 'missense_variant'" \
    > "$DIR_RESULTS/all_missense.vcf"
    
    bgzip -f "$DIR_RESULTS/all_missense.vcf"
    tabix -p vcf "$DIR_RESULTS/all_missense.vcf.gz"
    
    # 2. Split the joint VCF into individual VCFs to use 'bcftools isec'
    bcftools view -s "$STRAIN_1_ID" "$DIR_RESULTS/all_missense.vcf.gz" -Oz -o "$DIR_RESULTS/${STRAIN_1_ID}_missense.vcf.gz"
    bcftools view -s "$STRAIN_2_ID" "$DIR_RESULTS/all_missense.vcf.gz" -Oz -o "$DIR_RESULTS/${STRAIN_2_ID}_missense.vcf.gz"
    
    tabix -p vcf "$DIR_RESULTS/${STRAIN_1_ID}_missense.vcf.gz"
    tabix -p vcf "$DIR_RESULTS/${STRAIN_2_ID}_missense.vcf.gz"
    
    # 3. Find Intersection and Complements (Unique/Shared)
    # -p: prefix for output directory
    # -n=2: take inputs with 2 files
    log "Calculating intersections..."
    bcftools isec -p "$DIR_RESULTS/comparison" -n=2 \
        "$DIR_RESULTS/${STRAIN_1_ID}_missense.vcf.gz" \
        "$DIR_RESULTS/${STRAIN_2_ID}_missense.vcf.gz"
        
    # Rename outputs for clarity (bcftools outputs 0000.vcf, 0001.vcf, etc.)
    mv "$DIR_RESULTS/comparison/0000.vcf" "$DIR_RESULTS/${STRAIN_1_ID}_UNIQUE_missense.vcf"
    mv "$DIR_RESULTS/comparison/0001.vcf" "$DIR_RESULTS/${STRAIN_2_ID}_UNIQUE_missense.vcf"
    mv "$DIR_RESULTS/comparison/0002.vcf" "$DIR_RESULTS/SHARED_missense.vcf"
    rm "$DIR_RESULTS/comparison/0003.vcf" # We don't need the file representing records in neither (should be empty)
    
    # Generate simple stats
    COUNT_1=$(grep -v "^#" "$DIR_RESULTS/${STRAIN_1_ID}_UNIQUE_missense.vcf" | wc -l)
    COUNT_2=$(grep -v "^#" "$DIR_RESULTS/${STRAIN_2_ID}_UNIQUE_missense.vcf" | wc -l)
    COUNT_SHARED=$(grep -v "^#" "$DIR_RESULTS/SHARED_missense.vcf" | wc -l)
    
    echo "------------------------------------------------" > "$DIR_RESULTS/Summary_Report.txt"
    echo "Comparative Analysis Report" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Strain 1: $STRAIN_1_ID" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Strain 2: $STRAIN_2_ID" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "------------------------------------------------" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Unique Missense SNPs in $STRAIN_1_ID: $COUNT_1" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Unique Missense SNPs in $STRAIN_2_ID: $COUNT_2" >> "$DIR_RESULTS/Summary_Report.txt"
    echo "Shared Missense SNPs: $COUNT_SHARED" >> "$DIR_RESULTS/Summary_Report.txt"
    
    log "Analysis Complete. Check $DIR_RESULTS/Summary_Report.txt"
    mark_step "ANALYSIS"
fi

log "Pipeline finished successfully."
