# Bash Mini Project

## Comparative variant analysis between strains

### Objective

To write a bash script that performs all the key steps involved in variant analysis for priority microorganisms (raw sequence data), including quality checks, indexing, alignment, variant calling, variant annotation, and effect prediction. This project aims to provide hands-on experience in using Unix and shell scripting for one of the commonest bioinformatics workflows.  

### Bash Script Requirements

Write a single bash script that automates the comparative variant analysis process between two biological strains of the micro organism of your choice. Ensure that your script can do the following as the bare minimum:  

- Download the appropriate reference genome
- Download paired-end raw sequencing reads (FASTQ files) from NCBI for each of the two strains of your micro organism
- Perform quality checks on the raw sequencing data. Merge the quality checks results into a single html report
- Index the reference microorganism’s genome and align the sample reads to the reference genome.
- Convert SAM to BAM files, merge, sort, and index the sorted BAM files
- Call variants
- Annotate the variants and predict their effects
- Identify unique and shared missense SNPs between strains
- Store all the outputs from each step in properly named output directories
- Implement check points to ensure that the script can resume from the last completed step in case of an interruption.
- Run successfully on the High Performance Computing server without any errors

### Additional Guidelines
- Ensure that you create a conda environment with all the necessary tools required for your script to run successfully
- The script should be scalable and reusable for comparing any two strains of your chosen microorganism
- Use meaningful variable and directory names to improve the script's clarity.
- Test your script thoroughly to ensure it works without any errors before the final submission.
- Strain information file. Provide a text file containing the sample IDs for your two chosen strains
- Environment.yml file. Provide the configuration file with the name of the conda environment, channels, and all the tools required to run your bash script.
