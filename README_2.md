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

### How to run the script on a local computer 
1. Install Conda: Ensure you have `Miniconda` or `Anaconda` installed.
   # Move to the directory where the comparative_analysis.sh and environment.yml were downloaded to 
2. Create Environment:
   # Run this command from the directory where the environment.yml file has been put 
```Bash
conda env create -f environment.yml
```
3. Make Executable (optional):
```Bash
chmod +x comparative_analysis.sh
```
4. Activate the environment
   ```Bash
   conda activate variant_analysis_env
   ```
5. Run:
```Bash
./comparative_analysis.sh
```
Alternatively:  
```Bash
bash comparative_analysis.sh
```

# How to run the script on an HPC 
1. Log into the HPC
```Bash
ssh your_username@hpc_address
```
### Move to your project directory where the comparative_analysis.sh and environment.yml where placed.

2. Create Environment:
   ### load the anaconda module if the hpc uses modules 
   ```Bash
   module load anaconda
   ```
   ### Then acitivate conda
   ```Bash
   source $HOME/miniconda3/etc/profile.d/conda.sh
   ```
   ### Then run this command from the directory where the environment.yml file has been put 
```Bash
conda env create -f environment.yml
```
3. Make Executable:
```Bash
chmod +x comparative_analysis.sh
```
4. Activate the environment
   ```Bash
   conda activate variant_analysis_env
   ```
5. Run:
```Bash
sbatch comparative_analysis.sh
```
Check for the job status:  
```Bash
squeue -u $USER
```
View output:
```Bash
cat slurm-12345.out
```
