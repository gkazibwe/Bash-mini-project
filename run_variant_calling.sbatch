#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH --output=variant_calling.%j.out
#SBATCH --error=variant_calling.%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module purge

# Run pipeline INSIDE the environment (NO activation needed)
conda run -n variant_analysis_env bash variant_calling1.sh
