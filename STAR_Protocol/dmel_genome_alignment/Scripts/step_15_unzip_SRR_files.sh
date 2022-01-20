#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-egreenbl
#SBATCH --mem=32G
#SBATCH --job-name=Unzip_SSR
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

gunzip *
