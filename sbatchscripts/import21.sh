#!/bin/bash
#SBATCH --account=def-rbeiko
#SBATCH --mem=64GB
#SBATCH --time=04:00:00

module load apptainer

apptainer exec -B $PWD:/lustre06/project/6001026/dhaider/bb_data/split_pipeline/2021/02-PROKs /lustre06/project/6001026/dhaider/bb_data/split_pipeline/qiime2-2023.5.sif qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.tsv \
--output-path 16s.qza \
--input-format PairedEndFastqManifestPhred33V2
