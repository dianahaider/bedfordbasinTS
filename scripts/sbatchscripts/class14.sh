#!/bin/bash
#SBATCH --account=def-rbeiko
#SBATCH --mem=64GB
#SBATCH --time=30:00:00
#SBATCH --mail-user=diana.haider@dal.ca
#SBATCH --mail-type=END,FAIL

module load apptainer

for year in  {2017..2021..1}  ##outer loop
do
	apptainer exec -B $PWD:/lustre06/project/6001026/dhaider/bb_data/split_pipeline -B $year/02-PROKs/all_trims:/outputs -B $year/02-PROKs/all_trims:/inputs qiime2-2023.5.sif \
	qiime feature-classifier classify-sklearn --i-reads /inputs/F280R220/DADA2/representative_sequences.qza \
	--i-classifier /home/$USER/databases/qiime2-classification-db/silva-138.1-ssu-nr99-515FY-926R-classifier.qza \
	--output-dir /outputs/F280R220/class \
	--verbose
done
echo 'Generated all 16S 2014 classifications'

