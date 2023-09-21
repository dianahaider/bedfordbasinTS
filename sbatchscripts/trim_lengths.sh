#!/bin/bash
#SBATCH --account=def-rbeiko
#SBATCH --mem=64GB
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=diana.haider@dal.ca
#SBATCH --mail-type=END,FAIL

module load apptainer

for flen in {220..300..10} ##outer loop
do
	for rlen in {220..300..10} ##inner loop
	do
		cd 2014/02-PROKs
		mkdir -p all_trims/"F"$flen"R"$rlen
		cd ..
		echo $PWD
		cd ..
		echo $PWD
		apptainer exec -B $PWD:/lustre06/project/6001026/dhaider/bb_data/split_pipeline -B 2014/02-PROKs/all_trims:/outputs -B 2014/02-PROKs:/inputs qiime2-2023.5.sif \
		qiime dada2 denoise-paired --i-demultiplexed-seqs /inputs/16s.qza \
		--p-trunc-len-f $flen \
		--p-trunc-len-r $rlen \
		--p-n-threads 20 \
		--output-dir /outputs/"F"$flen"R"$rlen/DADA2
		--verbose
	done
done
echo 'Generated all 16S trim combinations'
