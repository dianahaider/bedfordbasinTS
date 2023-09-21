#1.0 Inspect read quality
mkdir fastqc_out
fastqc -t 1 raw_data/*.fastq.gz -o fastqc_out

cd fastqc_out
multiqc .
cd ..

#2.0 Download qiime environment and activate
#follow download instructions on https://docs.qiime2.org/2023.5/install/native/
source activate qiime2-2023.2

mkdir reads_qza
    
qiime tools import \
   --type SampleData[PairedEndSequencesWithQuality] \
   --input-path raw_data/ \
   --output-path reads_qza/reads.qza \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt

mkdir reads_qza
    
qiime tools import \
   --type SampleData[PairedEndSequencesWithQuality] \
   --input-path raw_data/ \
   --output-path reads_qza/reads.qza \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt

qiime cutadapt trim-paired \
   --i-demultiplexed-sequences reads_qza/reads.qza \
   --p-cores 1 \
   --p-front-f GTGYCAGCMGCCGCGGTAA \
   --p-front-r CCGYCAATTYMTTTRAGTTT \
   --p-discard-untrimmed \
   --p-no-indels \
   --o-trimmed-sequences reads_qza/reads_trimmed.qza

#Run dada2 and classification on compute canada
#from local: 
scp -r $year $user@computecanada.ca:$PATH

#**move to compute canada**
#from remote:
singularity exec -B $PWD:PATH -B 2017:/outputs \
  -B $year:/inputs qiime2-2023.2.sif \
  qiime dada2 denoise-paired --i-demultiplexed-seqs /inputs/reads_trimmed.qza --p-trunc-len-f 280 --p-trunc-len-r 280 --output-dir /outputs/dada2_output --verbose

wget https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza

singularity --mem=16GB exec -B $PWD:PATH -B 2014:/outputs \
  -B 2014/dada2_output:/inputs -B 2014:/class qiime2-2023.2.sif \
  qiime feature-classifier classify-sklearn --i-reads /inputs/representative_sequences.qza --i-classifier /class/silva-138-99-nb-classifier.qza --output-dir /outputs/taxa --verbose

#from local
qiime feature-table summarize --i-table dada2_output_270210/table.qza --o-visualization dada2_output_270210/dd2270210_table_summary.qzv

qiime feature-table filter-features --i-table dada2_output_270210/table.qza --p-min-frequency 34 --p-min-samples 1 --o-filtered-table dada2_output_270210/table_filt.qza

qiime taxa filter-table --i-table dada2_output_270210/table_filt.qza --i-taxonomy taxa_270210/classification.qza --p-exclude mitochondria,chloroplast --o-filtered-table dada2_output_270210/table_filt_contam.qza

qiime feature-table summarize --i-table dada2_output_270210/table_filt_contam.qza --o-visualization dada2_output_270210/table_filt_contam_summary.qzv


qiime diversity alpha-rarefaction --i-table dada2_output_270210/table_filt_contam.qza --p-max-depth 34 --p-steps 20 --p-metrics 'observed_features' --o-visualization rarefaction_curves_test_270210.qzv

qiime feature-table filter-samples --i-table dada2_output_270210/table_filt_contam.qza --p-min-frequency 7500 --o-filtered-table dada2_output_270210/table_filt_min.qza


qiime feature-table filter-seqs --i-data dada2_output_270210/representative_sequences.qza --i-table dada2_output_270210/table_filt_contam.qza --o-filtered-data dada2_output_270210/rep_seqs_filt_contam_final.qza


qiime fragment-insertion sepp --i-representative-sequences dada2_output_270210/rep_seqs_filt_contam_final.qza --i-reference-database sepp-refs-gg-13-8.qza --o-tree asvs-tree.qza --o-placements insertion-placements.qza


qiime diversity core-metrics-phylogenetic --i-table dada2_output_270210/table_filt_contam.qza --i-phylogeny asvs-tree.qza --p-sampling-depth 8000 --m-metadata-file METADATA_2.tsv --p-n-jobs-or-threads 1 --output-dir diversity --verbose


qiime composition ancom --i-table dada2_output_270210/table_filt_contam_pseudocount.qza --m-metadata-file METADATA_2.tsv --m-metadata-column 'Depth code' --output-dir ancom_output

#export biom table
qiime tools export    --input-path dada2_output/table_filt_contam.qza    --output-path dada2_output_exported

biom convert -i feature-table.biom -o feature-table.tsv --to-tsv






** SPLIT PIPELINE **
#Log in compute can, bbmap-env is already installed on the server
ssh -Y $user@computecanada.ca

#follow guidelines in README file to split 16S/18S
#from local, /ch2
mkdir split_pipeline

#updated mamba according to https://github.com/conda-forge/miniforge
#also add bioconda and pytorch channel to successfully install conda envs

#edited script setup-scripts/02 and 03- to change the path to databases (based it in ~/)
#change the path to the databases folder to the correct updated location
find . -type f | xargs perl -pi -e 's/\/home\/$USER/~/g'

#copy all data separated by years, move each raw files into /year/00-raw/
for file in *.fastq.gz; do mv "$file" "${file/_001.fastq.gz/.fastq.gz}"; done #rename so the script 

#run 01, libgcc-ng, libstdcxx-ng only exist for Linux, for mac
brew install gcc
#01 still fails so install independently:
conda create -n cutadaptenv cutadapt
#remove the '.' from .R1 and .R2 in 515FY926R.cfg


#############
git clone https://github.com#/jcmcnch/eASV-pipeline-for-515Y-926R.git
#if for first time; install all envs and prereqs
#cd eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/
source eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/00-trimming-sorting-scripts/00-run-cutadapt.sh
source eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/00-trimming-sorting-scripts/01-sort-16S-18S-bbsplit.sh 


#on compute canada
#1 transfer all files from my comp to dhaider@narval
#its too slow to use scp to transfer to computecan so use rsync
rsync -avP $year $user@narval.computecanada.ca:PATH
#log in to cc and cd to PATH
#here all years will be there
git clone https://github.com/jcmcnch/eASV-pipeline-for-515Y-926R.git

#create an env and install cutadapt in it
#updated mamba according to https://github.com/conda-forge/miniforge
#also add bioconda and pytorch channel to successfully install conda envs
virtualenv --no-download cutadapt-env
source cutadaptenv/bin/activate
pip install cutadapt

#edited script setup-scripts/02 and 03- to change the path to databases (based it in ~/)
#change the path to the databases folder to the correct updated location
find . -type f | xargs perl -pi -e 's/\/home\/$USER/~/g'

#copy all data separated by years, move each raw files into /year/00-raw/
for file in *.fastq.gz; do mv "$file" "${file/_001.fastq.gz/.fastq.gz}"; done #rename so the script 

#run 01, libgcc-ng, libstdcxx-ng only exist for Linux, for mac
brew install gcc
#01 still fails so install independently:
conda create -n cutadaptenv cutadapt
#remove the '.' from .R1 and .R2 in 515FY926R.cfg


#edit the script to add the variables and to specify
#rawFileEndingR1=R1.fastq.gz
#rawFileEndingR2=R2.fastq.gz
nano PATH/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/00-trimming-sorting-scripts/00-run-cutadapt.sh

sbatch PATH/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/00-trimming-sorting-scripts/01-sort-16S-18S-bbsplit.sh 

#for the yers i ran bbsplit without correcting the tsv outputs:
source PATH/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/02-utility-scripts/calc-EUK-fraction.sh > 230903-1126..EUKfrac-whole-dataset-after-bbpsplit.tsv 

#After bbsplit; cd into PROKS folder
#create manifest 
source ~PATH/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/01-prok-scripts/P00-create-manifest.sh
#import to qiime
#premierement build qiimes image
module load apptainer
apptainer build qiime2-2023.5.sif docker://quay.io/qiime2/core:2023.5

#### PROKARYOTES ####

#then write import.sh 
# apptainer has a hard time with symlinks so $PWD has to be the physical PWD (pwd -P)
*************************
#!/bin/bash
#SBATCH --account=$ACCOUNT
#SBATCH --mem=64GB
#SBATCH --time=04:00:00

module load apptainer

apptainer exec -B $PWD:PATH -B reads_qza:/outputs qiime2-2023.5.sif qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.tsv \
--output-path outputs/16s.qza \
--input-format PairedEndFastqManifestPhred33
****************************

#then run import for each 02-PROKs
sbatch import.sh










