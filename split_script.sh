** SPLIT PIPELINE **
#Log in compute can, bbmap-env is already installed on the server
ssh -Y dhaider@narval.computecanada.ca


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
rsync -avP 2018 dhaider@narval.computecanada.ca:~/projects/def-rbeiko/dhaider/bb_data/split_pipeline
#log in to cc and cd to ~/projects/def-rbeiko/dhaider/bb_data/split_pipeline
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
nano ~/projects/def-rbeiko/dhaider/bb_data/split_pipeline/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/00-trimming-sorting-scripts/00-run-cutadapt.sh


sbatch ~/projects/def-rbeiko/dhaider/bb_data/split_pipeline/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/00-trimming-sorting-scripts/01-sort-16S-18S-bbsplit.sh 

#for the yers i ran bbsplit without correcting the tsv outputs:
source ~/projects/def-rbeiko/dhaider/bb_data/split_pipeline/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/02-utility-scripts/calc-EUK-fraction.sh > 230903-1126..EUKfrac-whole-dataset-after-bbpsplit.tsv 

#After bbsplit; cd into PROKS
#create manifest 
source ~/projects/def-rbeiko/dhaider/bb_data/split_pipeline/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/01-prok-scripts/P00-create-manifest.sh
#import to qiime
#premierement build qiimes image
module load apptainer
apptainer build qiime2-2023.5.sif docker://quay.io/qiime2/core:2023.5

#### PROKARYOTES ####

#then write import.sh 
# apptainer has a hard time with symlinks so $PWD has to be the physical PWD (pwd -P)
*************************
#!/bin/bash
#SBATCH --account=def-rbeiko
#SBATCH --mem=64GB
#SBATCH --time=04:00:00

module load apptainer

apptainer exec -B $PWD:/lustre06/project/6001026/dhaider/bb_data/split_pipeline/2014/02-PROKs -B reads_qza:/outputs qiime2-2023.5.sif qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.tsv \
--output-path outputs/16s.qza \
--input-format PairedEndFastqManifestPhred33
****************************

#then run import for each 02-PROKs
sbatch import.sh


#from julies to local
rsync -av larochelabdhwani1@129.173.32.88:~/Diana/BB2014-2021a/raw_data/ rawdata
#from local to compute can.
rsync -av rawdata/ dhaider@narval.computecanada.ca:~/projects/def-rbeiko/dhaider/bb_data/split_pipeline/rawdata
