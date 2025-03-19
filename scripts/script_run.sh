#clone repo with pipeline from Furhman's lab
git clone https://github.com/jcmcnch/eASV-pipeline-for-515Y-926R.git

#make files executable
chmod a+x ./*

#change qiime version used in code to 2023.2
find . -type f | xargs perl -pi -e 's/qiime2-2023.5/qiime2-2023.2/g'

#setup from README from pipeline
cd eASV-pipeline-for-515Y-926R/
./setup-scripts/00-install-qiime2-2022.2.sh
./setup-scripts/01-install-conda-envs.sh
./setup-scripts/02-download-qiime2-classifiers-qiime2-2022.2.sh
./setup-scripts/03-make-bbsplit-db.sh




#run first few commands
source eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/00-trimming-sorting-scripts/00-run-cutadapt.sh


