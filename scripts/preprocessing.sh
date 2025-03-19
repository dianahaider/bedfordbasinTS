## Raw data location, cleaning and processing
### Data transfer to compute canada for analysis
Raw fastq reads  are located in the workstation + compute canada RSYNC. Make sure the data is formatted correctly, and fastq files are separated by year.
```bash
#from remote workspace with data to local
rsync -av remote@remote rawdata
#from local to compute can.
rsync -av rawdata/ user@remote:~/projectname/rawdata
```

### Script uploads from github

Log in to compute canada, bbmap-env is already installed on the server

Clone the repository for the analysis of the primers fastq files
```bash
git clone https://github.com/jcmcnch/eASV-pipeline-for-515Y-926R.git 
```

Update mamba according to https://github.com/conda-forge/miniforge
and add bioconda and pytorch channel to successfully install conda envs

#edited script setup-scripts/02 and 03- to change the path to databases (based it in ~/)
#change the path to the databases folder to the correct updated location


```bash
find . -type f | xargs perl -pi -e 's/\/home\/$USER/~/g'
```

#copy all data separated by years, move each raw files into /year/00-raw/