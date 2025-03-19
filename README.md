# Bedford Basin disturbance analysis

We determined the resistance and resilience of the microbial community in the Bedford Basin in response to rainfall events of varying frequency, and intensity.

### Structure of the starting repository

```Working directory  
root/
├── data/
│   ├── raw/
│   │  └── raw_reads.fastq.gz  
│   └── processed/
│      ├── METADATA.tsv
│      ├── MANIFEST.tsv
│      └── clean_and_processed_reads.qza
│
├── notebooks/
│   ├── 01_metadata_validation.ipynb
│   ├── 02_prelim_data_analysis.ipynb    #alpha, beta, composition analyses
│   └── ...
│
├── scripts/
│   ├── preprocessing.sh
│   └── utils.py                         #python functions for notebooks
│
├── outputs/                             #outputs generated
│
├── README.md
└── environment.yml




 ```
    
### Datasets
| Dataset       | Community     | Link | 
| ------------- | ------------- |------|
| Sequencing    | Even          |
| METADATA      | Staggered     |

### Generate the data
Run ```preprocessing.sh``` and ```generate_data.sh``` from the cloned repo, then open the scripts in order and run all.

### Final repository
