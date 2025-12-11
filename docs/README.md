# SignalScreen
A python implementation to find candidate compounds

## Overview   
- Mine literature/databases to find tool-compounds for  RAR, especially RARG.
- Need for retinoid-like activity with lower irritation; enabling primary screening; providing medicinal-chemistry-ready data.  
- Contains python code, data, and documentation

## Repository Structure  
```
project_root/
│
├── data/                           # raw & processed data
│   ├── targets/                    # initial targets of interest
│   ├── raw/                        # original downloads / dumps from external sources (ChEMBL export, PubChem JSON, BindingDB entries, etc.)
│   ├── processed/                  # cleaned / normalized data (e.g. merged bioactivity tables, deduplicated compounds, standardized IDs)
│   └── output/                     # final compound list(s) plus any summary tables or filtered “top candidates”
│
├── src/                            # source code / scripts
│   ├── fetch_gtopdb_targets.py     # script to fetch GtoP DB target
│   ├── fetch_chembl_for_targets.py # script to fetch general compound metadata (SMILES, identifiers, synonyms)  
│   ├── fetch_bioactivity.py        # script to fetch bioactivity data for specified targets (binding, functional assays)  
│   ├── normalize.py                # code to clean, normalize, dedupe metadata + bioactivity, unify identifiers/synonyms  
│   ├── filter_prioritize.py        # logic to apply potency/selectivity filters and choose “top candidates”  
│   └── export.py                   # script to export final lists / outputs (CSV/JSON/SQLite)  
│
├── docs/                         # documentation & notes  
│   ├── README.md                 # this file  
│   └── data_flow_diagram.drawio  # visualization of data flow from raw fetch → processing → output  
│
├── requirements.txt              # list of Python (or other) dependencies needed to run the scripts  
├── run_pipeline.sh               # convenience shell script to run the entire pipeline end-to-end  
└── .gitignore                    # standard ignore file (e.g. for large raw downloads, API keys, temporary files)  
```


