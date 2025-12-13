# SignalScreen
A python implementation to find candidate compounds

## Overview   
- Mine databases to find tool-compounds for  RAR, especially RARG.
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
│   ├── sort_rarg_ligands.py        # Rank ligands by activity and specificity toward retinoic acid receptor gamma (RARG)
│   ├── fetch_compound_info.py      # Enrich the top N sorted ligands with PubChem CID and synonyms
│   ├── fetch_target_info.py        # Find compound target info and call Reactome's AnalysisService 
│   ├── fetch_clinical_trials.py    # Queries CT.gov for name/synonym matches for each ligand
│   ├── fetch_literature.py         # Search for literature about ligands on pubmed, chembl, and bindingdb
│   ├── build_ligand_summary.py     # Build and save a summary DB
│   ├── open_targets_known_drugs.py # Query Open Targets for known drugs of a target
│   └── fetch_top_drugs.py          # Given df_ot_drugs, return each unique drug with its highest phase
│
├── docs/                           # documentation & notes  
│   ├── README.md                   # this file  
│   └── data_flow_diagram.drawio    # visualization of data flow from raw fetch → processing → output  
│
├── requirements.txt                # list of Python (or other) dependencies needed to run the scripts  
├── RAR_candidates.ipynb            # python notebook to run full analysis example 
└── .gitignore                      # standard ignore file (e.g. for large raw downloads, API keys, temporary files)  
```


