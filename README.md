# average interfacePAE score

## Description:
This Python code was developed to postprocess pdb and json files obtained from ColabFold / AlphaFold2 Multimer. 
It extracts coordinates from pdb files, identifies putative contact sites and calculates an average interface predicted aligned error (aiPAE) scores using PAE matrices. 
This aiPAE represents the average PAE of a given interface. It operates on a directory (./pdbs) containing folders with PDB and the respective JSON files, identifies relevant pairs, 
and then calculates (likely) interactions between chains in the protein structures. It writes results to an output file ans logs important messages throughout the process.<br>
## Table of Contents
- [Requirements](#requirements)
- [Running the code](#running-the-code)
- [List of Functions](#list-of-functions)
- [Output](#output)
- [Acknowledgements](#acknowledgements)
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)

## Requirements (tested versions):
- `Python`: 3.11.7  
- `numpy`: 1.26.4
- `pandas`: 2.1.4
  
:exclamation: **Identifiers:**
PDB and JSON file names must follow a certain pattern, as they are used as key for matching the correct PDB and JSON pairs. 
To obtain these file names, the following naming must be entered in the “jobname” field of Colabfold, when entering protein dimers: Number_ID1_ID2. 
In our project, we have assigned a consecutive number followed by the Uniprot ID of sequence A (ID1) and the Uniprot ID of sequence B (ID2). Example: 0001_P62258_P29310
