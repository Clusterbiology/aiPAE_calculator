# average interface PAE calculator

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
- [Funding] (#funding)
- [Contact](#contact)

## Requirements:
- `Python`: 3.11.7  
- `numpy`: 1.26.4
- `pandas`: 2.1.4
  
:exclamation: **Identifiers:**
PDB and JSON file names must follow a certain pattern, as they are used as key for matching the correct PDB and JSON pairs. 
To obtain these file names, the following naming must be entered in the “jobname” field of Colabfold, when entering protein dimers: Number_ID1_ID2. 
In our project, we have assigned a consecutive number followed by the Uniprot ID of sequence A (ID1) and the Uniprot ID of sequence B (ID2). Example: 0001_P62258_P29310


```python

""" The file names in this project follow a specific format, which is composed of several parts separated by underscores (_).
For correct matching of pdb and json files it is important to adapt the file names or to rewrite the code for generating unique file keys."""
        # 7) Generate a unique filename key
        base_name = os.path.basename(pdb_file_path).split('_')
        key = "_".join([base_name[i] for i in [0, 1, 2, 3, 5, 6] if i < len(base_name)])
        # These parts of the file name need to match between the pdb and the json file.
        
        # In this project the relevant parts of the file name are:
            # part0: Number
            # part1: Chain A ID (UniProt ID of Protein A)
            # part2: Chain B ID (UniProt ID of Protein B)
            # part3: ID added by alphafold
            # part5: 'rank'
            # part6: rank number
        # e.g. 0001_P63261_Q9NRH3_ab1cd_part4_rank_001_additional_info.pdb
```

## Running the code:

:exclamation: **Important:** After downloading AlphaFold results unzip them into a folder **./pdbs**.

| Script | Description | Input | Output |
|------|------|------|------|
| __aiPAE_calculator.py__ | Processes protein structure data, including extracting coordinates from PDB-files and calculating iPAE scores from PAE matrices. Then loops through the directory, identifies PDB and JSON file-pairs and calculates interactions between chains in the protein structures. | `./pdbs/*/*.pdb`<br>`./pdbs/*/*.json` | `./pdbs/summary_results.csv`<br>`./pdbs/*/*_interaction_results.txt`<br>`./pdbs/*/*_residue_pair_details.txt` |

### List of functions
| Function | Description |
|------|------|
| __ask_for_parallel_processing()__ | Prompts the user to choose whether to use parallel processing. It returns a boolean value:<br>- True if the user inputs 'yes'<br>- False if the user inputs 'no'. |
| __find_continuous_ranges_and_sequences(indices, residues)__ | Identifies continuous ranges of residues (as indices) and generates sequences of the corresponding one-letter amino acid code for these ranges. |
| __load_chain_coords_custom(pdb_file, chain_id)__ | Loads the atomic coordinates for the specified chain (chain_id) from the PDB file. It extracts only the alpha-carbaon (CA) atoms from the given chain. |
| __load_pae(file_path)__ | Loads a PAE matrix from a json file. |
| __find_multiple_interfaces(chain_a_coords, chain_b_coords, distance_threshold)__ | Identifies potential interaction interfaces between two chains based on atomic coordinates. It calculates the distance between each pair of residues from chains A and B, and groups them into interfaces if their distance is below a given threshold. Default distance_threshold=15.0 Å. |
| __calculate_weighted_ipae(interfaces, pae_matrix, chain_a_length, chain_b_length)__ | Calculates weighted iPAE scores for each interface. The iPAE score combines both the distance between residues and the PAE value for the residue pair. |
| __process_pdb_json_pair(pdb_file_path, json_file_path, output_file_path)__ | Processes a single pair of PDB and JSON files, calculates interaction interfaces, computes iPAE scores, and writes the results to output files. |
| __process_directory(directory_path, use_parallel)__ | Scans the pdbs directory for PDB and JSON files, pairs them up, and processes them by extracting interfaces and calculating iPAE scores. The results are summarized and written to a CSV file. |
| __get_relevant_parts(filename)__ | Extracts and returns relevant parts of the file names. File names are split by underscores and relevant parts selected based on their position in the split list. |


## Output: 
Returns a csv-file (summary_results.csv) containing IDs and iPAE Interfaces. IDs are resulting from pdb-file names and iPAE interfaces are calculated from residue distances and PAE values.<br>
Furthermore for each protein dimer two txt-files are created:<br>
- The *_interaction_results.txt contains all identified interfaces with their coordinates and the iPAE.<br>
- The *_residue_pair_details.txt summarises the per residue distance and the PAE of each Pair of residues in an interface.<br>

## Acknowledgements:

This project uses the following third-party libraries:

- [NumPy 1.26.4](https://numpy.org/) for numerical computing.
- [Pandas 2.1.4](https://pandas.pydata.org/) for data manipulation.

## License:

## Funding:

## Contact:
### Contributors
[Stefan Düsterhöft](https://github.com/ClusterBioSD): Conceptualization of the project,main development of the codebase, implementation of core functionalities, testing, and code review.<br>
[Sarah Knapp](https://github.com/knappsarah): Code documentation, addition of code components and code review.
