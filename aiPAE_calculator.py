import os
import json
import numpy as np
import logging
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Dictionary to convert three-letter residue codes to single-letter codes
three_to_one = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "ASX": "B",  # Ambiguous: Aspartate/Asparagine
    "GLX": "Z",  # Ambiguous: Glutamate/Glutamine
    "XAA": "X",  # Unknown amino acid
    "UNK": "X"   # Unknown residue
}

def find_continuous_ranges_and_sequences(indices, residues):
    # Initialize two lists: one for storing ranges (start, end) and one for storing amino acid sequences
    ranges = []
    sequences = []

    # If the indices list is empty, return empty lists for ranges and sequences
    if not indices:
        return ranges, sequences

    # Initialize 'start' and 'end' to the first index in the indices list.
    # 'start' represents the beginning of a new range, 'end' represents the current position in the range.
    start = end = indices[0]
    # Get the one-letter amino acid code for the residue at the 'start' position.
    # If the residue is not found in 'residues', default to 'X'.
    ## Use the dictionary `.get(key, fallback)` in case `start` isn't in `residues`
    seq = [three_to_one.get(residues.get(start, 'X'), 'X')]

    # Iterate through the remaining indices starting from the second index.
    for idx in indices[1:]:
        # Check if the current index is consecutive to the 'end' of the current range.
        if idx == end + 1:
            # If it is consecutive, extend the current range by updating 'end' and adding the new amino acid to 'seq'.
            end = idx
            seq.append(three_to_one.get(residues.get(idx, 'X'), 'X'))
        else:
            # If the current index is not consecutive, finalise the current range and sequence.
            ranges.append((start, end)) # Store the completed range (start, end)
            sequences.append(''.join(seq)) # Store the corresponding amino acid sequence
            # Start a new range starting at the current index 'idx', and initialize the new sequence.
            start = end = idx
            seq = [three_to_one.get(residues.get(start, 'X'), 'X')]

    # After finishing the loop, append the last range and sequence to the lists.
    ranges.append((start, end)) # Final range
    sequences.append(''.join(seq)) # Final sequence

    # Return the lists of ranges and sequences
    return ranges, sequences

# This function loads the coordinates of the CA atoms for a specific chain from a PDB file.
def load_chain_coords_custom(pdb_file, chain_id):
     # Initialize an empty list to store the coordinates for the given chain.
    chain_coords = []
    try:
        # Open the PDB file for reading.
        with open(pdb_file, 'r') as file:
            # Loop over each line in the PDB file.
            for line in file:
                # Check if the line represents an atom record (ATOM) and is of a valid length.
                # Also, ensure that the chain ID matches the given chain_id and that the atom is a C-alpha (CA) atom.
                if line.startswith("ATOM") and len(line) >= 54 and line[21] == chain_id and line[12:16].strip() == "CA":
                    # Extract the residue index (position in the sequence) from columns 22-26 and convert it to an integer.
                    residue_index = int(line[22:26].strip())
                    # Extract the residue name (three-letter code) from columns 17-20.
                    residue_name = line[17:20].strip()
                    # Extract the x, y, and z coordinates of the C-alpha atom from columns 30-38, 38-46, and 46-54, respectively.
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    # Append the residue index, residue name, and coordinates to the list `chain_coords` in the form of a dictionary.
                    # The coordinates are stored as a numpy array.
                    chain_coords.append({
                        'index': residue_index,
                        'residue': residue_name,
                        'coord': np.array([x, y, z])
                    })
    # Handle case if the specified PDB file is not found.
    except FileNotFoundError:
        logging.error(f"File {pdb_file} not found.")
    # Handle case if there is a parsing error (e.g., invalid data or incorrect formatting in the PDB file).
    except ValueError as e:
        logging.error(f"Issue with parsing PDB file {pdb_file}: {e}")
    # Return the list of chain coordinates.
    return chain_coords

# This function loads a PAE (predicted alignment error) matrix from a JSON file.
def load_pae(file_path):
    try:
        # Open the file specified by file_path in read mode.
        with open(file_path, 'r') as file:
            # Parse the JSON data from the file.
            data = json.load(file)

        # Check if the 'pae' key exists in the JSON data and if its value is a list.
        if 'pae' in data and isinstance(data['pae'], list):
            # Convert the list of PAE values into a numpy array for easier handling and manipulation.
            pae_matrix = np.array(data['pae'])
            # Ensure that the PAE matrix contains numeric values (integers or floats).
            if pae_matrix.dtype.kind in {'i', 'f'}:  # Must be integers or floats
                # If everything is correct, return the PAE matrix as a numpy array.
                return pae_matrix
            else:
                # If the PAE matrix contains non-numeric values, log an error and return None.
                logging.error(f"PAE matrix contains non-numeric values in {file_path}.")
                return None
        else:
            # If the 'pae' key is missing or the format is not correct, log an error and return None.
            logging.error(f"PAE data in {file_path} is not properly formatted or missing.")
            return None
    except (FileNotFoundError, json.JSONDecodeError) as e:
        # Handle errors in opening the file or parsing the JSON (e.g., file not found or invalid JSON format).
        logging.error(f"Error loading PAE file {file_path}: {e}")
        return None

def find_multiple_interfaces(chain_a_coords, chain_b_coords, distance_threshold=15.0):
    # A list to store the identified interfaces
    interfaces = []
    # A set to keep track of already visited residue pairs (to avoid duplicates)
    visited_pairs = set()

    # Log the information that the interface finding process has started with the given distance threshold
    logging.info(f"Finding multiple interfaces between chains with distance threshold {distance_threshold}Å.")

    # Iterate over all residues in chain_a
    for res_a in chain_a_coords:
        # Iterate over all residues in chain_b
        for res_b in chain_b_coords:
            # Calculate the Euclidean distance between the two residues
            distance = np.linalg.norm(res_a['coord'] - res_b['coord'])
            # If the distance is smaller than the threshold and the residue pair has not been visited yet
            if distance < distance_threshold and (res_a['index'], res_b['index']) not in visited_pairs:
                # Mark the pair as visited
                visited_pairs.add((res_a['index'], res_b['index']))
                # A flag to check if the pair was added to an existing interface
                added = False
                # Iterate over all previously found interfaces
                for interface in interfaces:
                    # If the current residue from chain_a is close to a residue in an existing interface
                    if np.linalg.norm(interface[-1]['res_a']['coord'] - res_a['coord']) < distance_threshold:
                        # Add the pair to this interface
                        interface.append({'res_a': res_a, 'res_b': res_b})
                        added = True
                        break
                # If the pair was not added to any existing interface, create a new interface
                if not added:
                    interfaces.append([{'res_a': res_a, 'res_b': res_b}])

    # Log the total number of identified interfaces
    logging.info(f"Identified {len(interfaces)} potential interfaces.")
    # Return the identified interfaces
    return interfaces

def calculate_weighted_ipae(interfaces, pae_matrix, chain_a_length, chain_b_length):
    distance_weights = {'veryhigh': 10.0, 'high': 3.0, 'medium': 1.0, 'low': 0.3, 'verylow': 0.1}
    ipaes = []

    logging.info("Calculating weighted iPAE scores for interfaces.")

    for interface_idx, interface in enumerate(interfaces):
        total_weighted_pae = 0
        total_weight = 0

        for pair_idx, pair in enumerate(interface):
            try:
                distance = np.linalg.norm(pair['res_a']['coord'] - pair['res_b']['coord'])
                distance_weight = (
                    distance_weights['veryhigh'] if distance <= 3 else
                    distance_weights['high'] if distance <= 5 else
                    distance_weights['medium'] if distance <= 8 else
                    distance_weights['low'] if distance <= 10 else
                    distance_weights['verylow'] if distance <= 15 else
                    0
                )

                if distance_weight == 0:
                    logging.debug(f"Skipping pair {pair_idx} in interface {interface_idx} due to distance > 15 Å.")
                    continue

                index_a = pair['res_a']['index'] - 1
                index_b = (pair['res_b']['index'] - 1) + chain_a_length

                # Ensure indices are within bounds of the PAE matrix
                if index_a < 0 or index_a >= pae_matrix.shape[0] or index_b < 0 or index_b >= pae_matrix.shape[1]:
                    logging.warning(f"Index out of bounds for PAE matrix: ({index_a}, {index_b}). Skipping this pair.")
                    continue

                pae_value = min(pae_matrix[index_a, index_b], pae_matrix[index_b, index_a])

                pae_weight = (
                    10 if 0 <= pae_value < 5 else
                    3 if 5 <= pae_value < 10 else
                    1 if 10 <= pae_value < 15 else
                    0.3
                )

                if pae_weight == 0.3 and pae_value > 15:
                    logging.debug(f"Assigning a low weight for pair {pair_idx} in interface {interface_idx} due to high PAE value (>15 Å).")

                combined_weight = distance_weight * pae_weight
                total_weighted_pae += pae_value * combined_weight
                total_weight += combined_weight
            except Exception as e:
                logging.error(f"Error calculating iPAE for pair {pair_idx} in interface {interface_idx}: {e}")
                continue

        if total_weight > 0:
            ipae = total_weighted_pae / total_weight
        else:
            ipae = float('inf')  # No valid pairs

        ipaes.append(ipae)
        logging.info(f"Interface {interface_idx} has iPAE score: {ipae:.2f}")

    return ipaes

def process_pdb_json_pair(pdb_file_path, json_file_path, output_file_path):
    """Process a single PDB-JSON pair to extract interfaces and calculate iPAE scores."""
    logging.info(f"Processing PDB file {pdb_file_path} and JSON file {json_file_path}.")

    try:
        # 1) Load the PAE matrix
        pae_matrix = load_pae(json_file_path)
        if pae_matrix is None:
            return None

        # 2) Load chain coordinates
        chain_a_coords = load_chain_coords_custom(pdb_file_path, 'A')
        chain_b_coords = load_chain_coords_custom(pdb_file_path, 'B')

        # 3) Find interfaces
        interfaces = find_multiple_interfaces(chain_a_coords, chain_b_coords)

        # 4) Calculate iPAE scores
        ipae_scores = calculate_weighted_ipae(
            interfaces,
            pae_matrix,
            len(chain_a_coords),
            len(chain_b_coords)
        )

        # 5) Build residue dictionaries
        chain_a_residues = {res['index']: res['residue'] for res in chain_a_coords}
        chain_b_residues = {res['index']: res['residue'] for res in chain_b_coords}

        # 6) Write main summary output
        with open(output_file_path, 'w') as output_file:
            for idx, ipae in enumerate(ipae_scores):
                confidence = (
                    "High confidence in interaction: likely biologically relevant." if ipae < 5 else
                    "Low confidence in interaction: may not be biologically relevant." if ipae > 10 else
                    "Moderate confidence in interaction."
                )
                output_file.write(f"Interface {idx}:\n")
                output_file.write(f"  iPAE Score: {ipae:.2f} ({confidence})\n")

                # Extract the residue indices for this interface
                chain_a_interface_indices = sorted({r['res_a']['index'] for r in interfaces[idx]})
                chain_b_interface_indices = sorted({r['res_b']['index'] for r in interfaces[idx]})

                # ----- A) Build local & global ranges for Chain A -----
                local_ranges_a, sequences_a = find_continuous_ranges_and_sequences(
                    chain_a_interface_indices,
                    chain_a_residues
                )
                global_range_a = (min(chain_a_interface_indices), max(chain_a_interface_indices))

                # Use .get(...) for safe fallback:
                global_sequence_a = ''.join([
                    three_to_one.get(chain_a_residues.get(i, 'X'), 'X')
                    for i in range(global_range_a[0], global_range_a[1] + 1)
                ])

                # ----- B) Build local & global ranges for Chain B -----
                local_ranges_b, sequences_b = find_continuous_ranges_and_sequences(
                    chain_b_interface_indices,
                    chain_b_residues
                )
                global_range_b = (min(chain_b_interface_indices), max(chain_b_interface_indices))

                # Same safe fallback for Chain B
                global_sequence_b = ''.join([
                    three_to_one.get(chain_b_residues.get(i, 'X'), 'X')
                    for i in range(global_range_b[0], global_range_b[1] + 1)
                ])

                # Write results for Chain A
                output_file.write("  Chain A Local Interface Coordinates (Start, End, Residue Sequence):\n")
                for r, seq in zip(local_ranges_a, sequences_a):
                    output_file.write(f"    {r[0]}, {r[1]}, \"{seq}\"\n")
                output_file.write(
                    f"  Chain A Global Interface Coordinates (Start, End, Residue Sequence): "
                    f"{global_range_a[0]}, {global_range_a[1]}, \"{global_sequence_a}\"\n"
                )

                # Write results for Chain B
                output_file.write("  Chain B Local Interface Coordinates (Start, End, Residue Sequence):\n")
                for r, seq in zip(local_ranges_b, sequences_b):
                    output_file.write(f"    {r[0]}, {r[1]}, \"{seq}\"\n")
                output_file.write(
                    f"  Chain B Global Interface Coordinates (Start, End, Residue Sequence): "
                    f"{global_range_b[0]}, {global_range_b[1]}, \"{global_sequence_b}\"\n"
                )

        # 7) Generate a unique filename key
        base_name = os.path.basename(pdb_file_path).split('_')
        key = "_".join([base_name[i] for i in [0, 1, 2, 3, 5, 6] if i < len(base_name)])

        # 8) Write additional per-residue file
        additional_output_file_path = os.path.join(
            os.path.dirname(pdb_file_path),
            f"{key}_residue_pair_details.txt"
        )
        with open(additional_output_file_path, 'w') as additional_file:
            for idx, interface in enumerate(interfaces):
                additional_file.write(f"Interface {idx}:\n")
                additional_file.write("  Per-Residue Distances and PAE Values:\n")
                for pair in interface:
                    try:
                        res_a = pair['res_a']
                        res_b = pair['res_b']
                        distance = np.linalg.norm(res_a['coord'] - res_b['coord'])
                        res_a_code = three_to_one.get(chain_a_residues.get(res_a['index'], 'X'), 'X')
                        res_b_code = three_to_one.get(chain_b_residues.get(res_b['index'], 'X'), 'X')

                        index_a = res_a['index'] - 1
                        index_b = (res_b['index'] - 1) + len(chain_a_coords)

                        if 0 <= index_a < pae_matrix.shape[0] and 0 <= index_b < pae_matrix.shape[1]:
                            pae_value = min(pae_matrix[index_a, index_b], pae_matrix[index_b, index_a])
                        else:
                            logging.warning(
                                f"PAE index out of bounds: A={res_a['index']}, B={res_b['index']}."
                            )
                            pae_value = None

                        if pae_value is not None:
                            additional_file.write(
                                f"    Residue A (index {res_a['index']}, {res_a_code}) - "
                                f"Residue B (index {res_b['index']}, {res_b_code}): "
                                f"Distance: {distance:.2f} Å, PAE: {pae_value:.2f}\n"
                            )
                        else:
                            additional_file.write(
                                f"    Residue A (index {res_a['index']}, {res_a_code}) - "
                                f"Residue B (index {res_b['index']}, {res_b_code}): "
                                f"Distance: {distance:.2f} Å, PAE: N/A\n"
                            )
                    except Exception as e:
                        logging.error(f"Error writing residue pair data for interface {idx}: {e}")



        logging.info(f"Finished processing {pdb_file_path}. "
                     f"Results written to {output_file_path} and {additional_output_file_path}.")
        return os.path.basename(pdb_file_path), ipae_scores

    except Exception as e:
        logging.error(f"Error processing {pdb_file_path} and {json_file_path}: {e}")
        return None
        
### added to let the user decide if he wants to do parallel processing
def ask_for_parallel_processing():
    """Prompt the user to choose whether to use parallel processing."""
    while True:
        choice = input("Do you want to use parallel processing? (yes/no): ").strip().lower()
        if choice == 'yes':
            return True
        elif choice == 'no':
            return False
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")

# Function renamed from process_directory_parallel 
def process_directory(directory_path, use_parallel=True):
    pdb_files = []
    json_files = []
    logging.info(f"Scanning directory {directory_path} for PDB and JSON files.")

    for root, _, files in os.walk(directory_path):
        for file in files:
            if file.endswith(".pdb"):
                pdb_files.append(os.path.join(root, file))
            elif file.endswith(".json"):
                json_files.append(os.path.join(root, file))

    def get_relevant_parts(filename):
        parts = filename.split("_")
        return "_".join([parts[i] for i in [0, 1, 2, 3, 5, 6] if i < len(parts)])

    pdb_dict = {get_relevant_parts(os.path.basename(f)): f for f in pdb_files}
    json_dict = {get_relevant_parts(os.path.basename(f)): f for f in json_files}

    pairs_to_process = [
        (
            pdb_dict[key],
            json_dict[key],
            os.path.join(os.path.dirname(pdb_dict[key]), f"{key}_interaction_results.txt")
        )
        for key in pdb_dict if key in json_dict
    ]

    logging.info(f"Identified {len(pairs_to_process)} PDB-JSON pairs to process.")

    summary_results = []
    if use_parallel:
        with ProcessPoolExecutor(max_workers=4) as executor:
            futures = {executor.submit(process_pdb_json_pair, *pair): pair for pair in pairs_to_process}
            for future in as_completed(futures):
                result = future.result()
                if result:
                    summary_results.append(result)


    else:
        for pair in pairs_to_process:
            result = process_pdb_json_pair(*pair)
            if result:
                summary_results.append(result)

    if summary_results:
        summary_output_path = os.path.join(directory_path, "summary_results.csv")
        with open(summary_output_path, 'w') as summary_file:
            max_interfaces = max(len(r[1]) for r in summary_results)
            header = "ID," + ",".join([f"iPAE Interface {i}" for i in range(max_interfaces)])
            summary_file.write(header + "\n")

            sorted_summary_results = sorted(summary_results, key=lambda x: x[0])
            for res in sorted_summary_results:
                key, ipae_scores = res
                summary_file.write(",".join(map(str, [key] + ipae_scores)) + "\n")

        logging.info(f"Summary results have been written to {summary_output_path}")
    else:
        logging.warning("No PDB-JSON pairs were processed successfully. No summary file was generated.")

if __name__ == "__main__":
    main_directory_path = './pdbs'
   
   
    # Ask the user if they want to use parallel processing
    use_parallel = ask_for_parallel_processing()
    # Process the directory accordingly
    process_directory(main_directory_path, use_parallel)

 
