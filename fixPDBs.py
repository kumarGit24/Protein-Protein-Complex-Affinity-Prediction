import os
import io
import multiprocessing
from collections import defaultdict
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Directory name with PDB files 
input_folder = "give_the_name_of_directory/"

# Directory where the final, correctly-chained PDBs will be saved.
output_folder = "give_the_name_of_directory/"

# All the PDBs should be arranged in the order: PDB_ID Antibody_Chain(s) Protein_Chain(s)
mapping_file = "antibody_antigen_chains.txt"


def create_pdb_chain_map(filepath):
    """Reads the mapping file and creates a dictionary.
    
    Returns:
        A dictionary like {'1fc2': ['C', 'D'], '2mta': ['H', 'L', 'A']}
    """
    pdb_map = {}
    try:
        with open(filepath, 'r') as f:
            next(f) # Skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                pdb_id = parts[0].lower() # Use lowercase for consistent matching
                
                chains = []
                # Expand multi-character chain strings like "HL" into "H", "L"
                for chain_str in parts[1:]:
                    chains.extend(list(chain_str))
                pdb_map[pdb_id] = chains
    except FileNotFoundError:
        print(f"FATAL ERROR: Mapping file not found at '{filepath}'")
        return None
    return pdb_map


def clean_pdb_file(pdb_path, output_path, chains_to_keep):
    """
    Splits a PDB by chain, fixes each one individually, and recombines them
    with their original chain IDs preserved. This is worker function.
    """
    # This print statement helps track progress in parallel execution
    process_id = multiprocessing.current_process().pid
    print(f"[PID {process_id}] Processing '{os.path.basename(pdb_path)}'")
    
    if not chains_to_keep:
        print(f"  -> WARNING: No chains specified for {os.path.basename(pdb_path)}. Skipping.")
        return

    # Step 1: Read the original PDB and split ATOM records by chain ID
    chain_atom_lines = defaultdict(list)
    try:
        with open(pdb_path, 'r') as f_in:
            for line in f_in:
                if line.startswith("ATOM"):
                    chain_id = line[21]
                    if chain_id in chains_to_keep:
                        chain_atom_lines[chain_id].append(line)
    except FileNotFoundError:
        print(f"  -> ERROR: Source file not found at {pdb_path}.")
        return
    except Exception as e:
        print(f"  -> ERROR: Could not read {os.path.basename(pdb_path)}. Reason: {e}")
        return

    final_pdb_text = ""
    # Use process and chain ID for unique temp file names
    temp_files_to_delete = []

    # Step 2: Process each chain individually
    for original_chain_id in chains_to_keep:
        if original_chain_id not in chain_atom_lines:
            print(f"  -> WARNING: Chain '{original_chain_id}' not found in {os.path.basename(pdb_path)}.")
            continue
        
        temp_pdb_path = f"temp_{process_id}_{original_chain_id}.pdb"
        temp_files_to_delete.append(temp_pdb_path)
        
        with open(temp_pdb_path, 'w') as f_temp:
            f_temp.writelines(chain_atom_lines[original_chain_id])
            f_temp.write("END\n")

        try:
            fixer = PDBFixer(filename=temp_pdb_path)
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.addMissingHydrogens(7.0)

            buffer = io.StringIO()
            PDBFile.writeFile(fixer.topology, fixer.positions, buffer)
            buffer.seek(0)

            for line in buffer:
                if line.startswith('ATOM'):
                    final_pdb_text += line[:21] + original_chain_id + line[22:]
            final_pdb_text += "TER\n"

        except Exception as e:
            print(f"  -> ERROR: PDBFixer failed for chain '{original_chain_id}' in {os.path.basename(pdb_path)}. Reason: {e}")
            continue
        finally:
            for temp_file in temp_files_to_delete:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            temp_files_to_delete.clear()
            
    # Step 3: Write the combined PDB file
    if final_pdb_text:
        with open(output_path, 'w') as f_out:
            f_out.write(final_pdb_text)
            f_out.write("END\n")
        # Quieter success message for parallel execution
        print(f"  -> Successfully cleaned and saved to {output_path}")
    else:
        print(f"  -> ERROR: No valid chains processed for {os.path.basename(pdb_path)}.")

# --- Parallel Execution ---
if __name__ == "__main__":
    # Set the number of CPU cores to use.
    # Or use multiprocessing.cpu_count() to use all available cores.
    num_processes = 16 

    # Ensure the output directory exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created output directory: {output_folder}")

    # Step 1: Create the PDB ID -> chain mapping ONCE in the main process
    pdb_chain_map = create_pdb_chain_map(mapping_file)
    if not pdb_chain_map:
        exit()

    # Step 2: Create a list of tasks for the process pool
    tasks = []
    all_pdb_files = [f for f in os.listdir(input_folder) if f.endswith(".pdb")]
    
    if not all_pdb_files:
        print(f"No .pdb files found in '{input_folder}'.")
    else:
        for filename in all_pdb_files:
            pdb_id = filename.split('_')[0].lower()
            if pdb_id in pdb_chain_map:
                # For each file, create a task tuple with ALL required arguments
                input_path = os.path.join(input_folder, filename)
                output_path = os.path.join(output_folder, filename)
                chains_to_keep = pdb_chain_map[pdb_id]
                tasks.append((input_path, output_path, chains_to_keep))
            else:
                print(f"Skipping '{filename}': PDB ID '{pdb_id}' not found in '{mapping_file}'.")

        print(f"\nStarting parallel cleaning for {len(tasks)} files using {num_processes} processes...\n")

        # Step 3: Use a multiprocessing Pool to run jobs in parallel
        with multiprocessing.Pool(processes=num_processes) as pool:
            # 'starmap' distributes the tasks to the worker function.
            # It automatically unpacks each 3-element tuple in 'tasks' 
            # into the three arguments for clean_pdb_file.
            pool.starmap(clean_pdb_file, tasks)

        print("\nBatch processing complete.")
