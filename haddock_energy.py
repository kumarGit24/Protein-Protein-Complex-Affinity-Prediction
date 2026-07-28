import os
import subprocess
import re
import pandas as pd
import multiprocessing

def process_single_pdb(pdb_path):
    """
    Worker function to run PRODIGY and extract a list of contact dictionaries.
    """
    filename = os.path.basename(pdb_path)
    base_name = filename[:-4]
    parts = base_name.split('_')
    
    if len(parts) < 3:
        return None # Skip files with wrong naming convention
    
    chains1_formatted = ",".join(list(parts[1]))
    chains2_formatted = ",".join(list(parts[2]))

    # --- Run PRODIGY ---
    command = [
        "prodigy",
        pdb_path,
        "--selection", chains1_formatted, chains2_formatted,
        "--contact_list",
        "-q"
    ]
    
    try:
        process = subprocess.run(command, capture_output=True, text=True, check=True)
        output = process.stdout.strip()

        # --- Parse Contacts ---
        # PRODIGY contact line format: ResName1 ResNum1 Chain1 ResName2 ResNum2 Chain2
        # Example: ASN 31 A TRP 107 B
        contacts = []
        for line in output.split('\n'):
            parts = line.split()
            if len(parts) == 6:
                contacts.append({
                    "PDB_File": filename,
                    "Res_1": parts[0],
                    "Idx_1": parts[1],
                    "Chain_1": parts[2],
                    "Res_2": parts[3],
                    "Idx_2": parts[4],
                    "Chain_2": parts[5]
                })
        
        if contacts:
            print(f"  {filename}: Found {len(contacts)} contacts.")
            return contacts
        else:
            return None

    except subprocess.CalledProcessError:
        print(f"  Error processing {filename}")
        return None

def run_contact_extraction(pdb_directory, output_tsv="all_contacts.tsv"):
    if not os.path.isdir(pdb_directory):
        print(f"Error: Directory '{pdb_directory}' not found.")
        return

    all_pdb_paths = [os.path.join(pdb_directory, f) for f in os.listdir(pdb_directory) if f.endswith(".pdb")]

    num_processes = 16 
    print(f"Starting parallel contact extraction on {len(all_pdb_paths)} files...\n")

    with multiprocessing.Pool(processes=num_processes) as pool:
        # result_lists will be a list of lists (each sub-list contains contacts for one PDB)
        result_lists = list(pool.imap_unordered(process_single_pdb, all_pdb_paths))

    # Flatten the list of lists into a single list of dictionaries
    flattened_contacts = [contact for sublist in result_lists if sublist for contact in sublist]

    if flattened_contacts:
        df = pd.DataFrame(flattened_contacts)
        
        # Save to TSV (sep='\t')
        df.to_csv(output_tsv, sep='\t', index=False)
        print(f"\n Complete! {len(flattened_contacts)} total contacts saved to '{output_tsv}'")
    else:
        print("\nNo contacts were extracted.")

if __name__ == "__main__":
    # Ensure you point to the folder containing the CLEANED PDBs
    pdb_folder = "FinalFix/" 
    run_contact_extraction(pdb_folder)
