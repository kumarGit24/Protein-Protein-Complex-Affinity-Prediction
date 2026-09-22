import os
import subprocess
import pandas as pd
import multiprocessing
import shutil

# ==========================================
# 1. CONFIGURATION
# ==========================================
PDB_DIR = "enter_the_name_of_cleaned_PDBs_folder/"
OUTPUT_TSV = "hbplus_output_file_name.tsv"
CHAINS_FILE = "antibody_antigen_chains.txt"

# Temporary directory for HBPLUS outputs (.hb2 files)
TMP_DIR = "hbplus_tmp/"
os.makedirs(TMP_DIR, exist_ok=True)

# ==========================================
# 2. WORKER FUNCTION
# ==========================================
def run_hbplus(task_info):
    pdb_path, pdb_id, ab_chains, ag_chains = task_info
    
    # Copy the PDB to the TMP_DIR to keep the working directory clean
    # HBPLUS names its output based on the input filename (e.g., target.pdb -> target.hb2)
    tmp_pdb = os.path.join(TMP_DIR, f"{pdb_id}.pdb")
    hb2_file = os.path.join(TMP_DIR, f"{pdb_id}.hb2")
    
    try:
        shutil.copy(pdb_path, tmp_pdb)
        
        # Run HBPLUS (use cwd=TMP_DIR so the .hb2 file is generated inside the temp folder)
        subprocess.run(["hbplus", f"{pdb_id}.pdb"], cwd=TMP_DIR, capture_output=True, text=True, check=True)
        
        if not os.path.exists(hb2_file):
            print(f"  - ? {pdb_id}: HBPLUS ran but did not generate a .hb2 file.")
            return None
            
        interface_hbond_count = 0
        
        # --- PARSE THE .HB2 FILE ---
        # HBPLUS .hb2 files have a strict fixed-width format. 
        # Data lines look like this: "A0011-VAL N   B0007-GLY O    3.02..."
        # Index 0 is the Donor Chain. Index 14 is the Acceptor Chain.
        # Hyphens are always at index 5 and 19 for valid data lines.
        with open(hb2_file, 'r') as f:
            for line in f:
                if len(line) > 20 and line[5] == '-' and line[19] == '-':
                    donor_chain = line[0]
                    acceptor_chain = line[14]
                    
                    # Check if the H-bond crosses the interface
                    # Scenario 1: Donor is Antibody, Acceptor is Antigen
                    if donor_chain in ab_chains and acceptor_chain in ag_chains:
                        interface_hbond_count += 1
                    # Scenario 2: Donor is Antigen, Acceptor is Antibody
                    elif donor_chain in ag_chains and acceptor_chain in ab_chains:
                        interface_hbond_count += 1
                        
        print(f"{pdb_id}: Found {interface_hbond_count} interface hydrogen bonds.")
        
        # Clean up the large output files to save hard drive space
        os.remove(tmp_pdb)
        os.remove(hb2_file)
        
        return {
            "PDB_ID": pdb_id,
            "Interface_H_Bonds": interface_hbond_count
        }

    except subprocess.CalledProcessError as e:
        print(f"{pdb_id}: HBPLUS crashed. Make sure 'hbplus' is in your system PATH.")
        return None
    except FileNotFoundError:
        print(f"Command Not Found: The 'hbplus' executable was not found.")
        return None
    except Exception as e:
        print(f"{pdb_id}: Unexpected error: {e}")
        return None

# ==========================================
# 3. PARALLEL EXECUTION
# ==========================================
def main():
    if not os.path.exists(CHAINS_FILE):
        print(f"Error: {CHAINS_FILE} not found.")
        return

    chains_df = pd.read_csv(CHAINS_FILE, sep='\s+')
    
    tasks = []
    for filename in os.listdir(PDB_DIR):
        if not filename.endswith(".pdb"): continue
        
        pdb_id = filename[:4].lower()
        match = chains_df[chains_df['PDB_ID'].str.lower() == pdb_id]
        
        if not match.empty:
            # Treat strings like "HL" as a list of characters ['H', 'L'] automatically
            # This makes the "chain in ab_chains" 
            ab_chains = str(match.iloc[0]['Antibody']).strip()
            ag_chains = str(match.iloc[0]['Protein']).strip()
            
            pdb_path = os.path.join(PDB_DIR, filename)
            tasks.append((pdb_path, pdb_id, ab_chains, ag_chains))

    num_processes = max(1, multiprocessing.cpu_count() - 2)
    print(f"\nStarting HBPLUS calculations on {len(tasks)} structures using {num_processes} cores...\n")

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = [res for res in pool.imap_unordered(run_hbplus, tasks) if res is not None]

    if results:
        df = pd.DataFrame(results)
        df = df.sort_values(by="Interface_H_Bonds", ascending=False)
        df.to_csv(OUTPUT_TSV, sep='\t', index=False)
        print(f"\n?? Complete! Interface H-bond counts saved to '{OUTPUT_TSV}'")
        
        # Remove the temporary directory if it's empty
        try:
            os.rmdir(TMP_DIR)
        except OSError:
            pass

if __name__ == "__main__":
    main()
