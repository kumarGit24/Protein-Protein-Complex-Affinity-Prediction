import os
import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from functools import partial

# ==========================================
# 1. CONFIGURATION
# ==========================================
FOLDX_PATH = "path_to_foldx_installation/foldx_20270131" 
PDB_DIR = "enter_the_name_of_cleaned_PDBs_folder/" 
OUTPUT_DIR = "foldx_results/"
CHAINS_FILE = "antibody_antigen_chains.txt"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================================
# 2. LOAD CHAIN MAPPINGS
# ==========================================
chains_df = pd.read_csv(CHAINS_FILE, sep='\s+', comment=None)
chains_df.columns = ['PDB_ID', 'Antibody_Chains', 'Protein_Chains']

chain_map = {}
for index, row in chains_df.iterrows():
    # Store as joined strings (e.g., "AB" and "G") for FoldX
    chain_map[row['PDB_ID'].lower()] = {
        'ab': str(row['Antibody_Chains']),
        'ag': str(row['Protein_Chains'])
    }

print(f"Loaded mappings for {len(chain_map)} complexes.")

# ==========================================
# 3. DEFINE THE FOLDX EXECUTION FUNCTION
# ==========================================
def run_foldx_analyse_complex(pdb_info):
    """
    Runs FoldX AnalyseComplex on a single PDB file.
    pdb_info is a tuple: (pdb_id, ab_str, ag_str)
    """
    pdb_id, ab_str, ag_str = pdb_info
    pdb_filename = f"{pdb_id}.pdb"
    pdb_path = os.path.join(PDB_DIR, pdb_filename)
    
    if not os.path.exists(pdb_path):
        return False, f"File not found: {pdb_filename}"

    # FoldX expects groups separated by a comma: e.g., "AB,G"
    complex_chains = f"{ab_str},{ag_str}"

    command = [
        FOLDX_PATH,
        "-c", "AnalyseComplex",
        "--pdb", pdb_filename,
        "--analyseComplexChains", complex_chains,
        "--pdb-dir", PDB_DIR,
        "--output-dir", OUTPUT_DIR
    ]
    
    try:
        # capture_output=True to keep the console clean during parallel runs
        subprocess.run(command, capture_output=True, check=True, text=True)
        return True, pdb_id
    except subprocess.CalledProcessError as e:
        return False, f"{pdb_filename} failed: {e.stderr}"

# ==========================================
# 4. RUN FOLDX IN PARALLEL
# ==========================================
print("Starting FoldX calculations...")

# Prepare list of tuples for the executor
tasks = [(pid, info['ab'], info['ag']) for pid, info in chain_map.items()]

with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
    # Use map to distribute the tasks list
    results = list(executor.map(run_foldx_analyse_complex, tasks))

for success, msg in results:
    if not success:
        print(f"Warning: {msg}")

print("FoldX runs completed. Parsing targeted interactions...")

# ==========================================
# 5. PARSE RESULTS
# ==========================================
parsed_data = []

for pdb_id, chains in chain_map.items():
    ab_list = list(chains['ab']) # For membership checking
    ag_list = list(chains['ag'])
    
    summary_file = os.path.join(OUTPUT_DIR, f"Summary_{pdb_id}_AC.fxout")
    
    if not os.path.exists(summary_file):
        continue
        
    complex_thermo = {
        'PDB_ID': pdb_id,
        'Interaction_Energy_kcal_mol': 0.0,
        'Backbone_Hbond': 0.0,
        'Sidechain_Hbond': 0.0,
        'Van_der_Waals': 0.0,
        'Electrostatics': 0.0,
        'Solvation_Polar': 0.0,
        'Solvation_Hydrophobic': 0.0,
        'Van_der_Waals_clashes': 0.0,
        'Entropy_Sidechain': 0.0,
        'Entropy_Mainchain': 0.0
    }
    
    pairs_found = 0
    
    with open(summary_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            # Standard FoldX output lines start with the PDB filename
            if line.startswith(pdb_id) or f"{pdb_id}.pdb" in line:
                parts = line.strip().split('\t')
                if len(parts) < 15: continue
                
                group1 = parts[1]
                group2 = parts[2]
                
                # Verify if this specific interaction row is the one we want
                is_ab_ag = (group1 == chains['ab'] and group2 == chains['ag'])
                is_ag_ab = (group1 == chains['ag'] and group2 == chains['ab'])
                
                if is_ab_ag or is_ag_ab:
                    pairs_found += 1
                    complex_thermo['Interaction_Energy_kcal_mol'] += float(parts[5])
                    complex_thermo['Backbone_Hbond'] += float(parts[6])
                    complex_thermo['Sidechain_Hbond'] += float(parts[7])
                    complex_thermo['Van_der_Waals'] += float(parts[8])
                    complex_thermo['Electrostatics'] += float(parts[9])
                    complex_thermo['Solvation_Polar'] += float(parts[10])
                    complex_thermo['Solvation_Hydrophobic'] += float(parts[11])
                    complex_thermo['Van_der_Waals_clashes'] += float(parts[12])
                    complex_thermo['Entropy_Sidechain'] += float(parts[13])
                    complex_thermo['Entropy_Mainchain'] += float(parts[14])
                    
    if pairs_found > 0:
        parsed_data.append(complex_thermo)

# ==========================================
# 6. EXPORT
# ==========================================
df_thermo = pd.DataFrame(parsed_data)
output_csv = "name_of_the_FoldX_output_file.csv"
df_thermo.to_csv(output_csv, index=False)

print(f"\nSuccess! Targeted properties extracted for {len(df_thermo)} complexes.")
