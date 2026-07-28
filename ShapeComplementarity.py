import os
import subprocess
import pandas as pd
import multiprocessing

# ==========================================
# 1. CONFIGURATION
# ==========================================
PDB_DIR = "enter_the_name_of_cleaned_PDBs_folder/"
OUTPUT_TSV = "rosetta_sc_output_file_name.tsv"
SCORE_DIR = "rosetta_sc_scores/"
CHAINS_FILE = "antibody_antigen_chains.txt"

ROSETTA_EXEC = "rosetta_scripts_jd3.cxx11threadserialization.linuxclangrelease" 
XML_FILE = "sc_packstat.xml"

os.makedirs(SCORE_DIR, exist_ok=True)

# ==========================================
# 2. GENERATE ROSETTASCRIPTS XML
# ==========================================
# Add packstat="true" to the InterfaceAnalyzer. 
xml_content = """<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="ref15" weights="ref2015.wts"/>
    </SCOREFXNS>
    <MOVERS>
        <InterfaceAnalyzerMover name="analyze_interface" scorefxn="ref15" packstat="true" pack_input="false" pack_separated="false" interface="%%interface_str%%"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover="analyze_interface"/>
    </PROTOCOLS>
</ROSETTASCRIPTS>
"""

with open(XML_FILE, "w") as f:
    f.write(xml_content)

# ==========================================
# 3. WORKER FUNCTION
# ==========================================
def run_rosetta_metrics(task_info):
    pdb_path, pdb_id, interface_str = task_info
    score_file = os.path.join(SCORE_DIR, f"{pdb_id}.sc")
    
    command = [
        ROSETTA_EXEC,
        "-s", pdb_path,
        "-parser:protocol", XML_FILE,
        "-parser:script_vars", f"interface_str={interface_str}",
        "-out:file:scorefile", score_file,
        "-jd3:nthreads", "1", 
        "-mute", "all" 
    ]
    
    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
        
        if os.path.exists(score_file):
            # --- CUSTOM ROSETTA SCORE PARSER ---
            with open(score_file, 'r') as f:
                lines = f.readlines()
            
            headers, data_values = [], []
            for line in lines:
                if line.startswith('SCORE:'):
                    parts = line.strip().split()
                    if len(parts) > 1 and parts[1] == 'total_score':
                        headers = parts
                    elif len(parts) > 1 and parts[1] != 'total_score':
                        data_values = parts
            
            if not headers or not data_values:
                print(f"  - {pdb_id}: Score file missing data rows.")
                return None
                
            score_dict = dict(zip(headers, data_values))
            
            # Extract Shape Complementarity and PackStat
            sc_value = float(score_dict.get('sc_value', 0.0))
            packstat = float(score_dict.get('packstat', 0.0))
            
            print(f"  -{pdb_id}: Sc = {sc_value:.3f} | PackStat = {packstat:.3f}")
            
            return {
                "PDB_ID": pdb_id,
                "Interface": interface_str,
                "Shape_Complementarity_Sc": sc_value,
                "Packing_Density_PackStat": packstat
            }
        else:
            print(f"  - {pdb_id}: No score file created.")
            return None
            
    except subprocess.CalledProcessError as e:
        print(f"  - Rosetta Crash on {pdb_id}: {e.stderr.strip()[-100:]}")
        return None

# ==========================================
# 4. PARALLEL EXECUTION
# ==========================================
def main():
    if not os.path.exists(CHAINS_FILE):
        print(f"Error: {CHAINS_FILE} not found.")
        return
        
    chains_df = pd.read_csv(CHAINS_FILE, sep=r'\s+')
    
    # Create the lookup dictionary mapping PDB to its interface (e.g., 'A_B')
    chain_map = {str(row['PDB_ID']).strip().lower(): f"{str(row['Antibody']).strip()}_{str(row['Protein']).strip()}" 
                 for _, row in chains_df.iterrows()}

    tasks = []
    for filename in os.listdir(PDB_DIR):
        if filename.endswith(".pdb"):
            pdb_id = filename[:4].lower() 
            if pdb_id in chain_map:
                tasks.append((os.path.join(PDB_DIR, filename), pdb_id, chain_map[pdb_id]))

    num_processes = max(1, multiprocessing.cpu_count() - 2)
    print(f"\nCalculating Sc and PackStat on {len(tasks)} files using {num_processes} cores...\n")

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = [res for res in pool.imap_unordered(run_rosetta_metrics, tasks) if res is not None]

    if results:
        df = pd.DataFrame(results)
        # Sort by best shape complementarity (highest is best)
        df = df.sort_values(by="Shape_Complementarity_Sc", ascending=False)
        df.to_csv(OUTPUT_TSV, sep='\t', index=False)
        print(f"\n Complete! Saved to '{OUTPUT_TSV}'")

if __name__ == "__main__":
    main()
