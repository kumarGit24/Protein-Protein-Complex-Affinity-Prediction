import os
import subprocess
import pandas as pd
import multiprocessing

# ==========================================
# 1. CONFIGURATION
# ==========================================
PDB_DIR = "enter_the_name_of_cleaned_PDBs_folder/"
OUTPUT_TSV = "rosetta_output_file_name.tsv"
SCORE_DIR = "rosetta_scores/"
CHAINS_FILE = "antibody_antigen_chains.txt"

ROSETTA_EXEC = "rosetta_scripts_jd3.cxx11threadserialization.linuxclangrelease" 
XML_FILE = "relax_and_analyze.xml"

os.makedirs(SCORE_DIR, exist_ok=True)

# ==========================================
# 2. GENERATE ROSETTASCRIPTS XML (WITH FASTRELAX)
# ==========================================
# Add FastRelax to the MOVERS and run it before the InterfaceAnalyzer
xml_content = """<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="ref15" weights="ref2015.wts"/>
    </SCOREFXNS>
    <MOVERS>
        <FastRelax name="relax_complex" scorefxn="ref15" repeats="1"/>
        
        <InterfaceAnalyzerMover name="analyze_interface" scorefxn="ref15" pack_input="true" pack_separated="true" interface="%%interface_str%%"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover="relax_complex"/>
        <Add mover="analyze_interface"/>
    </PROTOCOLS>
</ROSETTASCRIPTS>
"""

with open(XML_FILE, "w") as f:
    f.write(xml_content)

# ==========================================
# 3. WORKER FUNCTION
# ==========================================
def run_rosetta_interface(task_info):
    pdb_path, pdb_id, interface_str = task_info
    filename = os.path.basename(pdb_path)
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
            try:
                # --- CUSTOM ROSETTA SCORE PARSER ---
                with open(score_file, 'r') as f:
                    lines = f.readlines()
                
                headers = []
                data_values = []
                
                for line in lines:
                    if line.startswith('SCORE:'):
                        parts = line.strip().split()
                        # If the second word is 'total_score', this is a header line
                        if len(parts) > 1 and parts[1] == 'total_score':
                            headers = parts
                        # Otherwise, it is a data line. We keep overwriting to get the final output row.
                        elif len(parts) > 1 and parts[1] != 'total_score':
                            data_values = parts
                
                if not headers or not data_values:
                    print(f"  - Parse Error on {pdb_id}: Score file is missing valid data rows.")
                    return None
                    
                # Map the headers to the data values
                score_dict = dict(zip(headers, data_values))
                
                # Extract the metrics safely using .get() with a default of 0.0 if missing
                dg_separated = float(score_dict.get('dG_separated', 0.0))
                dsasa = float(score_dict.get('dSASA_int', 0.0))
                vdw_attractive = float(score_dict.get('fa_atr', 0.0))
                vdw_repulsive = float(score_dict.get('fa_rep', 0.0))
                electrostatics = float(score_dict.get('fa_elec', 0.0))
                # ---------------------------------------
                
                print(f"  - Success: {pdb_id} | dG: {dg_separated:.2f} | vdW(atr): {vdw_attractive:.2f}")
                
                return {
                    "PDB_ID": pdb_id,
                    "PDB_File": filename,
                    "Interface": interface_str,
                    "Total_dG_Binding": dg_separated,
                    "dSASA": dsasa,
                    "vdW_Attractive_Energy": vdw_attractive,
                    "vdW_Repulsive_Energy": vdw_repulsive,
                    "Electrostatic_Energy": electrostatics
                }
            except Exception as e:
                print(f"  - Parse Error on {pdb_id}: {e}")
        else:
            print(f"  - File Error on {pdb_id}. No score file created.")
            
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
        
    chains_df = pd.read_csv(CHAINS_FILE, sep='\s+')
    chains_df['Antibody'] = chains_df['Antibody'].astype(str)
    chains_df['Protein'] = chains_df['Protein'].astype(str)
    
    chain_map = {str(row['PDB_ID']).strip().lower(): f"{row['Antibody'].strip()}_{row['Protein'].strip()}" 
                 for _, row in chains_df.iterrows()}

    all_pdb_files = [f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")]
    tasks = []
    
    for filename in all_pdb_files:
        pdb_id = filename[:4].lower() 
        if pdb_id in chain_map:
            tasks.append((os.path.join(PDB_DIR, filename), pdb_id, chain_map[pdb_id]))

    num_processes = max(1, multiprocessing.cpu_count() - 2)
    print(f"\nStarting Minimization + Analysis on {len(tasks)} files using {num_processes} cores...\n")

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = list(pool.imap_unordered(run_rosetta_interface, tasks))

    successful_results = [res for res in results if res is not None]

    if successful_results:
        df = pd.DataFrame(successful_results)
        df = df.sort_values(by="Total_dG_Binding", ascending=True)
        df.to_csv(OUTPUT_TSV, sep='\t', index=False)
        print(f"\n Complete! Saved to '{OUTPUT_TSV}'")

if __name__ == "__main__":
    main()
