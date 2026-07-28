import csv

input_file = "antibody_antigen_chains.txt"   # tab-separated file with PDB_ID, Antibody, Protein
output_file = "script_chimerax.cxc" # ChimeraX script

input_folder = "enter_the_name_of_input_PDBs_directory"
output_folder = "enter_the_name_of_output_directory"

lines = []

# --- Process each PDB entry ---
with open(input_file, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        pdb_id = row["#PDB_ID"].strip()
        ab_chain = ",".join(list(row["Antibody"].strip()))   # HL -> H,L
        prot_chain = ",".join(list(row["Protein"].strip()))  # BC -> B,C
        chains = ",".join([ab_chain, prot_chain]) if ab_chain and prot_chain else ab_chain or prot_chain

        lines.append(f"# --- Process PDB: {pdb_id} ---")
        lines.append(f"echo Processing {pdb_id}...")
        lines.append(f"open {input_folder}/{pdb_id}_complex_cleaned.pdb")
        lines.append(f"select #1/{chains}")
        lines.append(f"save {output_folder}/{pdb_id}_complex_cleaned.pdb selectedOnly true models #1")
        lines.append(f"close #1\n")

# --- Save script ---
with open(output_file, "w") as f:
    f.write("\n".join(lines))

print(f"ChimeraX script saved as {output_file}")

