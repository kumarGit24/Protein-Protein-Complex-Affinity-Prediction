from Bio.PDB import PDBParser
import numpy as np
from glob import glob
import pandas as pd


CONTACT_CUTOFF = [5.5] # Cutoff distances for contacts (angstroms)

pdb_files = glob("enter_the_name_of_cleaned_PDBs_folder/*.pdb") 
df = pd.read_csv("antibody_antigen_chains.txt", delimiter="\t")

parser = PDBParser(QUIET=True)

for pdb_file in pdb_files:
    pdb_id = pdb_file.split("/")[-1].split("_")[0].lower()  # e.g. "1fc2" from "1fc2_cleaned_complex.pdb"

    # Find matching row in df
    row = df[df["#PDB_ID"].str.lower() == pdb_id]
    if row.empty:
        print(f"No entry for {pdb_id} in chains.txt, skipping")
        continue
    row = row.iloc[0]

    # Split multi-letter entries into individual chain IDs (e.g. "HL" -> ["H","L"])
    antibody_chains_expected = list(str(row["Antibody"]).strip())
    protein_chains_expected = list(str(row["Protein"]).strip())

    # Parse PDB structure
    structure = parser.get_structure("complex", pdb_file)

    # Extract actual chains present
    all_chains = [chain.get_id() for chain in structure.get_chains()]
    antibody_chains = [c for c in all_chains if c in antibody_chains_expected]
    protein_chains = [c for c in all_chains if c in protein_chains_expected]

    # Extract atom coordinates
    antibody_atoms = [atom for model in structure for chain in model if chain.get_id() in antibody_chains for residue in chain for atom in residue]
    protein_atoms = [atom for model in structure for chain in model if chain.get_id() in protein_chains for residue in chain for atom in residue]

    # Define residues
    charged_residues = {"ARG", "LYS", "ASP", "GLU", "HIS"} # Charged residues
    ApR = {"ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "MET", "TRP"}  # Apolar residues
    pR = {"ASN", "GLN", "SER", "THR", "TYR", "CYS"}  # Polar residues

    #Extract the aploar and polar residues
    antibody_apolar_atoms = [
        atom for model in structure for chain in model if chain.get_id() in antibody_chains
        for residue in chain if residue.get_resname() in ApR
        for atom in residue
    ]
    protein_apolar_atoms = [
        atom for model in structure for chain in model if chain.get_id() in protein_chains
        for residue in chain if residue.get_resname() in ApR
        for atom in residue
    ]

    #Extract the polar residues
    antibody_polar_atoms = [
        atom for model in structure for chain in model if chain.get_id() in antibody_chains
        for residue in chain if residue.get_resname() in pR
        for atom in residue
    ]
    protein_polar_atoms = [
        atom for model in structure for chain in model if chain.get_id() in protein_chains
        for residue in chain if residue.get_resname() in pR
        for atom in residue
    ]

    # Extract charged residues
    antibody_charged_atoms = [
        atom for model in structure for chain in model if chain.get_id() in antibody_chains
        for residue in chain if residue.get_resname() in charged_residues
        for atom in residue
    ]
    protein_charged_atoms = [
        atom for model in structure for chain in model if chain.get_id() in protein_chains
        for residue in chain if residue.get_resname() in charged_residues
        for atom in residue
    ]

    print(f"Number of charged antibody atoms: {len(antibody_charged_atoms)}")
    print(f"Number of charged protein atoms: {len(protein_charged_atoms)}")
    print(f"Number of apolar antibody atoms: {len(antibody_apolar_atoms)}")
    print(f"Number of apolar protein atoms: {len(protein_apolar_atoms)}")
    print(f"Number of polar antibody atoms: {len(antibody_apolar_atoms)}")  
    print(f"Number of polar protein atoms: {len(protein_apolar_atoms)}")

    # Compute interchain contacts for charged and charged residues
    for cutoff in CONTACT_CUTOFF:
        IC_CC = 0
        for atom1 in antibody_charged_atoms:
            for atom2 in protein_charged_atoms:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    IC_CC += 1
        # Print results for each cutoff
        print(f"Total number of charged-charged interchain contacts (cutoff {cutoff} Å): {IC_CC}")

        # Compute interchain contacts for charged and polar residues
        IC_CP = 0
        for atom1 in antibody_charged_atoms:
            for atom2 in protein_polar_atoms:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    IC_CP += 1
        # Print results for each cutoff
        print(f"Total number of charged-polar interchain contacts (cutoff {cutoff} Å): {IC_CP}")

        # Compute interchain contacts for charged and apolar residues
        IC_CA = 0
        for atom1 in antibody_charged_atoms:
            for atom2 in protein_apolar_atoms:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    IC_CA += 1
        # Print results for each cutoff
        print(f"Total number of charged-apolar interchain contacts (cutoff {cutoff} Å): {IC_CA}")

        # Compute interchain contacts for apolar and apolar residues
        IC_AA = 0
        for atom1 in antibody_apolar_atoms:
            for atom2 in protein_apolar_atoms:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    IC_AA += 1
        # Print results for each cutoff
        print(f"Total number of apolar-apolar interchain contacts (cutoff {cutoff} Å): {IC_AA}")

        # Compute interchain contacts for apolar and polar residues
        IC_AP = 0
        for atom1 in antibody_apolar_atoms:
            for atom2 in protein_polar_atoms:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    IC_AP += 1
        # Print results for each cutoff
        print(f"Total number of apolar-polar interchain contacts (cutoff {cutoff} Å): {IC_AP}")

        # Compute interchain contacts for polar and polar residues
        IC_PP = 0
        for atom1 in antibody_apolar_atoms:
            for atom2 in protein_apolar_atoms:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    IC_PP += 1
        # Print results for each cutoff
        print(f"Total number of polar-polar interchain contacts (cutoff {cutoff} Å): {IC_PP}")

        # Compute interchain contacts for hydrophilic and hydrophilic residues
        IC_hlhl = 0
        for atom1 in antibody_charged_atoms:
            for atom2 in protein_charged_atoms:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    IC_hlhl += 1
        # Print results for each cutoff
        print(f"Total number of hydrophilic-hydrophilic interchain contacts (cutoff {cutoff} Å): {IC_hlhl}")

        # Compute interchain contacts for hydrophilic and hydrophobic residues
        IC_hlhb = 0
        for atom1 in antibody_charged_atoms:
            for atom2 in protein_apolar_atoms:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    IC_hlhb += 1
        # Print results for each cutoff
        print(f"Total number of hydrophilic-hydrophobic interchain contacts (cutoff {cutoff} Å): {IC_hlhb}")

        # Compute interchain contacts for hydrophobic and hydrophobic residues
        IC_hbhb = 0
        for atom1 in antibody_apolar_atoms:
            for atom2 in protein_apolar_atoms:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    IC_hbhb += 1
        # Print results for each cutoff
        print(f"Total number of hydrophobic-hydrophobic interchain contacts (cutoff {cutoff} Å): {IC_hbhb}")

        # save the results to a file
        with open("output_file_name.txt", "a") as f:
            f.write(f"{pdb_file.split('/')[-1]}\t{cutoff}\t{IC_CC}\t{IC_CP}\t{IC_CA}\t{IC_AA}\t{IC_AP}\t{IC_PP}\t{IC_hlhl}\t{IC_hlhb}\t{IC_hbhb}\n")
        print(f"Results saved to interchain_contacts_results.txt")
