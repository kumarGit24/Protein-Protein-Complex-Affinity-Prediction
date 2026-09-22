import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import sasa
import pandas as pd
from glob import glob
from multiprocessing import Pool, cpu_count

# Fuction to calculate the SASA per residue
def calculate_sasa_for_residues(pdb, Rp, residue_list):
    u = mda.Universe(pdb_file)

    # Select residues of interest
    sel = u.select_atoms("resname " + " ".join(residue_list))
    if sel.n_atoms == 0:
        print(f"Warning: No residues found in {residue_list} for probe radius {Rp} nm.")
        return 0.0

    # Run Shrake-Rupley
    sr = sasa.ShrakeRupley(u, probe_radius=Rp * 10.0)  # convert nm → Å
    sr.run()

    # Get per-residue SASA and sum only for selected residues
    residue_sasa = sr.results.residue_sasa[sel.residues.indices].sum()
    return residue_sasa

pdb_files = glob("enter_the_name_of_cleaned_PDBs_folder/*.pdb") 
df = pd.read_csv("antibody_antigen_chains.txt", delimiter="\t")

def process_pdb_file(pdb_file):
    for i in pdb_files:
        pdb_id = i.split("/")[-1].split("_")[0]   # extract e.g. "1fc2" from "P-P/1fc2_complex.pdb"
        row = df[df["#PDB_ID"] == pdb_id].iloc[0]

        # Multi-letter entries into comma separated list
        ab_chain_ids = list(row["Antibody"].strip())   # "HL" -> ["H","L"]
        prot_chain_ids = list(row["Protein"].strip())  # "BC" -> ["B","C"]
        
        #print(f"Processing {pdb_id}: Antibody={ab_chain_ids}, Protein={prot_chain_ids}")

        complex_pdb = md.load(i)

        # residue classifications
        # Cr = ["ARG", "LYS", "HIS", "ASP", "GLU"]
        Cr = ["ARG", "LYS", "HIS", "HSD", "HSE", "HIP", "ASP", "ASH", "GLU", "GLH"]   # Charged residues
        ApR = ["ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "MET", "TRP"]  # Apolar residues
        pR = ["ASN", "GLN", "SER", "THR", "TYR", "CYS"]   # Polar residues

        # Create a lookup for PDB chain IDs 
        pdb_chain_to_internal_id = {}
        for chain in complex_pdb.topology.chains:
            pdb_chain_to_internal_id[chain.chain_id] = chain.index
            print(f"Chain {chain.chain_id} has internal index {chain.index}")

        prot_internal_ids = [pdb_chain_to_internal_id[c_id] for c_id in prot_chain_ids]
        ab_internal_ids = [pdb_chain_to_internal_id[c_id] for c_id in ab_chain_ids]

         # Select atoms for protein and antibody based on chain IDs
        protein_selection = ' or '.join([f'chainid {chain_id}' for chain_id in prot_internal_ids])
        antibody_selection = ' or '.join([f'chainid {chain_id}' for chain_id in ab_internal_ids])

        protein_atom_indices = complex_pdb.topology.select(protein_selection)
        antibody_atom_indices = complex_pdb.topology.select(antibody_selection)
        #print(antibody_atom_indices)
        #print(protein_atom_indices)

        #  new structures for protein and antibody separately
        protein_pdb = complex_pdb.atom_slice(protein_atom_indices)
        antibody_pdb = complex_pdb.atom_slice(antibody_atom_indices)
        # print(protein_pdb)
 
        print(f"Processing PDB file: {pdb_file.split('/')[-1]}")

        # Probe radii (in nm)
        Rps = [0.55]  # Probe radii (in nm)

        results = []
        for Rp in Rps:
            # Compute SASA for the entire complex, protein, and antibody
            sasa_complex_total = md.shrake_rupley(complex_pdb, probe_radius=Rp).sum(axis=1).sum()
            sasa_protein_total = md.shrake_rupley(protein_pdb, probe_radius=Rp).sum(axis=1).sum()
            sasa_antibody_total = md.shrake_rupley(antibody_pdb, probe_radius=Rp).sum(axis=1).sum()

            # Compute SASA for specific residue groups
            sasa_complex_charged = calculate_sasa_for_residues(complex_pdb, Rp, Cr)
            sasa_protein_charged = calculate_sasa_for_residues(protein_pdb, Rp, Cr)
            sasa_antibody_charged = calculate_sasa_for_residues(antibody_pdb, Rp, Cr)

            sasa_complex_apolar = calculate_sasa_for_residues(complex_pdb, Rp, ApR)
            sasa_protein_apolar = calculate_sasa_for_residues(protein_pdb, Rp, ApR)
            sasa_antibody_apolar = calculate_sasa_for_residues(antibody_pdb, Rp, ApR)

            sasa_complex_polar = calculate_sasa_for_residues(complex_pdb, Rp, pR)
            sasa_protein_polar = calculate_sasa_for_residues(protein_pdb, Rp, pR)
            sasa_antibody_polar = calculate_sasa_for_residues(antibody_pdb, Rp, pR)

            # Compute Buried Surface Area (BSA)
            BSA = ((sasa_protein_total + sasa_antibody_total) - sasa_complex_total) * 100
            BSA_CC = ((sasa_protein_charged + sasa_antibody_charged) - sasa_complex_charged) * 100
            BSA_AA = ((sasa_protein_apolar + sasa_antibody_apolar) - sasa_complex_apolar) * 100
            BSA_PP = ((sasa_protein_polar + sasa_antibody_polar) - sasa_complex_polar) * 100

            # Compute non-interacting surface areas
            non_interacting_charged = sasa_complex_charged - (sasa_protein_charged + sasa_antibody_charged - sasa_complex_charged)
            non_interacting_apolar = sasa_complex_apolar - (sasa_protein_apolar + sasa_antibody_apolar - sasa_complex_apolar)
            non_interacting_polar = sasa_complex_polar - (sasa_protein_polar + sasa_antibody_polar - sasa_complex_polar)

            # Compute fractions
            NIS_frac_PP = non_interacting_polar / sasa_complex_polar if sasa_complex_polar != 0 else 0
            NIS_frac_AA = non_interacting_apolar / sasa_complex_apolar if sasa_complex_apolar != 0 else 0
            NIS_frac_CC = non_interacting_charged / sasa_complex_charged if sasa_complex_charged != 0 else 0

            results.append((Rp, BSA, BSA_CC, BSA_AA, BSA_PP, NIS_frac_PP, NIS_frac_AA, NIS_frac_CC))

            return pdb_file, results
# Use multiprocessing to parallelize the computation
if  __name__ == "__main__":
    output_file = "name_of_output_sasa_file.txt"
    with open(output_file, "w") as f:
        f.write("Probe_Radius\tBSA\tBSA_CC\tBSA_AA\tBSA_PP\tNIS_frac_PP\tNIS_frac_AA\tNIS_frac_CC\n")

    with Pool(cpu_count()) as pool:
        results = pool.map(process_pdb_file, pdb_files)

    with open(output_file, "a") as f:
        for result in results:
            if result is not None:  # Ensure the result is not None
                pdb_file, pdb_results = result
                f.write(f"# Input PDB file: {pdb_file.split('/')[-1]}\n")
                for res in pdb_results:
                    f.write(f"{res[0]:.2f}\t{res[1]:.2f}\t{res[2]:.2f}\t{res[3]:.2f}\t{res[4]:.2f}\t{res[5]:.2f}\t{res[6]:.2f}\t{res[7]:.2f}\n")
