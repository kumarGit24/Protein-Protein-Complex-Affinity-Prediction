# Protein-Protein Complex Affinity Prediction #
Machine learning framework for predicting protein-protein complex binding affinities.

# Decriptor extraction Python Codes

# Requirements 
Python 3.10+

# Downloads and Installation Required

1. PDBFixer, PDBParser (https://github.com/openmm/pdbfixer)
2. Numpy, Pandas (https://numpy.org/, https://pandas.pydata.org/
3. MDAnalysis (https://www.mdanalysis.org/) 
4. Foldx (https://foldxsuite.crg.eu/academic-license-info; foldx-20270131)
5. Rosseta (https://rosettacommons.org/software/download)
6. HBPlus (https://github.com/RomanLas/HBPLUS)


# Step 1

All the PDBs requires to go thorugh the PDBFixer script as it cleans out any clashes among the AA.

# Step 2

Subsequently, all the cleaned 3D structure files (i.e. PDBs) (from Step 1) were processed to the chimerax script to keep only those chains from protein-protein complex that are required for the binding affinity (Kd).

# Step 3

Once, the cleaned structures were obtained from previous steps, the PDBs were processed to extract the structural and energetic desciptors (or features).

# Step 4

The obtained descriptors were arranged for each PDBID in a TSV file (e.g. test file in data folder), which was used to train the machine learning (ML) model.

# Step 5

Finally, a trained model was saved to predict the binding energy (or affinity) of the bNAb/HIV-1 Env trimer complex.
