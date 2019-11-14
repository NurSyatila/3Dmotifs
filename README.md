# 3Dmotifs
Sub-structural similarity searches that include searching for similarly arranged amino acid arrangements in the Protein Data Bank

This repository contain automated scripts to search and analyse substructures.
# AUTO_IMAAAGINE - Automated computational pipeline used to search for user-defined amino acid arrangement in the PDB
1. Search for user-defined hypothetical amino acid arrangement using IMAAAGINE algorithm (http://mfrlab.org/imaaagine/) with pre-defined inter-residue distance and distance tolerance 
2. Sequentially run multiple IMAAAGINE searches
3. Obtain output from the IMAAAGINE searches and generate PDB-formatted patterns for substructural comparison
4. Perform all-against-all superposition on substructures and cluster similarly arranged substructures (conserved substructure) based on RMSD values through hierarchical clustering
5. Perform fold comparison on proteins with conserved substructures for detection of homology using DaliLite
6. Perform multiple sequence alignment on proteins with conserved substructures for detection of homology using MUSCLE
7. Compile structural and biological informations of proteins with similarly arranged substructures useful for functional validation of conserved substructure

# AUTO_ASSAM - Automated computational pipeline used to search for substructures similar to known 3D sites in the PDB (tested on annotated drug binding sites)
1. Obtain protein-drug complexes from the PDB (https://www.rcsb.org/pdb/ligand/drugMapping.do)
2. Get sets of drug binding residues within distance of 4.0A and save PDB-formatted patterns
3. Sequentially run multiple ASSAM searches for sub-structural similarity searches in the PDB (http://mfrlab.org/assam/)
4. Obtain output from the ASSAM searches
5. Perform pairwise superposition on known substructure (drug binding site) and hit substructure (similarly arranged substructure) based on pairwise RMSD value
6. Perform fold comparison on proteins containing similarly arranged substructures for detection of homology using DaliLite
7. Perform multiple sequence alignment on proteins containing similarly arranged substructures for detection of homology using MUSCLE
8. Compile structural and biological informations of proteins with similarly arranged substructures useful for pairwise comparison and possible functional annotation transfer
