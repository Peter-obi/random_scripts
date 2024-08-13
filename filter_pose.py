import argparse
from Bio.PDB import PDBParser, NeighborSearch
from rdkit import Chem
import numpy as np
import os

def ligand_pose_filter(protein_structure_file, ligand_pose_file, key_residues, cutoff_distances):
    """
    This function filters ligand poses based on their proximity to specified residues using RDKit for handling the ligands,
    which can be in SDF, PDB, MOL or MOL2 formats (most popular formats), and BioPython for the proteins.

    Args:
    protein_structure_file (str): Path to the protein structure file in PDB format.
    ligand_pose_file (str): Path to the ligand pose file, supports formats like SDF, PDB, MOL or MOL2.
    key_residues (list): List of key residues in the format (chain_id, residue_id)  e.g., A1000 B3000.
    cutoff_distances (list): List of cutoff distances in angstroms, corresponding to each key residue.

    Returns:
    list of list of bool: Indicates for each ligand whether it is within the cutoff distance of the corresponding residues.
    """
    #load the protein structure from the specified file using biopython
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('Protein', protein_structure_file)
    
    #determine the file extension of the ligand file to decide the appropriate RDKit reader
    _, file_ext = os.path.splitext(ligand_pose_file)
    if file_ext.lower() in ['.sdf', '.mol', '.mol2']:  #check if SDF and Mol2
        ligands = Chem.SDMolSupplier(ligand_pose_file, removeHs=False)
    elif file_ext.lower() == '.pdb':  #check if pdb
        ligands = [Chem.MolFromPDBFile(ligand_pose_file, removeHs=False)]
    else:
        raise ValueError("Unsupported file format for ligands.")

    #initalize an empty list for the results
    results = []

    #iterate over the ligand
    for ligand in ligands:
        if ligand is None:
            # if RDKit failed to load a ligand, skip it.
            continue

        #get the 3D coordinates of each atom in the ligand. assumes since ligands are docked, conformers are ok
        conf = ligand.GetConformer()
        ligand_positions = [conf.GetAtomPosition(i) for i in range(ligand.GetNumAtoms())]

        #initialize the list to collect the proximity check results for the current ligand against all provided key residues
        proximity_results = []
        for (chain_id, residue_id), cutoff in zip(key_residues, cutoff_distances):
            try:
                #access the specified residue in the protein structure
                residue = structure[0][chain_id][residue_id]
                #get coordinates for all atoms in the provided key residue
                residue_positions = [atom.get_coord() for atom in residue.get_atoms()]

                #check if any atom in the ligand is within the cutoff distance of any atom in the residue
                is_within_cutoff = any(
                    np.linalg.norm(np.array(ligand_pos) - np.array(res_pos)) <= cutoff
                    for ligand_pos in ligand_positions
                    for res_pos in residue_positions
                )
                proximity_results.append(is_within_cutoff)
            except KeyError:
                #if the residue specified does not exist in the protein structure, log it as False.
                proximity_results.append(False)

        #append the proximity results of the current ligand to the results list.
        results.append(proximity_results)

    return results

if __name__ == "__main__":
    #set up command-line argument parsing.
    parser = argparse.ArgumentParser(description="Filter ligand poses based on proximity to key residues using RDKit")
    parser.add_argument("protein_structure_file", help="Path to the protein structure file in PDB format")
    parser.add_argument("ligand_pose_file", help="Path to the ligand pose file (e.g., SDF, PDB, or MOL2)")
    parser.add_argument("--key_residues", nargs='+', help="List of key residues, e.g., 'A100 B200'", required=True)
    parser.add_argument("--cutoff_distances", nargs='+', type=float, help="List of cutoff distances in angstroms", required=True)

    args = parser.parse_args()

    #convert the key residue inputs from string format to a tuple of (chain ID, residue number).
    key_residues = [(res[0], int(res[1:])) for res in args.key_residues]
    #execute the filtering function with the parsed arguments and print the results.
    result = ligand_pose_filter(args.protein_structure_file, args.ligand_pose_file, key_residues, args.cutoff_distances)
    print(result)

#example usage: python filter_ligand_pose.py 6kqi.pdb ligand.pdb --key_residues A241 A154 A200 --cutoff_distances 5.0 4.5 6.0
#install biopython - pip install biopython and rdkit - pip install rdkit before usage



