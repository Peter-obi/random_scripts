#Parse a pdb file and seperate atoms such as nucleic acid from heteroatomsa nad save the heteroatoms to another chain
#You need to only change 'your_pdb' to your pdb of interest.

from Bio.PDB import PDBParser, PDBIO, Select, Chain

# Select atom criteria
class RNASelect(Select):
    def accept_residue(self, residue):
        return residue.get_resname() in ["A", "C", "G", "U"]

class HeteroSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] != " "

# Parse the PDB file
parser = PDBParser()
structure = parser.get_structure("your_Structure", "your.pdb")

# Separate RNA and heteroatoms
rna_chain_X = [res for res in structure[0]["X"] if RNASelect().accept_residue(res)]
rna_chain_Y = [res for res in structure[0]["Y"] if RNASelect().accept_residue(res)]
hetero_chain_X = [res for res in structure[0]["X"] if HeteroSelect().accept_residue(res)]
hetero_chain_Y = [res for res in structure[0]["Y"] if HeteroSelect().accept_residue(res)]

# Make new structure
new_structure = structure.copy()

# Delete original chains
del new_structure[0]["X"]
del new_structure[0]["Y"]

# Create and add new chains
new_chain_X = Chain.Chain("X")
new_chain_Y = Chain.Chain("Y")
new_chain_Z = Chain.Chain("Z")
new_structure[0].add(new_chain_X)
new_structure[0].add(new_chain_Y)
new_structure[0].add(new_chain_Z)

# Add RNA residues to new chains X and Y
for res in rna_chain_X:
    new_chain_X.add(res.copy())
for res in rna_chain_Y:
    new_chain_Y.add(res.copy())

# Add heteroatom residues to new chain Z
for res in hetero_chain_X + hetero_chain_Y:
    new_chain_Z.add(res.copy())

# Renumber atoms
atom_number = 1
for atom in new_structure.get_atoms():
    atom.serial_number = atom_number
    atom_number += 1

# Write out new structure to a PDB file
io = PDBIO()
io.set_structure(new_structure)
io.save("new_structure.pdb")
