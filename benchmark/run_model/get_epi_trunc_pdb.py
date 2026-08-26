import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
from abnumber import Chain
from Bio.SeqUtils import seq1
from abnumber import Chain
from Bio import PDB
import copy
import argparse

def get_cb_or_ca(residue):
    """Return the Cβ atom if present, otherwise return the Cα atom (for Glycine)."""
    if "CB" in residue:
        return residue["CB"]
    elif "CA" in residue:  # Handle Glycine (no Cβ)
        return residue["CA"]
    return None

def calculate_shortest_beta_carbon_distances(pdb_file, chain1_id, chain2_id, res_ids):
    """
    Calculate the shortest Cβ distance between all residues in chain1 and 
    the specified residues in chain2.
    
    Parameters:
    - pdb_file: str, path to the PDB file
    - chain1_id: str, chain containing multiple residues to compare
    - chain2_id: str, chain containing specified residues
    - res_ids: list of int, residue numbers in chain2 to compare distances
    
    Returns:
    - List of tuples (residue in chain1, closest residue in chain2, shortest distance)
    """

    # Parse the PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    # Extract chains
    model = structure[0]  # Use the first model
    chain1 = model[chain1_id]
    chain2 = model[chain2_id]

    # Get Cβ atoms of specified residues in chain2
    target_atoms = []
    for residue in chain2:
        if residue.get_id()[1] in res_ids:
            atom = get_cb_or_ca(residue)
            if not(atom  is None) and atom.is_disordered() == 0:
                target_atoms.append((residue.get_id()[1], atom))  # Store (residue number, atom)
    
    if not target_atoms:
        raise ValueError(f"No valid Cβ or Cα atoms found in chain {chain2_id} for residues {res_ids}")

    # Compute shortest distances
    shortest_distances = []
    
    for residue in chain1:
        cb_atom = get_cb_or_ca(residue)
        if not(cb_atom  is None) and cb_atom.is_disordered() == 0:
            min_distance = float("inf")
            closest_residue = None
            
            # Compare against all specified residues in chain2
            for res_id, target_atom in target_atoms:
                distance = cb_atom - target_atom  # BioPython allows atom subtraction for distance
                if distance < min_distance:
                    min_distance = distance
                    closest_residue = res_id
            
            shortest_distances.append((residue.get_id()[1], closest_residue, min_distance))

    return shortest_distances

def find_consecutive_letters(chars, subscripts):
    result = []
    for index in subscripts:
        target = chars[list(chars.keys())[index]]  # Get the target character
        left, right = index, index
        # Expand to the left
        while left > 0 and chars[list(chars.keys())[left - 1]] == target:
            left -= 1
        # Expand to the right
        while right < len(chars) - 1 and chars[list(chars.keys())[right + 1]] == target:
            right += 1
        # Append the consecutive sequence
        result.extend(list(chars.keys())[left:right + 1])
    return set(result)
    
def extract_secondary_structure(input_pdb, output_pdb, chain_residues, need_chain):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    model = structure[0]  # Assuming single model
    io = PDB.PDBIO()
    io.set_structure(structure)
    
    # DSSP to calculate secondary structure
    dssp = PDB.DSSP(model, input_pdb, dssp='mkdssp')
    secondary_structures = {}
    
    # Map residues to their secondary structure
    for key in dssp.keys():
        (chain, res_id) = key[0], key[1][1]
        if chain == need_chain:
            secondary_structures[(chain, res_id)] = dssp[key][2]   # Secondary structure type
    ##idx
    idx = 0
    idx_list = []
    for i in secondary_structures:
        if i in res:
            idx_list.append(idx)
        idx += 1

    con_aa = find_consecutive_letters(secondary_structures, idx_list)
    filter_secondary_structures = copy.deepcopy(secondary_structures)
    for i in secondary_structures:
        if i not in con_aa:
            del filter_secondary_structures[i]
            
    # Select residues with matching secondary structures
    class SecStructSelect(PDB.Select):
        def accept_residue(self, residue):
            res_chain = residue.get_parent().id
            res_id = residue.id[1]
            return (res_chain, res_id) in filter_secondary_structures 
    
    io.save(output_pdb, SecStructSelect())
    print(f"Secondary structure segment saved to {output_pdb}")

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdb", help="input pdb file")
parser.add_argument("-a", "--hchain", help="h chain ID")
parser.add_argument("-b", "--agchain", help="antigen chain ID")
parser.add_argument("-o", "--output", help="output pdb file")
parser.add_argument("-e", "--epitope", help="output epitope information file")
args = parser.parse_args()
pdb_file = args.pdb
hchain = args.hchain
ag = args.agchain
out_pdb= args.output
out_info = args.epitope

p = PDBParser(QUIET=True)
struc = p.get_structure("protein", pdb_file)
##HCDR3残基的index
chains_seq = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in struc.get_chains()}
chain_resid = {chain.id:[residue.id[1] for residue in chain] for chain in struc.get_chains()}
chain_anno = Chain.multiple_domains(chains_seq[hchain], scheme='chothia')
for i in chain_anno:
    if i.chain_type == "H":
        hseq = i
hcdr3_start = chains_seq[hchain].find(hseq.cdr3_seq)
hcdr3_resid = chain_resid[hchain][hcdr3_start:hcdr3_start+len(hseq.cdr3_seq)]
###表位氨基酸残基
results = calculate_shortest_beta_carbon_distances(pdb_file, ag, hchain, hcdr3_resid)
epi = []
for i in results:
    if i[2] < 8: ##与HCDR3的距离小于8A
        epi.extend([i[0]])
if len(epi) == 0:
   print(f"{pdb_file}: No epitopes!") 
else:
    target_residues = [(ag, i) for i in epi]
    ###表位氨基酸10A附近的残基
    u = mda.Universe(pdb_file)
    cutoff_distance = 10.0  # Ångstrom
    # Get all atoms within cutoff distance
    selected_atoms = set()
    for chain_id, res_id in target_residues:
        sel_res = u.select_atoms(f"protein and segid {chain_id} and resid {res_id}")
        if len(sel_res) == 0:
            continue
        distances = distance_array(sel_res.positions, u.atoms.positions)
        close_atoms = u.atoms[np.any(distances < cutoff_distance, axis=0)]
        selected_atoms.update(close_atoms)
    # Convert to residue set
    selected_residues = {atom.residue for atom in selected_atoms}
    res = [(i.segid,i.resnum) for i in selected_residues if i.segid == ag]
    
    ####保留二级结构
    extract_secondary_structure(pdb_file, out_pdb, res, ag)

    ##保存表位信息和HCDR3长度信息
    epi_info = [f"{i[0]}{i[1]}" for i in target_residues]
    epi_info.extend([str(len(hseq.cdr3_seq))])
    with open(out_info, 'w') as f:
        f.write(f"{epi_info}\n")

