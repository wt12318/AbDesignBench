from igfold import IgFoldRunner
from igfold.refine.pyrosetta_ref import init_pyrosetta
import argparse
###参数
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdb", help="output pdb file")
parser.add_argument("-a", "--hchain", help="hchain")
parser.add_argument("-l", "--lchain", help="lchain")
args = parser.parse_args()
pdb = args.pdb
hchain = args.hchain
lchain = args.lchain

init_pyrosetta()

if lchain == "NA":
    sequences = {"H": hchain}
else:
    sequences = {"H": hchain, "L": lchain}

igfold = IgFoldRunner()
out = igfold.fold(
    pdb, # Output PDB file
    sequences=sequences, # Nanobody sequence
    do_refine=True, # Refine the antibody structure with PyRosetta
    do_renum=False, # Renumber predicted antibody structure (Chothia)
)
