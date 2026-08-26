import csv
import sys
import os
import numpy as np
import json
import argparse
from prody import *
from sidechainnet.utils.measure import *
from tqdm import tqdm

def tocdr(resseq):
    if 27 <= resseq <= 38:
        return '1'
    elif 56 <= resseq <= 65:
        return '2'
    elif 105 <= resseq <= 117:
        return '3'
    else:
        return '0'
    
def main(args):
    hchain = parsePDB(args.path + "/" + args.pdb, model=1, chain=args.hchain)
    _, _, hseq, _, _ = get_seq_coords_and_angles(hchain)
    N_atoms = hchain.select('backbone').select('name N')
    CA_atoms = hchain.select('backbone').select('name CA')
    C_atoms = hchain.select('backbone').select('name C')
    O_atoms = hchain.select('backbone').select('name O')
    N_coo = N_atoms.getCoords()
    CA_coo = CA_atoms.getCoords()
    C_coo = C_atoms.getCoords()
    O_coo = O_atoms.getCoords()
    hcdr = ''.join([tocdr(res.getResnum()) for res in hchain.iterResidues()])
    hcdr = hcdr[:len(hseq)]
    N_coo = eval(np.array2string(N_coo, separator=',', threshold=np.inf, precision=3, suppress_small=True))
    CA_coo = eval(np.array2string(CA_coo, separator=',', threshold=np.inf, precision=3, suppress_small=True))
    C_coo = eval(np.array2string(C_coo, separator=',', threshold=np.inf, precision=3, suppress_small=True))
    O_coo = eval(np.array2string(O_coo, separator=',', threshold=np.inf, precision=3, suppress_small=True))
    hcoords = {"N":N_coo,"CA":CA_coo,"C":C_coo,"O":O_coo}
    data = {"pdb": args.pdb, 
            "seq": hseq, "cdr": hcdr, 
            "coords": hcoords}
    with open(args.out, 'w') as f:
        json.dump(data, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, default='data')
    parser.add_argument("--pdb", type=str, default='test.pdb')
    parser.add_argument("--hchain", type=str, default='H')
    parser.add_argument("--out", type=str, default='out')
    
    args = parser.parse_args()
    main(args)