import pandas as pd
import subprocess
from tqdm import tqdm
import numpy as np

import argparse
###参数
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input file")
parser.add_argument("-p", "--path", help="input path")
parser.add_argument("-o", "--out", help="output file")
args = parser.parse_args()
in_file = args.input
in_path = args.path
out_path = args.out

dt = pd.read_csv(in_file,sep=" ",names=["gen_file","ori_file","pdb_id"])
res = []
for i in tqdm(range(len(dt))):
    f1 = in_path + "/extract_HCDR3/" + dt.ori_file[i]
    f2 = in_path + "/gen_extract_HCDR3/" + dt.gen_file[i]
    out = subprocess.run("/home/data/sdb/wt/miniconda3/envs/ab_rmsd/bin/calculate_rmsd"+
                 " --only-alpha-carbons " + f1 + " " + f2, shell=True,capture_output=True, text=True) 
    if "error" in out.stderr or "error" in out.stdout:
        tmp_res = np.nan
    else:
        tmp_res = float(out.stderr.replace("\n",""))
    res.extend([tmp_res])

dt["HCDR3_RMSD"] = res
dt.to_csv(out_path)
