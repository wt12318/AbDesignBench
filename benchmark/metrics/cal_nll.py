from antiberty import AntiBERTyRunner
import argparse
###参数
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input csv")
args = parser.parse_args()
input_file = args.input

import pandas as pd
dt = pd.read_csv(input_file)
antiberty = AntiBERTyRunner()
from tqdm import tqdm
score = []
batch_size_ = 100
for i in tqdm(range(0, len(dt), batch_size_)):
    batch_end = min(i + batch_size_, len(dt))
    seq = dt.gen_seq[i:batch_end].to_list()
    pll = antiberty.pseudo_log_likelihood(seq, batch_size=64)
    score.extend(pll.tolist())
dt["nll"] = score 
dt.to_csv(input_file)
