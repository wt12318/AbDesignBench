import tomllib
import argparse
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import shutil
import re
from Bio.PDB import PDBParser
from collections import defaultdict
###set project name
print("Start Generating...")
pro_name = str(int(time.time() * 1000))
print(f"The project name is {pro_name}")

def find_start_line(file_path, start_marker):
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if line.strip().startswith(start_marker):
                return i
    return 0 

def get_tfold_metric(pdb_path):
    line_numbers = [2,3,4]
    lDDT = []
    pTM = []
    ids = []
    for item in os.listdir(pdb_path):
        full_path = os.path.join(pdb_path, item)
        with open(full_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
            lines = [lines[i-1].strip() for i in line_numbers if 0 < i <= len(lines)]
            metrics = [float(re.sub(r'.+: ', r'', i)) for i in lines]
            lDDT.append(metrics[0])
            pTM.append(metrics[2])
            ids.append(re.sub(r'.pdb', r'', item))
    return pd.DataFrame({"combined_id":ids,"lDDT":lDDT,"pTM":pTM})

def get_igfold_metric(fasta_path,pro_name):
    commd = f"/home/ab/run/AbRSA -i {fasta_path} -o /home/ab/run/tmp/{pro_name}_out > /home/ab/run/tmp/{pro_name}_out.txt"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    sl = find_start_line(f"/home/ab/run/tmp/{pro_name}_out.txt","#similarity")
    df = pd.read_csv(f"/home/ab/run/tmp/{pro_name}_out.txt", skiprows=sl+1, header=None,
                     sep=": ", engine="python", names=["type","seq"])
    pdb_path = re.sub(r'.fasta', r'', i) + ".pdb"
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_path)
    residue_b_factors = defaultdict(list)
    for model in structure:
        for chain in model:
            if chain.id == "H":
                for residue in chain:
                    for atom in residue:
                        residue_id = f"{residue.get_resname()}_{residue.get_id()[1]}"
                        residue_b_factors[residue_id].append(atom.get_bfactor())
    all_rmsd = []
    for residue, bf in residue_b_factors.items():
        all_rmsd.extend(list(set(bf)))
    HCDR3 = np.cumsum([len(i) for i in df.seq[0:7]])[-3:-1].tolist() ###HCDR3
    os.remove(f"/home/ab/run/tmp/{pro_name}_out.txt")
    return np.round(np.mean(all_rmsd[HCDR3[0]:HCDR3[1]]),3)

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", help="config file")
args = parser.parse_args()
config_file = args.config
##read config file
with open(config_file, "rb") as f:
    config = tomllib.load(f)
###paras
print("Parsing configuration files ...")
fasta_file = config["input"]["fasta_file"]
ag_pdb_file = config["input"]["ag_pdb_file"]
design_method = config["settings"]["design_method"]
design_num = config["settings"]["design_num"]
pre_structure = config["settings"]["pre_structure"]
docking = config["settings"]["docking"]
out_path = config["settings"]["out_path"]
####
commd = f"/home/ab/run/AbRSA -i {fasta_file} -o /home/ab/run/tmp/{pro_name}_fasta_out > /home/ab/run/tmp/{pro_name}_fasta_out.txt"
result = subprocess.run(commd, shell=True, text=True, capture_output=True)
sl = find_start_line(f"/home/ab/run/tmp/{pro_name}_fasta_out.txt","#similarity")
df = pd.read_csv(f"/home/ab/run/tmp/{pro_name}_fasta_out.txt", skiprows=sl+1, header=None,
                 sep=": ", engine="python", names=["type","seq"])
os.remove(f"/home/ab/run/tmp/{pro_name}_fasta_out.txt")
###
records = SeqIO.parse(fasta_file, "fasta")
data = [(record.id, str(record.seq)) for record in records]
ag_seq = data[2][1] ##antigen sequence
h_seq = data[0][1] ##h chain sequence
l_seq = data[1][1] ##l chain sequence

se = np.cumsum([len(i) for i in df.seq[0:7]])[-3:-1].tolist()
res_all = pd.DataFrame({"ID":[],"Sequence":[],"Method":[]})
if "IgLM" in design_method:
    print("Run IgLM ...")
    commd = f"conda run -n IgLM iglm_infill {fasta_file} H {se[0]} {se[1]} --chain_token [HEAVY] --species_token [HUMAN] --num_seqs {str(design_num)} --output_dir {out_path}/iglm_gen_{pro_name}"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ###read res
    records = SeqIO.parse(f"{out_path}/iglm_gen_{pro_name}/generated_seqs.fasta", "fasta")
    data = [(record.id, str(record.seq)) for record in records]
    res_iglm = pd.DataFrame(data, columns=["ID", "Sequence"])
    res_iglm["Method"] = "IgLM"
    ###del
    commd = f"rm -rf {out_path}/iglm_gen_{pro_name}"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    res_all = pd.concat([res_all,res_iglm], axis=0, ignore_index=True)
    print("Done !")

if "AbGPT" in design_method:
    print("Run AbGPT ...")
    start_res = df.seq[0]+df.seq[1]+df.seq[2]+df.seq[3]+df.seq[4]
    commd = f"cd /home/ab/abgpt; conda run -n IgLM abgpt_generate --chain_type heavy --starting_residue {start_res} --num_seqs {design_num} --out {out_path}/abgpt_gen_{pro_name}"
    res = subprocess.run(commd, shell=True, text=True, capture_output=True)
    ###read res
    tmp_lines = []
    lines = []
    with open(f'{out_path}/abgpt_gen_{pro_name}/heavy_{start_res}.txt', 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                tmp_lines.append(line)
            else:
                lines.append("".join(tmp_lines))
                tmp_lines = []
    res_abgpt = pd.DataFrame({"ID":["seq_"+str(i) for i in list(range(len(lines)))],"Sequence":lines})
    res_abgpt["Method"] = "AbGPT"
    ###del tmp
    commd = f"rm -rf {out_path}/abgpt_gen_{pro_name}/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    res_all = pd.concat([res_all,res_abgpt], axis=0, ignore_index=True)
    print("Done !")

if "ReprogBERT" in design_method:
    print("Run ReprogBERT ...")
    mask_input = df.seq[0]+df.seq[1]+df.seq[2]+df.seq[3]+df.seq[4] + "*"*len(df.seq[5]) + df.seq[6]
    commd = f"cd /home/ab/ReprogBERT; conda run -n ReprogBERT python main.py  --run inference --single_input {mask_input} --model_type base --exp_dir gen_{pro_name} --checkpoint ./output/base_cd3/checkpoints/last.ckpt --progen_dir ./progen --num_samples {design_num}"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ###read res
    records = SeqIO.parse(f"/home/ab/ReprogBERT/output/gen_{pro_name}/inference_smpl.fasta", "fasta")
    data = [(record.id, str(record.seq)) for record in records]
    res_reporgbert = pd.DataFrame(data, columns=["ID", "Sequence"])
    res_reporgbert["ID"] = res_reporgbert['ID'].str.replace('sample', 'seq')
    res_reporgbert["Method"] = "ReprogBERT"
    ##del gen
    commd = f"rm -rf /home/ab/ReprogBERT/output/gen_{pro_name}"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    res_all = pd.concat([res_all,res_reporgbert], axis=0, ignore_index=True)
    print("Done !")

if "PALM" in design_method:
    print("Run PALM ...")
    commd = f"cd /home/ab/PALM/Code/; conda run -n PALM bash run_test.sh {ag_seq} {h_seq} {se[0]} {se[1]} {len(df.seq[5])} {design_num} {out_path}/palm_gen_{pro_name}"
    print(commd)
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ##read res
    res = pd.read_csv(f"{out_path}/palm_gen_{pro_name}/result.csv")
    res["ID"] = ["seq_"+str(i) for i in range(len(res))]
    res_palm = res[["ID","Heavy_Chain"]].rename(columns={"Heavy_Chain":"Sequence"})
    res_palm["Method"] = "PALM"
    ###del tmp
    commd = f"rm -rf {out_path}/palm_gen_{pro_name}/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    res_all = pd.concat([res_all,res_palm], axis=0, ignore_index=True)
    print("Done !")

res_all["combined_id"] = [res_all.Method[i]+"_"+res_all.ID[i] for i in range(len(res_all))]
if pre_structure:

    ##tfold
    print("Predicting structure using tFold ...")
    os.makedirs(f"{out_path}/tfold_{pro_name}/ids/", exist_ok=True)
    os.makedirs(f"{out_path}/tfold_{pro_name}/fasta/", exist_ok=True)
    os.makedirs(f"{out_path}/tfold_{pro_name}/pdbs/", exist_ok=True)
    for _, row in res_all.iterrows():
        record = []
        record.append(SeqRecord(
            Seq(row['Sequence']),
            id="H",
            description=""
        ))
        record.append(SeqRecord(
            Seq(l_seq),
            id="L",
            description=""
        ))
        filename = f"{out_path}/tfold_{pro_name}/fasta/{row['combined_id']}.fasta"
        SeqIO.write(record, filename, "fasta")
    with open(f"{out_path}/tfold_{pro_name}/ids/seqid.txt", 'w') as f:
        for item in res_all["combined_id"]:
            f.write("%s\n" % item)

    commd = f"cd /home/ab/tFold/; conda run -n tfold python projects/tfold_ab/predict.py --pid_fpath={out_path}/tfold_{pro_name}/ids/seqid.txt --fas_dpath={out_path}/tfold_{pro_name}/fasta/ --pdb_dpath={out_path}/tfold_{pro_name}/pdbs/ --mdl_dpath params/params"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ###add scores
    scores = get_tfold_metric(f'{out_path}/tfold_{pro_name}/pdbs/')
    res_all = pd.merge(res_all, scores, on='combined_id')

    ###zip files
    shutil.make_archive(f"{out_path}/tfold_pdb_{pro_name}", "zip", f"{out_path}/tfold_{pro_name}/pdbs/")
    commd = f"rm -rf {out_path}/tfold_{pro_name}/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    print("Done !")

    ####lgfold
    print("Predicting structure using igFold ...")
    os.makedirs(f"{out_path}/lgfold_{pro_name}/", exist_ok=True)
    for _, row in res_all.iterrows():
        commd = f"conda run -n igfold python /home/ab/software/igfold_run.py -p {out_path}/lgfold_{pro_name}/{row['combined_id']}.pdb -a {row['Sequence']} -l {l_seq}"
        subprocess.run(commd, shell=True, text=True, capture_output=False)
    ###read metric
    fasta_files = [f"{out_path}/lgfold_{pro_name}/"+res_all.combined_id[i]+".fasta" for i in range(len(res_all))]
    all_hcdr3_rmsd = []
    for i in fasta_files:
        all_hcdr3_rmsd.append(get_igfold_metric(i, pro_name))
    res_all["HCDR3_RMSD"] = all_hcdr3_rmsd
    ##zip
    subprocess.run(f"rm {out_path}/lgfold_{pro_name}/*.fasta", shell=True, text=True, capture_output=True)
    shutil.make_archive(f"{out_path}/lgfold_pdb_{pro_name}", "zip", f"{out_path}/lgfold_{pro_name}/")
    commd = f"rm -rf {out_path}/lgfold_{pro_name}/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    print("Done !")

if docking:
    print("Performing Antibody-Antigen Docking ...")
    ####docking
    record = []
    for _, row in res_all.iterrows():
        record.append(SeqRecord(
            Seq(row['Sequence']+":"+l_seq),
            id=row['combined_id'],
            description=""
        ))
    os.makedirs(f"{out_path}/docking_{pro_name}/", exist_ok=True)
    SeqIO.write(record, f"{out_path}/docking_{pro_name}/fd_input.fasta", "fasta")
    commd = f"cd /home/ab/Fold-Dock; conda run -n fd python fold_dock.py {out_path}/docking_{pro_name}/fd_input.fasta -a {ag_pdb_file} -o {out_path}/docking_{pro_name}/ -t 1"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ##read scores
    all_docking_score = []
    for _, row in res_all.iterrows():
        score_file = f"{out_path}/docking_{pro_name}/{row['combined_id']}/scores.csv"
        tmp_score = pd.read_csv(score_file)
        all_docking_score.append(round(tmp_score.score[0],3))
    res_all["Docking_score"] = all_docking_score
    ##zip files
    shutil.make_archive(f"{out_path}/docking_res", "zip", f"{out_path}/docking_{pro_name}/")
    commd = f"rm -rf {out_path}/docking_{pro_name}/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    print("Done !")

###output
res_all.to_csv(f"{out_path}/seq_model_res.csv")
print("All Done !")
###clean tmp file
commd = f"rm -rf /home/ab/run/tmp/{pro_name}*"
subprocess.run(commd, shell=True, text=True, capture_output=False)

