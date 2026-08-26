import tomllib
import argparse
import subprocess
import pandas as pd
import numpy as np
import time
import os
import shutil
import re
from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1
import glob

def find_start_line(file_path, start_marker):
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if line.strip().startswith(start_marker):
                return i
    return 0

def remove_zero_coords(input_pdb, output_pdb):
    parser = PDBParser()
    structure = parser.get_structure('protein', input_pdb)

    atoms_to_remove = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if all(coord == 0 for coord in atom.get_coord()):
                        atoms_to_remove.append((chain.id, residue.id, atom.name))

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in list(residue):
                    if (chain.id, residue.id, atom.name) in atoms_to_remove:
                        residue.detach_child(atom.get_id())

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

###set project name
print("Start Generating...")
pro_name = str(int(time.time() * 1000))
print(f"The project name is {pro_name}")

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", help="config file")
args = parser.parse_args()
config_file = args.config
##read config file
with open(config_file, "rb") as f:
    config = tomllib.load(f)
###paras
print("Parsing configuration files ...")
pdb_file = config["input"]["pdb_file"]
design_method = config["settings"]["design_method"]
design_num = config["settings"]["design_num"]
hchain = config["settings"]["Heavy_Chain"]
agchain = config["settings"]["Antigen_Chain"]
out_path = config["settings"]["out_path"]
binding_affinity = config["settings"]["binding_affinity"]
binding_energy = config["settings"]["binding_energy"]

pdbparser = PDBParser(QUIET=True)
structure = pdbparser.get_structure("protein", f"{pdb_file}")
chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}
hseq = chains[f"{hchain}"]
agseq = chains[f"{agchain}"]

res_all = pd.DataFrame({"ID":[],"Sequence":[],"Method":[]})
if "DiffAb" in design_method:
    print("Run DiffAb ...")
    commd = f"cd /home/ab/diffab; conda run -n diffab bash run_test.sh {design_num} {hchain} {out_path}/diffab_{pro_name}/ {pdb_file}"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    gen_pdbs = glob.glob(f"{out_path}/diffab_{pro_name}/tmp/*/H_CDR3/")[0]
    os.remove(gen_pdbs + "/REF1.pdb")
    seqid = []
    sequence = []
    pdbparser = PDBParser()
    for item in os.listdir(gen_pdbs):
        full_path = os.path.join(gen_pdbs, item)
        structure = pdbparser.get_structure("protein", full_path)
        chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}
        sequence.append(chains[hchain])
        seqid.append(re.sub(r'.pdb', r'', item))
    res_diffab = pd.DataFrame({"ID":seqid,"Sequence":sequence,"Method":"DiffAb"})
    res_all = pd.concat([res_all,res_diffab], axis=0, ignore_index=True)
    ###zip pdbs
    shutil.move(gen_pdbs, f'{out_path}/diffab_{pro_name}/gen_pdb')
    if binding_affinity:
        print("Predicting affinity ...")
        all_pdbs = os.listdir(f"{out_path}/diffab_{pro_name}/gen_pdb")
        os.makedirs(f"{out_path}/diffab_{pro_name}/process_pdb/", exist_ok=True)
        affinity_res = pd.DataFrame({"id":[],"affinity":[]})
        for i in all_pdbs:
            full_path = f"{out_path}/diffab_{pro_name}/gen_pdb/" + i
            tmp_id = re.sub(r'.pdb',r'',i)
            commd = f"python sort_renumber.py -i {full_path} -o {out_path}/diffab_{pro_name}/process_pdb/{tmp_id}"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ##pre affinity
            commd = f"cd /home/ab/ANTIPASTI-main/notebooks; conda run -n antipasti-env python pre_affinity.py -i {out_path}/diffab_{pro_name}/process_pdb/{tmp_id}_chothia.pdb -n {tmp_id} -o {out_path}/diffab_{pro_name}/process_pdb/diffab_affinity.csv -d 0000"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ###read res
            tmp_aff = pd.read_csv(f"{out_path}/diffab_{pro_name}/process_pdb/diffab_affinity.csv")
            tmp_aff = tmp_aff.iloc[:, -2:]
            affinity_res = pd.concat([affinity_res, tmp_aff], axis=0, ignore_index=True)
        #affinity_res.to_csv(f"{out_path}/diffab_{pro_name}/diffab_affinity.csv")

    if binding_energy:
        print("Predicting binding energy ...")
        all_pdbs = os.listdir(f"{out_path}/diffab_{pro_name}/gen_pdb")
        os.makedirs(f"{out_path}/diffab_{pro_name}/rosetta_pdb/", exist_ok=True)
        ids = []
        dg = []
        for i in all_pdbs:
            full_path = f"{out_path}/diffab_{pro_name}/gen_pdb/" + i
            tmp_id = re.sub(r'.pdb',r'',i)
            commd = f"cd /home/ab/rosetta/run; ./score_jd2.static.linuxgccrelease -s {full_path} -no_optH false -ignore_unrecognized_res -out:pdb -out:path:all {out_path}/diffab_{pro_name}/rosetta_pdb/"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            commd = f"cd /home/ab/rosetta/run; ./InterfaceAnalyzer.static.linuxgccrelease -s {out_path}/diffab_{pro_name}/rosetta_pdb/{tmp_id}_0001.pdb @./pack_input_options.txt -interface {hchain}_{agchain} -out:path:all {out_path}/diffab_{pro_name}/rosetta_pdb/"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ###read score
            commd = "awk '{print $6}' " + f" {out_path}/diffab_{pro_name}/rosetta_pdb/pack_input_score.sc | sed '1d' > {out_path}/diffab_{pro_name}/rosetta_pdb/outscore"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            tmp_score = pd.read_csv(f"{out_path}/diffab_{pro_name}/rosetta_pdb/outscore")
            ids.append(tmp_id)
            dg.append(tmp_score.dG_separated.tolist()[0])
            ##clean
            commd = f"rm {out_path}/diffab_{pro_name}/rosetta_pdb/*"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
        energy_res = pd.DataFrame({"id":ids,"dg":dg})
        energy_res.to_csv(f"{out_path}/diffab_{pro_name}/gen_pdb/diffab_rosetta.csv")
    if binding_affinity:
        affinity_res.to_csv(f"{out_path}/diffab_{pro_name}/gen_pdb/diffab_affinity.csv")
    ###zip
    shutil.make_archive(f"{out_path}/diffab_gen_{pro_name}", "zip", f"{out_path}/diffab_{pro_name}/gen_pdb")
    commd = f"rm -rf {out_path}/diffab_{pro_name}"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    print("Done !")

if "AbOpt" in design_method:
    print("Run AbOpt ...")
    os.makedirs(f"{out_path}/abopt_{pro_name}/", exist_ok=True)
    commd1 = f"cd /home/ab/ab_opt/AbDock/; conda run -n antibody_dock python dock_pdb.py --heavy {hchain} --pdb_path {pdb_file} --config configs/test/dock_cdr.yml --ckpt reproduction/dock_single_cdr/250000.pt -n {design_num} -b 10 -d cuda -o {out_path}/abopt_{pro_name}/dock_single"
    subprocess.run(commd1, shell=True, text=True, capture_output=False)
    s1_path = glob.glob(f"{out_path}/abopt_{pro_name}/dock_single/dock_cdr/*/H_CDR3")[0]
    commd2 = f"cd /home/ab/ab_opt/AbDock/; conda run -n antibody_dock python optimize_ab.py --num_gpus 1 --process_per_gpu 1 --docked_pose_dir {s1_path} --seq_design_dir {out_path}/abopt_{pro_name}/seq_design --design_contig '' --screen_dir {out_path}/abopt_{pro_name}/screen --heavy_chain_id {hchain} --nums 1"
    subprocess.run(commd2, shell=True, text=True, capture_output=False)
    ##read res
    all_files = glob.glob(f"{out_path}/abopt_{pro_name}/seq_design/seq_design/*/aa.csv", recursive=True)
    ids = []
    sequences = []
    origin_seq = []
    for i in all_files:
        tmp = pd.read_csv(i)
        tmp_seq = tmp.sampled_aa[0]
        tmp_id = re.sub(r'.+/',r'',re.sub(r'_rosetta.pdb_/aa.csv', r'',i))
        if tmp_id != "REF1":
            ids.append(tmp_id)
            sequences.append(tmp_seq)
            origin_seq.append(tmp.native_aa[0])
    full_seqs = [re.sub(rf"{list(set(origin_seq))[0]}",rf"{i}",hseq) for i in sequences]
    full_seqs = [re.sub(r"XXXXXXXXXXXXX",r"",i) for i in full_seqs]
    res_abopt = pd.DataFrame({"ID":ids,"Sequence":full_seqs,"Method":"AbOpt"})
    res_all = pd.concat([res_all,res_abopt], axis=0, ignore_index=True)
    all_pdbs = glob.glob(f"{out_path}/abopt_{pro_name}/seq_design/seq_design/*/*/0000.pdb", recursive=True)
    os.makedirs(f"{out_path}/abopt_{pro_name}/gen_pdb")
    for i in all_pdbs:
        tmp_id = re.sub(r'.+/',r'',re.sub(r'_rosetta.+', r'',i))
        if tmp_id != "REF1":
            shutil.move(i, f"{out_path}/abopt_{pro_name}/gen_pdb/gen_{tmp_id}.pdb")
    if binding_affinity:
        print("Predicting Affinity ...")
        all_pdbs = os.listdir(f"{out_path}/abopt_{pro_name}/gen_pdb")
        os.makedirs(f"{out_path}/abopt_{pro_name}/process_pdb/", exist_ok=True)
        affinity_res = pd.DataFrame({"id":[],"affinity":[]})
        for i in all_pdbs:
            full_path = f"{out_path}/abopt_{pro_name}/gen_pdb/" + i
            tmp_id = re.sub(r'.pdb',r'',i)
            commd = f"python fill_pdb.py -o {pdb_file} -g {full_path} -a {hchain} -b {hchain} -s chothia -p {out_path}/abopt_{pro_name}/process_pdb/{tmp_id}"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ##pre affinity
            commd = f"cd /home/ab/ANTIPASTI-main/notebooks; conda run -n antipasti-env python pre_affinity.py -i {out_path}/abopt_{pro_name}/process_pdb/{tmp_id}_chothia.pdb -n {tmp_id} -o {out_path}/abopt_{pro_name}/process_pdb/abopt_affinity.csv -d 0000"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ###read res
            tmp_aff = pd.read_csv(f"{out_path}/abopt_{pro_name}/process_pdb/abopt_affinity.csv")
            tmp_aff = tmp_aff.iloc[:, -2:]
            affinity_res = pd.concat([affinity_res, tmp_aff], axis=0, ignore_index=True)
        #affinity_res.to_csv(f"{out_path}/abopt_{pro_name}/gen_pdb/abopt_affinity.csv")

    if binding_energy:
        print("Predicting binding energy ...")
        all_pdbs = os.listdir(f"{out_path}/abopt_{pro_name}/gen_pdb")
        os.makedirs(f"{out_path}/abopt_{pro_name}/rosetta_pdb/", exist_ok=True)
        ids = []
        dg = []
        for i in all_pdbs:
            full_path = f"{out_path}/abopt_{pro_name}/gen_pdb/" + i
            tmp_id = re.sub(r'.pdb',r'',i)
            commd = f"cd /home/ab/rosetta/run; ./score_jd2.static.linuxgccrelease -s {full_path} -no_optH false -ignore_unrecognized_res -out:pdb -out:path:all {out_path}/abopt_{pro_name}/rosetta_pdb/"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            commd = f"cd /home/ab/rosetta/run; ./InterfaceAnalyzer.static.linuxgccrelease -s {out_path}/abopt_{pro_name}/rosetta_pdb/{tmp_id}_0001.pdb @./pack_input_options.txt -interface {hchain}_{agchain} -out:path:all {out_path}/abopt_{pro_name}/rosetta_pdb/"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ###read score
            commd = "awk '{print $6}' " + f"{out_path}/abopt_{pro_name}/rosetta_pdb/pack_input_score.sc | sed '1d' > {out_path}/abopt_{pro_name}/rosetta_pdb/outscore"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            tmp_score = pd.read_csv(f"{out_path}/abopt_{pro_name}/rosetta_pdb/outscore")
            ids.append(tmp_id)
            dg.append(tmp_score.dG_separated.tolist()[0])
            ##clean
            commd = f"rm {out_path}/abopt_{pro_name}/rosetta_pdb/*"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
        energy_res = pd.DataFrame({"id":ids,"dg":dg})
        energy_res.to_csv(f"{out_path}/abopt_{pro_name}/gen_pdb/abopt_rosetta.csv")
    if binding_affinity:
        affinity_res.to_csv(f"{out_path}/abopt_{pro_name}/gen_pdb/abopt_affinity.csv")
    ###zip
    shutil.make_archive(f"{out_path}/abopt_gen_{pro_name}", "zip", f"{out_path}/abopt_{pro_name}/gen_pdb/")
    commd = f"rm -rf {out_path}/abopt_{pro_name}"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    print("Done !")

if "AbDockgen" in design_method:
    print("Run Abdockgen ...")
    os.makedirs(f"{out_path}/abdockgen_{pro_name}/", exist_ok=True)
    ###igmt number
    from renum import renumber
    pdb_id = re.sub(r'.+/',r'',re.sub(r'.pdb', r'',pdb_file))
    h,l,ag = renumber(f"{pdb_file}", f"{out_path}/abdockgen_{pro_name}/{pdb_id}_imgt.pdb", return_other_chains=True, scheme="imgt")
    commd = f"""
    cd /home/ab/abdockgen &&
    source activate abdockgen &&
    python process_data.py {out_path}/abdockgen_{pro_name}/{pdb_id}_imgt.pdb {hchain} {agchain} > {out_path}/abdockgen_{pro_name}/{pdb_id}.json &&
    conda deactivate
    """
    subprocess.run(commd, shell=True, executable="/bin/bash", check=True)
    ####
    commd = f"""
    cd /home/ab/abdockgen &&
    source activate abdockgen &&
    python gen_seq.py --json {out_path}/abdockgen_{pro_name}/{pdb_id}.json --batch_size 2 --num {design_num} --out {out_path}/abdockgen_{pro_name}/ &&
    conda deactivate
    """
    subprocess.run(commd, shell=True, executable="/bin/bash", check=True)
    ##read res
    abdockgen_res = pd.read_csv(f"{out_path}/abdockgen_{pro_name}/gen_seq.csv")
    gen_full_seq = [re.sub(rf"{abdockgen_res.origin_cdr[0]}",rf"{abdockgen_res.gen_cdr[i]}",hseq) for i in range(len(abdockgen_res))]
    gen_full_seq = [re.sub(r"XXXXXXXXXXXXX",r"",i) for i in gen_full_seq]

    abdockgen_res = pd.DataFrame({"ID":["Seq_"+str(i) for i in range(len(abdockgen_res))],"Sequence":gen_full_seq,"Method":"AbDockgen"})
    res_all = pd.concat([res_all, abdockgen_res], axis=0, ignore_index=True)
    os.makedirs(f"{out_path}/abdockgen_{pro_name}/gen_pdb")
    for i, row in abdockgen_res.iterrows():
        shutil.move(f"{out_path}/abdockgen_{pro_name}/{pdb_id}_imgt.pdb_pred_{i}.pdb",
                    f"{out_path}/abdockgen_{pro_name}/gen_pdb/{row['ID']}.pdb")
    if binding_affinity:
        print("Predicting affinity ...")
        all_pdbs = os.listdir(f"{out_path}/abdockgen_{pro_name}/gen_pdb")
        os.makedirs(f"{out_path}/abdockgen_{pro_name}/process_pdb/", exist_ok=True)
        affinity_res = pd.DataFrame({"id":[],"affinity":[]})
        for i in all_pdbs:
            full_path = f"{out_path}/abdockgen_{pro_name}/gen_pdb/" + i
            tmp_id = re.sub(r'.pdb',r'',i)
            commd = f"python fill_pdb.py -o {pdb_file} -g {full_path} -a {hchain} -b H -s imgt -p {out_path}/abdockgen_{pro_name}/process_pdb/{tmp_id}"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ##pre affinity
            commd = f"cd /home/ab/ANTIPASTI-main/notebooks; conda run -n antipasti-env python pre_affinity.py -i {out_path}/abdockgen_{pro_name}/process_pdb/{tmp_id}_chothia.pdb -n {tmp_id} -o {out_path}/abdockgen_{pro_name}/process_pdb/abdockgen_affinity.csv -d 0000"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ###read res
            tmp_aff = pd.read_csv(f"{out_path}/abdockgen_{pro_name}/process_pdb/abdockgen_affinity.csv")
            tmp_aff = tmp_aff.iloc[:, -2:]
            affinity_res = pd.concat([affinity_res, tmp_aff], axis=0, ignore_index=True)
        #affinity_res.to_csv(f"{out_path}/abdockgen_{pro_name}/gen_pdb/abdockgen_affinity.csv")
    if binding_energy:
        print("Predicting binding energy ...")
        all_pdbs = os.listdir(f"{out_path}/abdockgen_{pro_name}/gen_pdb")
        os.makedirs(f"{out_path}/abdockgen_{pro_name}/rosetta_pdb/", exist_ok=True)
        ids = []
        dg = []
        for i in all_pdbs:
            full_path = f"{out_path}/abdockgen_{pro_name}/gen_pdb/" + i
            tmp_id = re.sub(r'.pdb',r'',i)
            ##remove 0
            remove_zero_coords(full_path, f"{out_path}/abdockgen_{pro_name}/rosetta_pdb/{tmp_id}.pdb")
            commd = f"cd /home/ab/rosetta/run; ./score_jd2.static.linuxgccrelease -s {out_path}/abdockgen_{pro_name}/rosetta_pdb/{tmp_id}.pdb -no_optH false -ignore_unrecognized_res -out:pdb -out:path:all {out_path}/abdockgen_{pro_name}/rosetta_pdb/"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            commd = f"cd /home/ab/rosetta/run; ./InterfaceAnalyzer.static.linuxgccrelease -s {out_path}/abdockgen_{pro_name}/rosetta_pdb/{tmp_id}_0001.pdb @./pack_input_options.txt -interface H_A -out:path:all {out_path}/abdockgen_{pro_name}/rosetta_pdb/"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ###read score
            commd = "awk '{print $6}' " + f"{out_path}/abdockgen_{pro_name}/rosetta_pdb/pack_input_score.sc | sed '1d' > {out_path}/abdockgen_{pro_name}/rosetta_pdb/outscore"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            tmp_score = pd.read_csv(f"{out_path}/abdockgen_{pro_name}/rosetta_pdb/outscore")
            ids.append(tmp_id)
            dg.append(tmp_score.dG_separated.tolist()[0])
            ##clean
            commd = f"rm {out_path}/abdockgen_{pro_name}/rosetta_pdb/*"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
        energy_res = pd.DataFrame({"id":ids,"dg":dg})
        energy_res.to_csv(f"{out_path}/abdockgen_{pro_name}/gen_pdb/abdockgen_rosetta.csv")
    if binding_affinity:
        affinity_res.to_csv(f"{out_path}/abdockgen_{pro_name}/gen_pdb/abdockgen_affinity.csv")
    ##zip res
    shutil.make_archive(f"{out_path}/abdockgen_gen_{pro_name}", "zip", f"{out_path}/abdockgen_{pro_name}/gen_pdb/")
    commd = f"rm -rf {out_path}/abdockgen_{pro_name}"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    print("Done !")

if "RFantibody" in design_method:
    print("Run RFantibody ...")
    ##pre-process
    pdb_id = re.sub(r'.+/',r'',re.sub(r'.pdb', r'',pdb_file))
    os.makedirs(f"{out_path}/rfantibody_{pro_name}/", exist_ok=True)
    commd = f"python get_epi_trunc_pdb.py -p {pdb_file} -a {hchain} -b {agchain} -o {out_path}/rfantibody_{pro_name}/{pdb_id}_trunc.pdb -e {out_path}/rfantibody_{pro_name}/{pdb_id}_epi.txt"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ##read epi info and hcdr3 info
    with open(f'{out_path}/rfantibody_{pro_name}/{pdb_id}_epi.txt', 'r') as file:
        content = file.read()
        my_list = eval(content)
    epi_info = str(my_list[:-1])
    epi_info = "ppi.hotspot_res="+re.sub(r"'",r"",epi_info)
    hcdr3_info = "antibody.design_loops=[H3:"+my_list[-1]+"]"

    commd = f"pdb_selchain -{hchain} {pdb_file} > {out_path}/rfantibody_{pro_name}/target_H.pdb"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ##run
    os.makedirs(f"{out_path}/rfantibody_{pro_name}/gen_pdb/", exist_ok=True)
    os.makedirs(f"{out_path}/rfantibody_{pro_name}/gen_seq/", exist_ok=True)
    os.makedirs(f"{out_path}/rfantibody_{pro_name}/gen_af2/", exist_ok=True)
    commd = f"docker cp {out_path}/rfantibody_{pro_name}/ rfantibody_wt:/data/"
    subprocess.run(commd, shell=True, text=True, capture_output=False) ###move data to docker
    commd = f"docker exec rfantibody_wt poetry run python /home/scripts/util/chothia2HLT.py --input_pdb /data/rfantibody_{pro_name}/target_H.pdb --heavy {hchain} --output /data/rfantibody_{pro_name}/target_H_HLT.pdb"
    subprocess.run(commd, shell=True, text=True, capture_output=False) ###annote h chain pdb
    ###run step 1
    commd = f"docker exec rfantibody_wt poetry run python  /home/scripts/rfdiffusion_inference.py --config-name antibody antibody.target_pdb=/data/rfantibody_{pro_name}/{pdb_id}_trunc.pdb antibody.framework_pdb=/data/rfantibody_{pro_name}/target_H_HLT.pdb inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt '{epi_info}' '{hcdr3_info}' inference.num_designs={design_num} inference.final_step=48 inference.deterministic=True diffuser.T=50 inference.output_prefix=/data/rfantibody_{pro_name}/gen_pdb/gen"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ###run step2
    commd = f"docker exec rfantibody_wt poetry run python /home/scripts/proteinmpnn_interface_design.py -pdbdir /data/rfantibody_{pro_name}/gen_pdb/ -outpdbdir /data/rfantibody_{pro_name}/gen_seq/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ###run step3
    commd = f"docker exec rfantibody_wt poetry run python /home/scripts/rf2_predict.py input.pdb_dir=/data/rfantibody_{pro_name}/gen_seq/ output.pdb_dir=/data/rfantibody_{pro_name}/gen_af2/ model.model_weights=/home/weights/RF2_ab.pt"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ###move out results
    commd = f"docker cp rfantibody_wt:/data/rfantibody_{pro_name}/gen_pdb {out_path}/rfantibody_{pro_name}/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    commd = f"docker cp rfantibody_wt:/data/rfantibody_{pro_name}/gen_seq {out_path}/rfantibody_{pro_name}/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    commd = f"docker cp rfantibody_wt:/data/rfantibody_{pro_name}/gen_af2 {out_path}/rfantibody_{pro_name}/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ###read res
    all_gen = os.listdir(f"{out_path}/rfantibody_{pro_name}/gen_seq/")
    pdbparser = PDBParser(QUIET=True)
    sequence = []
    seqid = []
    for i in all_gen:
        full_path = os.path.join(f"{out_path}/rfantibody_{pro_name}/gen_seq/", i)
        structure = pdbparser.get_structure("protein", full_path)
        chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}
        sequence.append(chains["H"])
        seqid.append(re.sub(r'_dldesign_0.pdb', r'', i))
    rfantibody_res = pd.DataFrame({"ID":seqid,"Sequence":sequence,"Method":"RFantibody"})
    res_all = pd.concat([res_all, rfantibody_res], axis=0, ignore_index=True)
    ###score
    print("Collect AF2 scores ...")
    all_af2 = os.listdir(f"{out_path}/rfantibody_{pro_name}/gen_af2/")
    seqid = []
    seq_pae = []
    seq_rmsd = []
    for i in all_af2:
        full_path = os.path.join(f"{out_path}/rfantibody_{pro_name}/gen_af2/", i)
        ll = find_start_line(full_path,"SCORE")
        ll_score = pd.read_csv(full_path, skiprows=ll,sep=": ",engine='python', names=["type","score"])
        pae = ll_score[ll_score['type'].str.contains('SCORE pae')].score.tolist()[0]
        rmsd = ll_score[ll_score['type'].str.contains('H3_rmsd')].score.tolist()[0]
        seqid.append(re.sub(r'_dldesign_0_best.pdb', r'', i))
        seq_pae.append(pae)
        seq_rmsd.append(rmsd)
    rfantibody_score = pd.DataFrame({"ID":seqid, "PAE":seq_pae, "RMSD":seq_rmsd})
    if binding_energy:
        print("Predicting binding energy ...")
        all_pdbs = os.listdir(f"{out_path}/rfantibody_{pro_name}/gen_seq")
        os.makedirs(f"{out_path}/rfantibody_{pro_name}/rosetta_pdb/", exist_ok=True)
        ids = []
        dg = []
        for i in all_pdbs:
            full_path = f"{out_path}/rfantibody_{pro_name}/gen_seq/" + i
            tmp_id = re.sub(r'.pdb',r'',i)
            commd = f"cd /home/ab/rosetta/run; ./score_jd2.static.linuxgccrelease -s {full_path} -no_optH false -ignore_unrecognized_res -out:pdb -out:path:all {out_path}/rfantibody_{pro_name}/rosetta_pdb/"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            commd = f"cd /home/ab/rosetta/run; ./InterfaceAnalyzer.static.linuxgccrelease -s {out_path}/rfantibody_{pro_name}/rosetta_pdb/{tmp_id}_0001.pdb @./pack_input_options.txt -interface H_T -out:path:all {out_path}/rfantibody_{pro_name}/rosetta_pdb/"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            ###read score
            commd = "awk '{print $6}' " + f"{out_path}/rfantibody_{pro_name}/rosetta_pdb/pack_input_score.sc | sed '1d' > {out_path}/rfantibody_{pro_name}/rosetta_pdb/outscore"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
            tmp_score = pd.read_csv(f"{out_path}/rfantibody_{pro_name}/rosetta_pdb/outscore")
            ids.append(tmp_id)
            dg.append(tmp_score.dG_separated.tolist()[0])
            ##clean
            commd = f"rm {out_path}/rfantibody_{pro_name}/rosetta_pdb/*"
            subprocess.run(commd, shell=True, text=True, capture_output=False)
        energy_res = pd.DataFrame({"id":ids,"dg":dg})
        energy_res.to_csv(f"{out_path}/rfantibody_{pro_name}/gen_seq/rfantibody_rosetta.csv")
    rfantibody_score.to_csv(f"{out_path}/rfantibody_{pro_name}/gen_seq/rfantibody_scores.csv")
    ###zip res
    shutil.make_archive(f"{out_path}/rfantibody_gen_{pro_name}", "zip", f"{out_path}/rfantibody_{pro_name}/gen_seq/")
    commd = f"rm -rf {out_path}/rfantibody_{pro_name}/"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    ###rm tmp
    commd = f"docker exec -it rfantibody_wt rm -rf /data/rfantibody_{pro_name}"
    subprocess.run(commd, shell=True, text=True, capture_output=False)
    print("Done !")

###
res_all.to_csv(f"{out_path}/struc_gen_seq.csv")
