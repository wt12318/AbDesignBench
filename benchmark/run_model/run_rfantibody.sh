docker run -itd --name rfantibody_wt --gpus all --memory 100g rfantibody_last bash
docker exec rfantibody_wt poetry run python /home/scripts/util/chothia2HLT.py --input_pdb /data/rfantibody_{pro_name}/target_H.pdb --heavy {hchain} --output /data/rfantibody_{pro_name}/target_H_HLT.pdb
##step1
docker exec rfantibody_wt poetry run python  /home/scripts/rfdiffusion_inference.py --config-name antibody antibody.target_pdb=/data/rfantibody_{pro_name}/{pdb_id}_trunc.pdb antibody.framework_pdb=/data/rfantibody_{pro_name}/target_H_HLT.pdb inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt '{epi_info}' '{hcdr3_info}' inference.num_designs={design_num} inference.final_step=48 inference.deterministic=True diffuser.T=50 inference.output_prefix=/data/rfantibody_{pro_name}/gen_pdb/gen
##step2
docker exec rfantibody_wt poetry run python /home/scripts/proteinmpnn_interface_design.py -pdbdir /data/rfantibody_{pro_name}/gen_pdb/ -outpdbdir /data/rfantibody_{pro_name}/gen_seq/
##step3
docker exec rfantibody_wt poetry run python /home/scripts/rf2_predict.py input.pdb_dir=/data/rfantibody_{pro_name}/gen_seq/ output.pdb_dir=/data/rfantibody_{pro_name}/gen_af2/ model.model_weights=/home/weights/RF2_ab.pt
