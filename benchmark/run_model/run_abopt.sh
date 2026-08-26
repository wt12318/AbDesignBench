while IFS=" " read -r pdbid hchain agchain
do
	cd ~/ab_opt/AbDock
	mkdir /home/data/sdb/wt/antibody_design/benchmark/ab_opt/gen/gen_${pdbid}_${hchain}_${agchain}
	python dock_pdb.py --heavy $hchain --pdb_path /home/data/sdb/wt/antibody_design/benchmark/ab_opt/processed_pdb/${pdbid}_${hchain}_${agchain}.pdb --config configs/test/dock_cdr.yml --ckpt reproduction/dock_single_cdr/250000.pt -n 100 -b 100 -d cuda -o /home/data/sdb/wt/antibody_design/benchmark/ab_opt/gen/gen_${pdbid}_${hchain}_${agchain}/dock_single
	python optimize_ab.py --num_gpus 1 --process_per_gpu 8 --docked_pose_dir /home/data/sdb/wt/antibody_design/benchmark/ab_opt/gen/gen_${pdbid}_${hchain}_${agchain}/dock_single/dock_cdr/${pdbid}_${hchain}_${agchain}.pdb_/H_CDR3/ --seq_design_dir /home/data/sdb/wt/antibody_design/benchmark/ab_opt/gen/gen_${pdbid}_${hchain}_${agchain}/seq_design --design_contig '' --screen_dir /home/data/sdb/wt/antibody_design/benchmark/ab_opt/gen/gen_${pdbid}_${hchain}_${agchain}/screen --heavy_chain_id "${hchain}"  --nums 1000
done < "xaa"
