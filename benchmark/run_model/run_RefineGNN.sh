while IFS=" " read -r pdbid hchain agchain
do
	mkdir /home/data/sdb/wt/antibody_design/benchmark/RefineGNN/gen/${pdbid}
	cd ~/RefineGNN
	python gen.py --data_path /home/data/sdb/wt/antibody_design/benchmark/RefineGNN/processed_pdb/${pdbid}.jsonl --load_model ckpts/RefineGNN-rabd/model.best --out /home/data/sdb/wt/antibody_design/benchmark/RefineGNN/gen/${pdbid}
done < "pdbs.txt"
