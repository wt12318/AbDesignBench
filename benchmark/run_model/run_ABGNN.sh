while IFS=" " read -r pdbid
do
	cd /home/data/sdb/wt/antibody_design/ABGNN/ABGNN
	mkdir data/gen
	cp /home/data/sdb/wt/antibody_design/benchmark/ABGNN/processed_pdb/$pdbid.json  data/gen/test_data.jsonl
	python gen_structure_seq.py --cdr_type 3 --cktpath checkpoints/exp2/pflen5_iter5_loss1_1_2_lr0.0001_bsz2_seed128/checkpoint_best.pt --data_path data/gen --savepdb data/gen/${pdbid}_gen.pdb --savemetrics data/gen/${pdbid}_gen.csv --num_decode 1000
	mv data/gen/${pdbid}_gen.pdb /home/data/sdb/wt/antibody_design/benchmark/ABGNN/gen/
	mv data/gen/${pdbid}_gen.csv /home/data/sdb/wt/antibody_design/benchmark/ABGNN/gen/
	\rm -rf data/gen
done < "pdbs_ids.txt"
