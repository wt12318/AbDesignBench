while IFS=" " read -r pdbid hchain agchain
do
	cd ~/abdockgen
	python process_data_for_RefineGNN.py --path /home/data/sdb/wt/antibody_design/benchmark/pdb_files/ --pdb ${pdbid}.pdb --hchain ${hchain} --out /home/data/sdb/wt/antibody_design/benchmark/RefineGNN/processed_pdb/${pdbid}.jsonl
done < "pdbs.txt"
