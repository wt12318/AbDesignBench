while IFS=" " read -r pdbid
do
	mkdir /home/data/sdb/wt/antibody_design/benchmark/abdockgen/gen/${pdbid}
	cd /home/wt/abdockgen
	python gen_seq.py --json /home/data/sdb/wt/antibody_design/benchmark/abdockgen/processed_pdb/${pdbid}.json --batch_size 200 --num 1000 --out /home/data/sdb/wt/antibody_design/benchmark/abdockgen/gen/${pdbid}/
done < "remain_ids"
