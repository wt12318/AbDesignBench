while IFS=" " read -r pdbid hchain agchain
do
	cd /home/wt/diffab
	python design_pdb.py --heavy $hchain -c configs/test/codesign_single.yml -o /home/data/sdb/wt/antibody_design/benchmark/DiffAb/gen/ /home/data/sdb/wt/antibody_design/benchmark/DiffAb/processed_pdb/${pdbid}_${hchain}_${agchain}.pdb
done < "pdbs.txt"
