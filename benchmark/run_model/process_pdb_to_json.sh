while IFS=" " read -r pdbid hchain agchain
do
	python /home/wt/abdockgen/process_data.py ../pdb_files/$pdbid.pdb $hchain $agchain > processed_pdb/$pdbid.json 
done < "pdbs.txt"
