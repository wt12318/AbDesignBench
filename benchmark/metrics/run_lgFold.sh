while IFS=" " read -r ids gen_seq lseq
do
	python pre_pdb_lgFold.py -p ./pdbs/${ids}.pdb -a $gen_seq -l $lseq
done < "xab"
