while IFS=" " read -r seqs sp start end seqid
do
	iglm_infill fasta/$seqid.fasta :H $start $end --chain_token [HEAVY] --species_token [$sp] --num_seqs 1000 --output_dir ./gen/;
	mv ./gen/generated_seqs.fasta ./gen/$seqid.gen
done < "allseqs.txt"
