while IFS=" " read -r origin_seq leading_seq seqid min_len
do
	abgpt_generate --chain_type heavy --starting_residue $leading_seq --num_seqs 200 --min_len $min_len;
	mv /home/data/sdb/wt/antibody_design/abgpt/bcr_design/heavy_${leading_seq}.txt /home/data/sdb/wt/antibody_design/benchmark/AbGPT/gen2/$seqid.gen1;
	abgpt_generate --chain_type heavy --starting_residue $leading_seq --num_seqs 200 --min_len $min_len;
	mv /home/data/sdb/wt/antibody_design/abgpt/bcr_design/heavy_${leading_seq}.txt /home/data/sdb/wt/antibody_design/benchmark/AbGPT/gen2/$seqid.gen2;
	abgpt_generate --chain_type heavy --starting_residue $leading_seq --num_seqs 200 --min_len $min_len;
	mv /home/data/sdb/wt/antibody_design/abgpt/bcr_design/heavy_${leading_seq}.txt /home/data/sdb/wt/antibody_design/benchmark/AbGPT/gen2/$seqid.gen3;
	abgpt_generate --chain_type heavy --starting_residue $leading_seq --num_seqs 200 --min_len $min_len;
	mv /home/data/sdb/wt/antibody_design/abgpt/bcr_design/heavy_${leading_seq}.txt /home/data/sdb/wt/antibody_design/benchmark/AbGPT/gen2/$seqid.gen4;
	abgpt_generate --chain_type heavy --starting_residue $leading_seq --num_seqs 200 --min_len $min_len;
	mv /home/data/sdb/wt/antibody_design/abgpt/bcr_design/heavy_${leading_seq}.txt /home/data/sdb/wt/antibody_design/benchmark/AbGPT/gen2/$seqid.gen5;
	cd /home/data/sdb/wt/antibody_design/benchmark/AbGPT;
	cat gen2/$seqid.gen1 gen2/$seqid.gen2 gen2/$seqid.gen3 gen2/$seqid.gen4 gen2/$seqid.gen5 > gen2/$seqid.gen;
	\rm gen2/$seqid.gen1 gen2/$seqid.gen2 gen2/$seqid.gen3 gen2/$seqid.gen4 gen2/$seqid.gen5;
	cd /home/data/sdb/wt/antibody_design/abgpt/
done < "/home/data/sdb/wt/antibody_design/benchmark/AbGPT/allseqs.txt"
