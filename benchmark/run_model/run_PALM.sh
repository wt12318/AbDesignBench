while IFS=" " read -r agseq Hseq start end cdr3len seqid
do
	sed 's/HEAVY/'$Hseq'/g' config/common/seq2seq_generate.json > tmp.json
	sed -i 's/CDH3_B/'$start'/g' tmp.json
	sed -i 's/CDH3_E/'$end'/g' tmp.json
	sed -i 's/ANTIGEN/'$agseq'/g' tmp.json
	sed -i 's/CDR3LEN/'$cdr3len'/g' tmp.json
	python generate_antibody.py --config tmp.json
	mv ../Result_seq2seq_gen/datasplit/CoV_AbDab-Seq2seq-Evaluate-Common/*/result.csv /home/data/sdb/wt/antibody_design/benchmark/PALM/gen/$seqid.gen
	\rm -rf ../Result_seq2seq_gen/datasplit/CoV_AbDab-Seq2seq-Evaluate-Common/*
	\rm -rf tmp.json
done < "/home/data/sdb/wt/antibody_design/benchmark/PALM/allseqs.txt"
