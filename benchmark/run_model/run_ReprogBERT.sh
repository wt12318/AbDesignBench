while IFS=" " read -r origin_seq mask_seq seqid
do
	python main.py  --run inference --single_input $mask_seq --model_type base --exp_dir base_cd3 --checkpoint /home/wt/ReprogBERT/output/base_cd3/checkpoints/last.ckpt --progen_dir /home/data/sdb/wt/antibody_design/ReprogBERT --num_samples 1000;
	mv /home/wt/ReprogBERT/output/base_cd3/inference_smpl.fasta /home/data/sdb/wt/antibody_design/benchmark/ReprogBERT/gen/$seqid.gen
done < "/home/data/sdb/wt/antibody_design/benchmark/ReprogBERT/allseqs.txt"
