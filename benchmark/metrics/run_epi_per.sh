while IFS=" " read -r gen_file ori_file hchain agchain
do
	python cal_epi_per.py extract_HCDR3/${ori_file} gen_extract_HCDR3/${gen_file} ${hchain} ${agchain} ${hchain} ${agchain} -c 8 -o epi_per_results.csv
done < "mapping_ids"
