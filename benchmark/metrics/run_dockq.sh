while IFS=" " read -r gen_file ori_file hchain agchain
do
	DockQ --mapping ${hchain}${agchain}:${hchain}${agchain} ../gen_extract_HCDR3/${gen_file} ../extract_HCDR3/${ori_file} --json DockQ_out/${gen_file}.json
done < "xaa"
