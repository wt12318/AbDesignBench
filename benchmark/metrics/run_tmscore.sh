while IFS=" " read -r gen_file ori_file
do
	TMscore ${ori_file} ${gen_file} -outfmt 2 > TMscore_out/${gen_file}.out
done < "xaa"
