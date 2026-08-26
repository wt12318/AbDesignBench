#!/usr/bin/bash
ROSETTA3=/home/data/sdb/wt/antibody_design/rosetta/rosetta.binary.ubuntu.release-371/main/source
while IFS=" " read -r fileid pdbid hchain agchain
do
        mkdir ${fileid}
        cd ${fileid}
        $ROSETTA3/bin/score_jd2.static.linuxgccrelease -s ../../gen_extract_HCDR3/${fileid}.pdb -no_optH false -ignore_unrecognized_res -out:pdb
        mv ${fileid}_0001.pdb ${fileid}_scored.pdb
        $ROSETTA3/bin/InterfaceAnalyzer.static.linuxgccrelease -s ${fileid}_scored.pdb @../pack_input_options.txt -interface ${hchain}_${agchain}
        mv pack_input_score.sc ../ROSETTA_energy/${fileid}_pack_input_score.sc
        mv score.sc ../ROSETTA_energy/${fileid}_score.sc
        mv ${fileid}_scored.pdb ../ROSETTA_energy/
        cd ../
	\rm -rf ${fileid}
done < "xaa"
