#!/bin/bash

###################################################################################
#	Binn QC using CheckM
##!!!!! Run after "conda activate checkm"
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"

	clean_dir="${out_root}/clean_fastq"
	unzip_dir="${clean_dir}/unzip"
        assem_dir="${out_root}/assembly"
        assemmap_dir="${assem_dir}/map_reads2assembly"
	assembam_dir="${assemmap_dir}/bam"
	
	# Output
	bin_dir="${out_root}/bin/nre100_bins"
        assemdepth_dir="/home4/sjeong6/Paerl/output2/checkm_nre100/depth"
	binqc_dir="/home4/sjeong6/Paerl/output2/checkm_nre100/QC"
	mkdir -p $assemdepth_dir $binqc_dir
	
        # Tool
	#checkm="/home4/sjeong6/miniconda3/envs/checkm/bin/checkm"
#       checkm data setRoot /home4/sjeong6/tools/checkm_DB
	


	# Bin QC using checkM
	cmd="checkm lineage_wf -x fa ${bin_dir} ${binqc_dir} --tab_table -f ${binqc_dir}/MAGs_checkm.tab --reduced_tree "

	echo $cmd
	eval $cmd

	
