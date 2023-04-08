#!/bin/bash

###################################################################################
#	unzip fastqc.gz files
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"
	clean_dir="${out_root}/clean_fastq"

        # Outputs
	unzip_dir="${clean_dir}/unzip"
	mkdir -p $unzip_dir

	cd ${clean_dir}
	fqs=(`ls | grep '.fastq.gz$'`)

	# unzip files
	#fq="BF10_S16_cont_removed_R1.fastq.gz"
	for fq in ${fqs[@]}; do
	    cmd="gunzip -k $fq"
	    echo $cmd
	    eval $cmd
	done

