#!/bin/bash

###################################################################################
#	fastqc for raw reads
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"
	clean_dir="${out_root}/clean_fastq"

        # Outputs
	qc2_dir="${out_root}/QC2"

        # Tool
	fastqc="/opt/fastqc/0.11.9/fastqc"

	cd ${clean_dir}
	fqs=(`ls | grep '.fastq.gz'`)
	
	# run fastqc
	#fq="BF10_S16_cont_removed_R1.fastq.gz"
	for fq in ${fqs[@]}; do
	    cmd="$fastqc $fq -o ${qc2_dir}"
	    echo $cmd
	    eval $cmd
	done
