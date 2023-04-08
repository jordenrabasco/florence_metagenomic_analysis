#!/bin/bash

###################################################################################
#	fastqc for raw reads
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"

        # Outputs
	qc_dir="${out_root}/QC"

        # Tool
	fastqc="/opt/fastqc/0.11.9/fastqc"

	cd ${dat_root}/NVS139B_fastq
	fqs=(`ls | grep '.fastq.gz'`)

	# run fastqc
	#fq="BF10_S16_L001_R1_001.fastq.gz"
	for fq in ${fqs[@]}; do
	    cmd="$fastqc $fq -o ${qc_dir}"
	    echo $cmd
	    eval $cmd
	done
