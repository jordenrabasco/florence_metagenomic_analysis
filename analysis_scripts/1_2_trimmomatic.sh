#!/bin/bash

###################################################################################
#	Trim adapters and filter the raw reads using Trimmomatic
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"

        # Outputs
	trim_dir="${out_root}/trim"
	logstrim_dir="${out_root}/logs_trim"

        # Tool
	trimmomatic="/opt/Trimmomatic/0.39/trimmomatic-0.39.jar"

	cd ${dat_root}/NVS139B_fastq
	r1s=(`ls | grep '_R1_001.fastq.gz$'`)

	# run fastqc
	#r1="BF10_S16_L001_R1_001.fastq.gz"
	for r1 in ${r1s[@]}; do
	    pre=`echo $r1 | sed 's/_L001_R1_001.fastq.gz//'`
	    r2=`echo $r1 | sed 's/_R1_/_R2_/'`
	    echo $pre
	    ls $r1
	    ls $r2
	    
	    cmd="java -jar ${trimmomatic} PE -threads 4 -trimlog ${logstrim_dir}/${pre}_trimmed.log ${r1} ${r2} \
	    	           ${trim_dir}/${pre}_R1_trimmed.fastq.gz ${trim_dir}/${pre}_R1_unpaired.fastq.gz \
			   ${trim_dir}/${pre}_R2_trimmed.fastq.gz ${trim_dir}/${pre}_R2_unpaired.fastq.gz \
			   ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50"
	    echo $cmd
	    eval $cmd
	done
