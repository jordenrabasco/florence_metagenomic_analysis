#!/bin/bash

###################################################################################
#	Remove host contamination (human) using Bowtie2
#       1) Build host genome NCBI db (GRCh38)
#       2) Map reads against host database
#       3) Convert sam to bam file using samtools
#       4) Extract unmapped reads (non-host)
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"
        trim_dir="${out_root}/trim"
	
        # Outputs
	filter_dir="${out_root}/rm_cont"
	bam_dir="${filter_dir}/bam"
	mkdir -p $bam_dir
	clean_dir="${out_root}/clean_fastq"
	single_dir="${clean_dir}/singleton"
	mkdir -p $single_dir
		
        # Tool
	bowtie2="/opt/bowtie2/2.4.4/bowtie2"
	bowtie2_build="/opt/bowtie2/2.4.4/bowtie2-build"
	samtools="/opt/samtools/1.13/bin/samtools"
	
	# 1) Download GRCh38 db
	cd ${filter_dir}
#	wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
#	unzip GRCh38_noalt_as.zip
	
	cd ${trim_dir}
       	r1s=(`ls | grep '_R1_trimmed.fastq.gz$'`)

	#r1="BF10_S16_R1_trimmed.fastq.gz"
	for r1 in ${r1s[@]}; do
	    pre=`echo $r1 | sed 's/_R1_trimmed.fastq.gz//'`
	    r2=`echo $r1 | sed 's/_R1_/_R2_/'`
	    echo $pre
	    ls $r1
	    ls $r2

	    # 2) Map reads against GRCh38 db  
	    cmd="${bowtie2} --sensitive-local -p 8 -x ${filter_dir}/GRCh38_noalt_as/GRCh38_noalt_as -1 $r1 -2 $r2 -S ${filter_dir}/${pre}_mapped_2human.sam 2> ${filter_dir}/${pre}_mapped_2human.log"
	    echo $cmd
	    eval $cmd
	    
            # 3) Convert sam to bam file using samtools 
	    cmd="${samtools} view -@ 8 -bS ${filter_dir}/${pre}_mapped_2human.sam > ${bam_dir}/${pre}_mapped_2human.bam"
	    echo $cmd
	    eval $cmd
	    
	    # 4) Extract unmapped reads (non-host)
	    cmd="${samtools} view -@ 8 -b -f 12 -F 256 ${bam_dir}/${pre}_mapped_2human.bam > ${bam_dir}/${pre}_unmapped_2human.bam"
	    echo $cmd
	    eval $cmd

	    # 5) Sort BAM file to organize paired reads
	    cmd="${samtools} sort -n ${bam_dir}/${pre}_unmapped_2human.bam -o ${bam_dir}/${pre}_unmapped_sorted.bam"
	    echo $cmd
	    eval $cmd

	    # 6) Convert BAM to fastq
	    cmd="${samtools} fastq -@ 8 ${bam_dir}/${pre}_unmapped_sorted.bam \
	    		     -1 ${clean_dir}/${pre}_cont_removed_R1.fastq.gz -2 ${clean_dir}/${pre}_cont_removed_R2.fastq.gz \
			     -0 ${clean_dir}/${pre}_other.fastq.gz -s ${single_dir}/${pre}_singleton.fastq.gz -n"
	    echo $cmd
	    eval $cmd

	done
