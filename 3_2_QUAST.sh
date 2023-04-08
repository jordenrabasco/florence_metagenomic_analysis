#!/bin/bash

###################################################################################
#	QC Assembly 
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"

	clean_dir="${out_root}/clean_fastq"
	unzip_dir="${clean_dir}/unzip"
        assem_dir="${out_root}/assembly"
	
        # Outputs
	assemqc_dir="${assem_dir}/assembly_QC"
	mkdir -p $assemqc_dir
	
        # Tool
	module add python
	quast="/home4/sjeong6/tools/quast-5.2.0/quast.py"
#	quast="/home4/sjeong6/miniconda3/pkgs/quast-5.0.2-py27pl5262h8eb80aa_5/bin/metaquast.py"

	nre=("nre30" "nre70" "nre100" "nre180")
	
	# Run metaquast for assembly files together
	cd ${assem_dir}
	fa="_assembled.contigs.fa"

	cmd="${quast} -m 50 -t 10 ./${nre[0]}/${nre[0]}${fa} ./${nre[1]}/${nre[1]}${fa} ./${nre[2]}/${nre[2]}${fa} ./${nre[3]}/${nre[3]}${fa} -o $assemqc_dir"
	echo $cmd
	eval $cmd

