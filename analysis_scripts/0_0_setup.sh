#!/bin/bash

###################################################################################
#	set up folder structure
###################################################################################

        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/output2"

	######
	# 1. set up output sub-folders, in the order of output
	trim_dir="trim"
	logstrim_dir="logs_trim"
	filter_dir="rm_cont"
	clean_dir="clean_fastq"
	pool_dir="pooled" # all samples combined
	qc_dir="QC"     # fastqc for raw read fastq files
	qc2_dir="QC2"   # fastqc for cleaned fastq files
	kdb_dir="kaiju_db"
	pretaxa_dir="pre_taxa" # taxonomic assignment for reads 
	assem_dir="assembly"   # using metahit
#	assem2_dir="assembly2" # using metaSPAdes
	bin_dir="bin"
	binclass_dir="bin_class"
	bintaxa_dir="bin_taxa"
		
#	mkdir -p $out_root && cd $out_root
#	mkdir -p $log_dir $logstrim_dir $merge_dir $trim_dir $filter_dir $clean_dir $pool_dir $qc_dir $qc2_dir $kdb_dir $pretaxa_dir $assem_dir $assem2_dir $bin_dir $binclass_dir $bintaxa_dir


	################# INSTALLATION #################
	tool_dir="/home4/sjeong6/tools"
	cd ${tool_dir}

	# metabat2
#	conda install -c bioconda metabat2
#	conda install -c bioconda/label/cf201901 metabat2

	# spades
#	wget http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Linux.tar.gz
#	tar -xzf SPAdes-3.15.4-Linux.tar.gz
#	cd SPAdes-3.15.4-Linux/bin/

	# for checkm
#	conda create -n checkm python=2.7
#	conda activate checkm
#	conda install -c bioconda numpy matplotlib pysam
#	conda install -c bioconda hmmer prodigal pplacer pysam
	pip3 install checkm-genome

	conda install -c bioconda checkm-genome
#	mkdir -p ${tool_dir}/checkm_DB
#	tar -xvzf checkm_data_2015_01_16.tar.gz -C "${tool_dir}/checkm_DB" --strip 1 > /dev/null
#	checkm data setRoot "${tool_dir}/checkm_DB"
	export CHECKM_DATA_PATH=/home4/sjeong6/tools/checkm_DB
	
	# GTDB-Tk
#        conda create -n gtdbtk-2.1.1 -c conda-forge -c bioconda gtdbtk=2.1.1
#	conda activate gtdbtk-2.1.1

	cd $tool_dir
#	wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz
#	mkdir -p ${tool_dir}/gtdb_tk_DB
#	tar -xvzf gtdbtk_r207_v2_data.tar.gz -C "/home4/sjeong6/tools/gtdb_tk_DB" --strip 1 > /dev/null
#	rm gtdbtk_r207_v2_data.tar.gz
#	conda env config vars set GTDBTK_DATA_PATH="/home4/sjeong6/tools/gtdb_tk_DB"
	
	
