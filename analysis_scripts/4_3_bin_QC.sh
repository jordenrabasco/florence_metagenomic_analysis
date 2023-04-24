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
	bin_dir="${out_root}/bin"
        assemdepth_dir="/home4/sjeong6/Paerl/output2/test_checkm/depth"
	binqc_dir="/home4/sjeong6/Paerl/output2/test_checkm/QC"
	mkdir -p $assemdepth_dir $binqc_dir
	
        # Tool
	#checkm="/home4/sjeong6/miniconda3/envs/checkm/bin/checkm"
#       checkm data setRoot /home4/sjeong6/tools/checkm_DB
	#export CHECKM_DATA_PATH=/home4/sjeong6/tools/checkm_DB
	
        # column -s, -t < ${dat_root}/meta/NGS_WorkOrder_NRE_MetaG_12_2_21_JS_Sample_Names.csv # | less -#2 -N -S
        nre30=("N411_S21" "BF11_S17" "BF3_S9" "N591_S26" "BF15_S35" "BF7_S13" "N824_S31")
        nre70=("N416_S22" "BF12_S18" "BF4_S10" "N596_S27" "N613_S30" "BF8_S14" "N830_S32")
        nre100=("N419_S23" "BF13_S19" "BF5_S11" "N599_S28" "BF17_S36" "BF9_S15" "N833_S33")
        nre180=("N424_S24" "N425_S25" "BF14_S20" "BF6_S12" "N605_S29" "BF18_S37" "BF10_S16" "N839_S34")
	
	nres=("nre30" "nre70" "nre100" "nre180")

	# Bin QC using checkM
	cmd="checkm lineage_wf -x fa ${bin_dir} ${binqc_dir} --tab_table -f ${binqc_dir}/MAGs_checkm.tab --reduced_tree "


	#cmd="checkm lineage_wf -g -x fna -e 1e-10 -l 0.7 -f ${binqc_dir}/MAGs_checkm.tab --tab_table --reduced_tree ${bin_dir} ${binqc_dir}"
#	cmd="${checkm} lineage_wf -x fa ${bin_dir} ${binqc_dir} --tab_table -f ${binqc_dir}/MAGs_checkm.tab --reduced_tree"
	echo $cmd
	eval $cmd

	
