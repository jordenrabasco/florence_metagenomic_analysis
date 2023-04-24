#!/bin/bash

###################################################################################
#	MAG taxonomic classification using GTDB-Tk
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
        bin_dir="${out_root}/bin"
	
	# Output
	binclass_dir="${out_root}/bin_class"

	
        # Tool
	gtdbtk="/home4/sjeong6/miniconda3/envs/gtdbtk-2.1.1/bin/gtdbtk"

        # column -s, -t < ${dat_root}/meta/NGS_WorkOrder_NRE_MetaG_12_2_21_JS_Sample_Names.csv # | less -#2 -N -S
        nre30=("N411_S21" "BF11_S17" "BF3_S9" "N591_S26" "BF15_S35" "BF7_S13" "N824_S31")
        nre70=("N416_S22" "BF12_S18" "BF4_S10" "N596_S27" "N613_S30" "BF8_S14" "N830_S32")
        nre100=("N419_S23" "BF13_S19" "BF5_S11" "N599_S28" "BF17_S36" "BF9_S15" "N833_S33")
        nre180=("N424_S24" "N425_S25" "BF14_S20" "BF6_S12" "N605_S29" "BF18_S37" "BF10_S16" "N839_S34")
	
	nres=("nre30" "nre70" "nre100" "nre180")

	cd ${bin_dir}
	cmd="${gtdbtk} classify_wf --extension fa --genome_dir ${bin_dir} --out_dir ${binclass_dir}"
	echo $cmd
	eval $cmd
	
