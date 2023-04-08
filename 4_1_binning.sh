#!/bin/bash

###################################################################################
#	Binning for MAGs 
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
        assemdepth_dir="${assem_dir}/depth"
	mkdir -p $assemdepth_dir	
	
        # Tool
	metabat2="/home4/sjeong6/miniconda3/pkgs/metabat2-2.15-h137b6e9_0/bin/metabat2"
	jgi_summarize_bam_contig_depths="/home4/sjeong6/miniconda3/pkgs/metabat2-2.15-h137b6e9_0/bin/jgi_summarize_bam_contig_depths"

        # column -s, -t < ${dat_root}/meta/NGS_WorkOrder_NRE_MetaG_12_2_21_JS_Sample_Names.csv # | less -#2 -N -S
        nre30=("N411_S21" "BF11_S17" "BF3_S9" "N591_S26" "BF15_S35" "BF7_S13" "N824_S31")
        nre70=("N416_S22" "BF12_S18" "BF4_S10" "N596_S27" "N613_S30" "BF8_S14" "N830_S32")
        nre100=("N419_S23" "BF13_S19" "BF5_S11" "N599_S28" "BF17_S36" "BF9_S15" "N833_S33")
        nre180=("N424_S24" "N425_S25" "BF14_S20" "BF6_S12" "N605_S29" "BF18_S37" "BF10_S16" "N839_S34")
	
	nres=("nre30" "nre70" "nre100" "nre180")

	
        # 1) Generate a depth file from BAM files to calculate abundance
	cd ${assembam_dir}
	bams=(`ls | grep '_sorted_contigs.bam$'`)
	cmd="${jgi_summarize_bam_contig_depths} ${bams[*]} --outputDepth ${assemdepth_dir}/nres_depth.txt --pairedContigs ${assemdepth_dir}/nres_paired.txt --minContigLength 1000 --minContigDepth 2"
	echo $cmd
	eval $cmd

	# 2) Bin contigs
	# When assembly is poor quality or from highly complex community,
	# this one has greatest number of large contigs.
	# It appears that there are significant contaminations in terms of both strain level or above, so it is not advised to use lower minContig cutoff.
	# (reference : https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices)
	for nre in ${nres[@]}; do
	    cmd="${metabat2} -i ${assem_dir}/${nre}/${nre}_assembly.contigs.fa -a ${assemdepth_dir}/nres_depth.txt -o ${bin_dir} -v --sameTNF saved.tnf --saveDistance saved.dist"
#	    cmd="${metabat2} -i ${assem_dir}/${nre}/${nre}_assembly.contigs.fa -a ${assemdepth_dir}/nres_depth.txt -o ${bin_dir} -v --sameTNF saved.tnf --saveDistance saved.dist -B 20"    #Ensemble binning? ????
	    echo $cmd
	    eval $cmd

	done
	
	
