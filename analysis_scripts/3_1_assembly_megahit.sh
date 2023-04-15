#!/bin/bash

###################################################################################
#	Assembly of reads to contigs using MEGAHIT
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"

	clean_dir="${out_root}/clean_fastq"
	unzip_dir="${clean_dir}/unzip"

        # Outputs
	assem_dir="${out_root}/assembly"

        # Tool
	module add python
	megahit="/home4/sjeong6/tools/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit"

	# For co-assembly of samples within the same location(nre)
	# column -s, -t < ${dat_root}/meta/NGS_WorkOrder_NRE_MetaG_12_2_21_JS_Sample_Names.csv # | less -#2 -N -S

	nre30=("N411_S21" "BF11_S17" "BF3_S9" "N591_S26" "BF15_S35" "BF7_S13" "N824_S31")
	nre70=("N416_S22" "BF12_S18" "BF4_S10" "N596_S27" "N613_S30" "BF8_S14" "N830_S32")
	nre100=("N419_S23" "BF13_S19" "BF5_S11" "N599_S28" "BF17_S36" "BF9_S15" "N833_S33")
	nre180=("N424_S24" "N425_S25" "BF14_S20" "BF6_S12" "N605_S29" "BF18_S37" "BF10_S16" "N839_S34")

	nre=("nre30" "nre70" "nre100" "nre180")
	
	# Run MEGAHIT
	kmers="55,75,95"
	
	cd ${unzip_dir}
	r1="_cont_removed_R1.fastq"
	r2="_cont_removed_R2.fastq"

	# NRE30
	cmd="${megahit} --k-list ${kmers} -t 20 -o ${assem_dir}/${nre[0]} --out-prefix ${nre[0]}_assembled \
	     		-1 ${nre30[0]}${r1},${nre30[1]}${r1},${nre30[2]}${r1},${nre30[3]}${r1},${nre30[4]}${r1},${nre30[5]}${r1},${nre30[6]}${r1} \
	     		-2 ${nre30[0]}${r2},${nre30[1]}${r2},${nre30[2]}${r2},${nre30[3]}${r2},${nre30[4]}${r2},${nre30[5]}${r2},${nre30[6]}${r2}"
        echo $cmd
	eval $cmd

        # NRE70
        cmd="${megahit} --k-list ${kmers} -t 20 -o ${assem_dir}/${nre[1]} --out-prefix ${nre[1]}_assembled \
                        -1 ${nre70[0]}${r1},${nre70[1]}${r1},${nre70[2]}${r1},${nre70[3]}${r1},${nre70[4]}${r1},${nre70[5]}${r1},${nre70[6]}${r1} \
                        -2 ${nre70[0]}${r2},${nre70[1]}${r2},${nre70[2]}${r2},${nre70[3]}${r2},${nre70[4]}${r2},${nre70[5]}${r2},${nre70[6]}${r2}"
        echo $cmd
        eval $cmd

	
	# NRE100
	cmd="${megahit} --k-list ${kmers} -t 20 -o ${assem_dir}/${nre[2]} --out-prefix ${nre[2]}_assembled \
			-1 ${nre100[0]}${r1},${nre100[1]}${r1},${nre100[2]}${r1},${nre100[3]}${r1},${nre100[4]}${r1},${nre100[5]}${r1},${nre100[6]}${r1} \
			-2 ${nre100[0]}${r2},${nre100[1]}${r2},${nre100[2]}${r2},${nre100[3]}${r2},${nre100[4]}${r2},${nre100[5]}${r2},${nre100[6]}${r2}"
       echo $cmd
	eval $cmd

	# NRE180
	cmd="${megahit} --k-list ${kmers} -t 20 -o ${assem_dir}/${nre[3]} --out-prefix ${nre[3]}_assembled \
			-1 ${nre180[0]}${r1},${nre180[1]}${r1},${nre180[2]}${r1},${nre180[3]}${r1},${nre180[4]}${r1},${nre180[5]}${r1},${nre180[6]}${r1},${nre180[7]}${r1} \
			-2 ${nre180[0]}${r2},${nre180[1]}${r2},${nre180[2]}${r2},${nre180[3]}${r2},${nre180[4]}${r2},${nre180[5]}${r2},${nre180[6]}${r2},${nre180[7]}${r2}"
	echo $cmd
	eval $cmd
