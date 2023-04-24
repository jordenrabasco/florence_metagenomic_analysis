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
	
        # Outputs
        assem_dir="${out_root}/assembly"
	assemqc_dir="${assem_dir}/assembly_QC"
	assemdb_dir="${assem_dir}/assem_db"
	assemmap_dir="${assem_dir}/map_reads2assembly"
	assembam_dir="${assemmap_dir}/bam"
	assemstat_dir="${assemmap_dir}/stat"

	mkdir -p $assemdb_dir $assemmap_dir $assembam_dir $assemstat_dir
	
        # Tool
        bowtie2="/opt/bowtie2/2.4.4/bowtie2"
        bowtie2_build="/opt/bowtie2/2.4.4/bowtie2-build"
        samtools="/opt/samtools/1.13/bin/samtools"

        # For co-assembly of samples within the same location(nre)
        # column -s, -t < ${dat_root}/meta/NGS_WorkOrder_NRE_MetaG_12_2_21_JS_Sample_Names.csv # | less -#2 -N -S
        
        nre30=("N411_S21" "BF11_S17" "BF3_S9" "N591_S26" "BF15_S35" "BF7_S13" "N824_S31")
        nre70=("N416_S22" "BF12_S18" "BF4_S10" "N596_S27" "N613_S30" "BF8_S14" "N830_S32")
        nre100=("N419_S23" "BF13_S19" "BF5_S11" "N599_S28" "BF17_S36" "BF9_S15" "N833_S33")
        nre180=("N424_S24" "N425_S25" "BF14_S20" "BF6_S12" "N605_S29" "BF18_S37" "BF10_S16" "N839_S34")

	nres=("nre30" "nre70" "nre100" "nre180")

        # 1) Build assembly database
	cd ${assem_dir}
	fa="_assembled.contigs.fa"
	#nre="nre30"
	for nre in ${nres[@]}; do
	    cmd="${bowtie2_build} ./${nre}/${nre}${fa} ${assemdb_dir}/${nre}"
	    echo $cmd
	    eval $cmd
	done
	
	# 2) Map input reads to assembly
	cd ${unzip_dir}
	r1="_cont_removed_R1.fastq"
	r2="_cont_removed_R2.fastq"

        # NRE30
        for samp in ${nre30[@]}; do
            cmd="${bowtie2} --sensitive-local -x ${assemdb_dir}/nre30 \
                            -1 ${samp}${r1} -2 ${samp}${r2} --no-unal -p 4 -S ${assemmap_dir}/${samp}.sam 2> ${assemmap_dir}/${samp}_aligned2nre30.log"
            echo $cmd
            eval $cmd
        done
    
        # NRE70
        for samp in ${nre70[@]}; do
            cmd="${bowtie2} --sensitive-local -x ${assemdb_dir}/nre70 \
                            -1 ${samp}${r1} -2 ${samp}${r2} --no-unal -p 4 -S ${assemmap_dir}/${samp}.sam 2> ${assemmap_dir}/${samp}_aligned2nre70.log"
            echo $cmd
            eval $cmd
        done

        # NRE100
        for samp in ${nre100[@]}; do
            cmd="${bowtie2} --sensitive-local -x ${assemdb_dir}/nre100 \
                            -1 ${samp}${r1} -2 ${samp}${r2} --no-unal -p 4 -S ${assemmap_dir}/${samp}.sam 2> ${assemmap_dir}/${samp}_aligned2nre100.log"
            echo $cmd
            eval $cmd
        done

        # NRE180
        for samp in ${nre180[@]}; do
            cmd="${bowtie2} --sensitive-local -x ${assemdb_dir}/nre180 \
                            -1 ${samp}${r1} -2 ${samp}${r2} --no-unal -p 4 -S ${assemmap_dir}/${samp}.sam 2> ${assemmap_dir}/${samp}_aligned2nre180.log"
            echo $cmd
            eval $cmd
        done

	# SAM to BAM index
        cd ${assemmap_dir}
	sams=(`ls | grep '.sam'`)
	
	#sam="BF10_S16.sam"
        for sam in ${sams[@]}; do
	    pre=`echo $sam | sed 's/.sam//'`
	    echo ${sam}
	    echo ${pre}
	    
            # 3) Convert SAM to BAM file
	    cmd1="${samtools} view ${sam} -b -o ${assembam_dir}/${pre}_assembled_contigs.bam"
            echo $cmd1
            eval $cmd1
	    
	    # 4) Sort BAM file
	    cmd2="${samtools} sort -@ 10 ${assembam_dir}/${pre}_assembled_contigs.bam -o ${assembam_dir}/${pre}_sorted_contigs.bam"
	    echo $cmd2
	    eval $cmd2
	    
	    # 5) Index BAM file
	    cmd3="${samtools} index ${assembam_dir}/${pre}_sorted_contigs.bam"
	    echo $cmd3
	    eval $cmd3
        done

	# 6) Generate assembly stats
	cd ${assembam_dir}
	bams=(`ls | grep '_sorted_contigs.bam$'`)
	

	for bam in ${bams[@]}; do
	    pre=`echo $bam | sed 's/_sorted_contigs.bam//'`
	    cmd="${samtools} idxstats ${bam} > ${assemstat_dir}/${pre}_idxstats.txt"
	    echo $cmd
	    eval $cmd
	done

