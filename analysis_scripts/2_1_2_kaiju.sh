#!/bin/bash

###################################################################################
#	Taxonomic assignment of reads using Kaiju
#       DB :  refseq or nr_euk
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"
	clean_dir="${out_root}/clean_fastq"
        kdb_dir="${out_root}/kaiju_db"

	# Output
#	db="refseq"
	db="nr_euk"
	pretaxa_dir="${out_root}/pre_taxa"
	cd ${pretaxa_dir}
	mkdir -p $db
	
        # Tool
        kaiju="/home4/sjeong6/tools/kaiju/bin/kaiju"
        kaiju_makedb="/home4/sjeong6/tools/kaiju/bin/kaiju-makedb"
	
	cd ${clean_dir}
	r1s=(`ls | grep '_cont_removed_R1.fastq.gz$'`)
#	r1s="BF11_S17_cont_removed_R1.fastq.gz"
	
	# run Kaiju
	#r1="BF10_S16_cont_removed_R1.fastq.gz"
	for r1 in ${r1s[@]}; do
	    pre=`echo $r1 | sed 's/_cont_removed_R1.fastq.gz//'`
	    r2=`echo $r1 | sed 's/_R1./_R2./'`
	    echo $pre
	    ls $r1
	    ls $r2
	    
	    cmd="${kaiju} -t ${kdb_dir}/${db}/nodes.dmp -f ${kdb_dir}/${db}/kaiju_db_${db}.fmi -i $r1 -j $r2 -o ${pretaxa_dir}/${db}/${pre}_pre_taxa.out -z 8 -E 1e-05 -v"  
	    echo $cmd
	    eval $cmd
	done

	
