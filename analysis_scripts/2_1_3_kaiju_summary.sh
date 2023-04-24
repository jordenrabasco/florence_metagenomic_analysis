#!/bin/bash

###################################################################################
#	Taxonomic assignment of reads using Kaiju
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"
	clean_dir="${out_root}/clean_fastq"
        kdb_dir="${out_root}/kaiju_db"

	# Output
	pretaxa_dir="${out_root}/pre_taxa"
#	db="refseq"
        db="nr_euk"

	cd ${pretaxa_dir}/${db}
	mkdir -p taxa_summary taxa_plot
	
        # Tool
        kaiju="/home4/sjeong6/tools/kaiju/bin/kaiju"
        kaiju_makedb="/home4/sjeong6/tools/kaiju/bin/kaiju-makedb"
	kaiju2table="/home4/sjeong6/tools/kaiju/bin/kaiju2table"

	cd ${pretaxa_dir}/${db}
	outs=(`ls | grep '_pre_taxa.out'`)

	# A. ALL taxa levels  => Create stacked bars using 2_1_4_taxa_summary.R
	cmd="${kaiju2table} -t ${kdb_dir}/${db}/nodes.dmp -n ${kdb_dir}/${db}/names.dmp -r species ${outs[*]} \
			    -o ${pretaxa_dir}/${db}/taxa_summary/kaiju_${db}_summary.tsv -l superkingdom,phylum,class,order,family,genus,species"
	echo $cmd
	eval $cmd
	
	# B. Specific level  => Create stacked bars using 2_1_5_taxa_summary.R

	levels=("phylum", "class", "order", "family", "genus")   # "species" is the same with the above run for ALL taxa levels

	for level in ${levels[@]}; do
	    cmd="${kaiju2table} -t ${kdb_dir}/${db}/nodes.dmp -n ${kdb_dir}/${db}/names.dmp -r ${level} ${outs[*]} \
	    			-o ${pretaxa_dir}/${db}/taxa_summary/kaiju_${db}_${level}_summary.tsv -l superkingdom,phylum,class,order,family,genus,species"
	    echo $cmd
	    eval $cmd
	done
