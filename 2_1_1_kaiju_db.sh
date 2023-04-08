#!/bin/bash

###################################################################################
#	Build Kaiju DB
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"

        # Outputs
	kdb_dir="${out_root}/kaiju_db"

        # Tool
	kaiju="/home4/sjeong6/tools/kaiju/bin/kaiju"
	kaiju_makedb="/home4/sjeong6/tools/kaiju/bin/kaiju-makedb"

	cd ${kdb_dir}
	# refseq : Completely assembled and annotated reference genomes of Archaea, Bacteria, and viruses from the NCBI RefSeq database
#	${kaiju_makedb} -s refseq
#	mv ./*.dmp ./refseq

	# nr_euk : Subset of NCBI BLAST nr database containing all proteins belonging to Archaea, Bacteria and Viruses + proteins from fungi and microbial eukaryotes
	${kaiju_makedb} -s nr_euk
#	mv ./*.dmp ./nr_euk
