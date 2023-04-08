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
	ntdb_dir="${dat_root}/nt_db"
	mkdir -p $ntdb_dir

        # Tool
	kraken2="/home4/sjeong6/tools/kraken2/kraken2"
	kraken2_build="/home4/sjeong6/tools/kraken2/kraken2-build"
	
	cd ${ntdb_dir}
	${kraken2_build} --download-library nt --db nt_db
