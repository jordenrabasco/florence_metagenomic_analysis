# Florence_Metagenomic_Analysis
In the analysis documented below we analyzed shotgun sequencing data sampled from various locations before and after the landfall of hurricane florence. The sequencing data was quality controlled and taxonomically analyzed to detect taxonomic shifts, in a local estuary, produced by Florence's landfall and the subsequence ecological perturbation that followed. A detailed summary of the goals of this analysis can be found in the included Plan of Work document. 

This markdown is intented as an analysis summary and sample workflow from which our procedure can be replicated. The aforementioned procedure is broken up into various steps seperated out in different bash scripts to emphasize each step in the process. The scripts are intended to be run sequentially in the order in which they are presented.

![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/project_overveiw.png)

## Setup

All tools utilized and their versions are included below in list format. All programs used are included in the `tools` folder along with their accompanying databses, which can be found with the bulk data on the harddrive provided.

1) Fastqc, 0.11.9 - Read quality
2) Trimmomatic, 0.39 - Read trimming
3) Samtools, 1.13 - File conversions and managment
4) Bowtie, 2.4.4 - Sequence alignment and sequence analysis 
5) Kaiju, 1.9.2 - taxonomic classification
6) Megahit, 1.2.9 - NGS de novo assembly 
7) Quast, 5.2.0 - Quality Assessment for Genome Assemblies
8) Metabat2, 2.15 - Binning and genome reconstruction
9) CheckM, 1.0.7 - Binning Quality Control

This script sets up the work folders and directories needed in the downstream analysis. 

```shell
###################################################################################
#	set up folder structure
###################################################################################

        root="root"
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
	bin_dir="bin"
	binclass_dir="bin_class"
	bintaxa_dir="bin_taxa"
		
	mkdir -p $out_root && cd $out_root
	mkdir -p $log_dir $logstrim_dir $merge_dir $trim_dir $filter_dir $clean_dir $pool_dir $qc_dir $qc2_dir $kdb_dir $pretaxa_dir $assem_dir $assem2_dir $bin_dir $binclass_dir $bintaxa_dir


	################# INSTALLATION #################
	tool_dir="/home4/sjeong6/tools"
	cd ${tool_dir}

	# metabat2
	conda install -c bioconda metabat2
	conda install -c bioconda/label/cf201901 metabat2


	# for checkm
	conda create -n checkm python=2.7
	conda activate checkm
	conda install -c bioconda numpy matplotlib pysam
	conda install -c bioconda hmmer prodigal pplacer pysam


	mkdir -p ${tool_dir}/checkm_DB
	tar -xvzf checkm_data_2015_01_16.tar.gz -C "${tool_dir}/checkm_DB" --strip 1 > /dev/null
	checkm data setRoot "${tool_dir}/checkm_DB"
	export CHECKM_DATA_PATH=/home4/sjeong6/tools/checkm_DB
	
	# GTDB-Tk
        conda create -n gtdbtk-2.1.1 -c conda-forge -c bioconda gtdbtk=2.1.1
	conda activate gtdbtk-2.1.1

	cd $tool_dir
	wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz
	mkdir -p ${tool_dir}/gtdb_tk_DB
	tar -xvzf gtdbtk_r207_v2_data.tar.gz -C "/home4/sjeong6/tools/gtdb_tk_DB" --strip 1 > /dev/null
	rm gtdbtk_r207_v2_data.tar.gz
	conda env config vars set GTDBTK_DATA_PATH="/home4/sjeong6/tools/gtdb_tk_DB"
	
	
```

## Seqeuence QC

The Sequnce QC portion of the workflow utilized fastQC, trimmomatic, samtools, and bowtie2. This section is intended to trim, decontaminate, and seqeunce analyze the reads before subsequent data analysis. Initally all reads were processed thorugh fastqc before trimming to assess if the trimming and intial decontamination maintained the integrity of the data.

FastQC Before Trimming
```shell
###################################################################################
#	fastqc for raw reads
###################################################################################

        output="output2"

        # Inputs
        root="root"
        dat_root="${root}/data"
        out_root="${root}/${output}"

        # Outputs
	qc_dir="${out_root}/QC"

        # Tool
	fastqc="/opt/fastqc/0.11.9/fastqc"

	cd ${dat_root}/NVS139B_fastq
	fqs=(`ls | grep '.fastq.gz'`)

	# run fastqc
	for fq in ${fqs[@]}; do
	    cmd="$fastqc $fq -o ${qc_dir}"
	    echo $cmd
	    eval $cmd
	done
```

An example quailty score graph from the archetype sampel `BF3_S9` is included;

![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/fastqc_example_before_trimming.png)

As you can see there is a decreace in sequence quality torwards the end of the sequence, but overall the quality of the reads were good. 

Trimmomatic was used to trim any residual Illumina specific sequences such as residual primers, seqeuncing bar codes, etc. found in the IlluminaClip sequence library. 

Trimming via trimmomatic
```shell
###################################################################################
#	Trim adapters and filter the raw reads using Trimmomatic
###################################################################################

        output="output2"

        # Inputs
        root="root"
        dat_root="${root}/data"
        out_root="${root}/${output}"

        # Outputs
	trim_dir="${out_root}/trim"
	logstrim_dir="${out_root}/logs_trim"

        # Tool
	trimmomatic="/opt/Trimmomatic/0.39/trimmomatic-0.39.jar"

	cd ${dat_root}/NVS139B_fastq
	r1s=(`ls | grep '_R1_001.fastq.gz$'`)

	# run fastqc
	for r1 in ${r1s[@]}; do
	    pre=`echo $r1 | sed 's/_L001_R1_001.fastq.gz//'`
	    r2=`echo $r1 | sed 's/_R1_/_R2_/'`
	    echo $pre
	    ls $r1
	    ls $r2
	    
	    cmd="java -jar ${trimmomatic} PE -threads 4 -trimlog ${logstrim_dir}/${pre}_trimmed.log ${r1} ${r2} \
	    	           ${trim_dir}/${pre}_R1_trimmed.fastq.gz ${trim_dir}/${pre}_R1_unpaired.fastq.gz \
			   ${trim_dir}/${pre}_R2_trimmed.fastq.gz ${trim_dir}/${pre}_R2_unpaired.fastq.gz \
			   ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50"
	    echo $cmd
	    eval $cmd
	done
```


Contamination reads, sourced from the human genome, were identified and removed utlizing the GRCh38_noalt_as database.

Contaminant Analysis
```shell
###################################################################################
#	Remove host contamination (human) using Bowtie2
#       1) Build host genome NCBI db (GRCh38)
#       2) Map reads against host database
#       3) Convert sam to bam file using samtools
#       4) Extract unmapped reads (non-host)
###################################################################################

        output="output2"

        # Inputs
        root="root"
        dat_root="${root}/data"
        out_root="${root}/${output}"
        trim_dir="${out_root}/trim"
	
        # Outputs
	filter_dir="${out_root}/rm_cont"
	bam_dir="${filter_dir}/bam"
	mkdir -p $bam_dir
	clean_dir="${out_root}/clean_fastq"
	single_dir="${clean_dir}/singleton"
	mkdir -p $single_dir
		
        # Tool
	bowtie2="/opt/bowtie2/2.4.4/bowtie2"
	bowtie2_build="/opt/bowtie2/2.4.4/bowtie2-build"
	samtools="/opt/samtools/1.13/bin/samtools"
	
	# 1) Download GRCh38 db
	cd ${filter_dir}
	wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
	unzip GRCh38_noalt_as.zip
	
	cd ${trim_dir}
       	r1s=(`ls | grep '_R1_trimmed.fastq.gz$'`)

	for r1 in ${r1s[@]}; do
	    pre=`echo $r1 | sed 's/_R1_trimmed.fastq.gz//'`
	    r2=`echo $r1 | sed 's/_R1_/_R2_/'`
	    echo $pre
	    ls $r1
	    ls $r2

	    # 2) Map reads against GRCh38 db  
	    cmd="${bowtie2} --sensitive-local -p 8 -x ${filter_dir}/GRCh38_noalt_as/GRCh38_noalt_as -1 $r1 -2 $r2 -S ${filter_dir}/${pre}_mapped_2human.sam 2> ${filter_dir}/${pre}_mapped_2human.log"
	    echo $cmd
	    eval $cmd
	    
            # 3) Convert sam to bam file using samtools 
	    cmd="${samtools} view -@ 8 -bS ${filter_dir}/${pre}_mapped_2human.sam > ${bam_dir}/${pre}_mapped_2human.bam"
	    echo $cmd
	    eval $cmd
	    
	    # 4) Extract unmapped reads (non-host)
	    cmd="${samtools} view -@ 8 -b -f 12 -F 256 ${bam_dir}/${pre}_mapped_2human.bam > ${bam_dir}/${pre}_unmapped_2human.bam"
	    echo $cmd
	    eval $cmd

	    # 5) Sort BAM file to organize paired reads
	    cmd="${samtools} sort -n ${bam_dir}/${pre}_unmapped_2human.bam -o ${bam_dir}/${pre}_unmapped_sorted.bam"
	    echo $cmd
	    eval $cmd

	    # 6) Convert BAM to fastq
	    cmd="${samtools} fastq -@ 8 ${bam_dir}/${pre}_unmapped_sorted.bam \
	    		     -1 ${clean_dir}/${pre}_cont_removed_R1.fastq.gz -2 ${clean_dir}/${pre}_cont_removed_R2.fastq.gz \
			     -0 ${clean_dir}/${pre}_other.fastq.gz -s ${single_dir}/${pre}_singleton.fastq.gz -n"
	    echo $cmd
	    eval $cmd

	done
```

Post trimming and decontamination another FastQC round was run on each sample in order to assess if there has been any loss of data integrity thorughout the QC process.

FastQC after trimming
```shell
###################################################################################
#	fastqc for raw reads
###################################################################################

        output="output2"

        # Inputs
        root="root"
        dat_root="${root}/data"
        out_root="${root}/${output}"
	clean_dir="${out_root}/clean_fastq"

        # Outputs
	qc2_dir="${out_root}/QC2"

        # Tool
	fastqc="/opt/fastqc/0.11.9/fastqc"

	cd ${clean_dir}
	fqs=(`ls | grep '.fastq.gz'`)
	
	# run fastqc
	for fq in ${fqs[@]}; do
	    cmd="$fastqc $fq -o ${qc2_dir}"
	    echo $cmd
	    eval $cmd
	done
```

An example quailty score graph, the partner of the one included in the first step of the QC process BF3_S9, is included;

![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/fastqc_example_after_trimming.png)

As evident from the graph above the QC process removed the lower quality bases torwards the end of the read. 

## Raw Read Taxonomic Assignment

The read taxonomic assignment utlizes Kaiju to match raw reads with sequences from a predownloaded database either refseq or nr_euk. 
The first step in the workflow is to download and setup the local database. 

Building Kaiju DB
```shell
###################################################################################
#	Build Kaiju DB
###################################################################################

        output="output2"

        # Inputs
        root="root"
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

	 nr_euk : Subset of NCBI BLAST nr database containing all proteins belonging to Archaea, Bacteria and Viruses + proteins from fungi and microbial eukaryotes
	${kaiju_makedb} -s nr_euk
	mv ./*.dmp ./nr_euk
```

After setting up the local databases we then assigned the raw reads a taxonomy database.

Taxonomic Assingment using Kaiju
```shell
###################################################################################
#	Taxonomic assignment of reads using Kaiju
#       DB :  refseq or nr_euk
###################################################################################

        output="output2"

        # Inputs
        root="root"
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
	
	# run Kaiju
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
```

The kaiju assignments were then systematically changed into table format and graphed through an imbeded R script. The Rscript will produce stacked bar charts of the taxonomic assingments of the raw reads.

Kaiju data conversion to table and graph generation 
```shell
###################################################################################
#	Taxonomic assignment of reads using Kaiju
###################################################################################

        output="output2"

        # Inputs
        root="root"
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

```

As we can see there was a large number of unclassified seqeunces in the raw reads. 

![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/Read_tax_unclass.png)

Now lets take a closer look at the taxonomic abundances not including the unclassified reads

![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/Read_tax_no_unclass.png)

The taxonomic abundances are more easily understood without the unclassifed reads interfering with the signal from the data. 

## Assembly

The assembly steps of the workflow involve using megahit to resolve the raw reads into contigs and metaquast to assess contig quality.

A detailed version of the Assembly process can be seen here:
![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/assembly_procedure.png)

Initally the cleaned and trimmed files were unzipping in preparation for the binning procedure. 

Unzip files in preparation for binning
```shell
###################################################################################
#	unzip fastqc.gz files
###################################################################################

        output="output2"

        # Inputs
        root="root"
        dat_root="${root}/data"
        out_root="${root}/${output}"
	clean_dir="${out_root}/clean_fastq"

        # Outputs
	unzip_dir="${clean_dir}/unzip"
	mkdir -p $unzip_dir

	cd ${clean_dir}
	fqs=(`ls | grep '.fastq.gz$'`)

	# unzip files
	for fq in ${fqs[@]}; do
	    cmd="gunzip -k $fq"
	    echo $cmd
	    eval $cmd
	done
```

The binning procedure was then carried out via Megahit to produce the contigs from the raw reads

Assembly using Megahit
```shell
###################################################################################
#	Assembly of reads to contigs using MEGAHIT
###################################################################################

        output="output2"

        # Inputs
        root="root"
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
```

After assembly the contigs were checked for quality, both combined and seperated out by location

Contig Quality Check - all locations together
```shell
###################################################################################
#	QC Assembly 
###################################################################################

        output="output2"

        # Inputs
        root="root"
        dat_root="${root}/data"
        out_root="${root}/${output}"

	clean_dir="${out_root}/clean_fastq"
	unzip_dir="${clean_dir}/unzip"
        assem_dir="${out_root}/assembly"
	
        # Outputs
	assemqc_dir="${assem_dir}/assembly_QC"
	mkdir -p $assemqc_dir
	
        # Tool
	module add python
	quast="/home4/sjeong6/tools/quast-5.2.0/quast.py"

	nre=("nre30" "nre70" "nre100" "nre180")
	
	# Run metaquast for assembly files together
	cd ${assem_dir}
	fa="_assembled.contigs.fa"

	cmd="${quast} -m 50 -t 10 ./${nre[0]}/${nre[0]}${fa} ./${nre[1]}/${nre[1]}${fa} ./${nre[2]}/${nre[2]}${fa} ./${nre[3]}/${nre[3]}${fa} -o $assemqc_dir"
	echo $cmd
	eval $cmd
```

Contig Quality Check - locations seperated
```shell
###################################################################################
#	QC Assembly 
###################################################################################

        output="output2"

        # Inputs
        root="root"
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
```

The metrics produced are best seen from when the locations are seperated as thats how they were assembled. The output metrics can be seen below for further inspection:

![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/assembly_qc_stats.png)

## Metagenomic Binning

The binning steps of the workflow involve using megabat2 resolve contigs into bins, jgi to assess bin depth, and checkm to assess the quality of the resolved bins. After the bins have been generated they were taxonomically assinged using GTDB-Tk.

A detailed version of the Binning process can be seen here:
![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/binning_procedure.png)

The first step of the binning process is to bin the contigs themselves. This was done seperated by location to better assess the differences between before and after rather than between locations. 

Binning Procedure using Metabat2
```shell
###################################################################################
#	Binning for MAGs 
###################################################################################

        output="output2"

        # Inputs
        root="root"
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
	binqc_dir="${bin_dir}/QC"
	mkdir -p $assemdepth_dir $binqc_dir
	
        # Tool
	metabat2="/home4/sjeong6/miniconda3/pkgs/metabat2-2.15-h137b6e9_0/bin/metabat2"
	jgi_summarize_bam_contig_depths="/home4/sjeong6/miniconda3/pkgs/metabat2-2.15-h137b6e9_0/bin/jgi_summarize_bam_contig_depths"
	checkm="/home4/sjeong6/miniconda3/pkgs/checkm-genome-1.0.12-py27_0/bin/checkm"

        nre30=("N411_S21" "BF11_S17" "BF3_S9" "N591_S26" "BF15_S35" "BF7_S13" "N824_S31")
        nre70=("N416_S22" "BF12_S18" "BF4_S10" "N596_S27" "N613_S30" "BF8_S14" "N830_S32")
        nre100=("N419_S23" "BF13_S19" "BF5_S11" "N599_S28" "BF17_S36" "BF9_S15" "N833_S33")
        nre180=("N424_S24" "N425_S25" "BF14_S20" "BF6_S12" "N605_S29" "BF18_S37" "BF10_S16" "N839_S34")
	
	nres=("nre30" "nre70" "nre100" "nre180")

	
        # 1) Generate a depth file from BAM files to calculate abundance
	# For 4 different NRE files, BAM files of 7 or 8 samples within the same NRE location together -> 4 depth.txt files
	cd ${assembam_dir}
	bam="_sorted_contigs.bam"

	# concatenate the sample ID with bam filename(_sorted_contigs.bam)
	nres30=( "${nre30[@]/%/${bam}}" )
	nres70=( "${nre70[@]/%/${bam}}" )
	nres100=( "${nre100[@]/%/${bam}}" )
	nres180=( "${nre180[@]/%/${bam}}" )

	# Run
	# Ref : https://bitbucket.org/berkeleylab/metabat/wiki/Example_Large_Data
	
	cmd="${jgi_summarize_bam_contig_depths} ${nres30[*]} --outputDepth ${assemdepth_dir}/nre30_depth.txt --pairedContigs ${assemdepth_dir}/nre30_paired.txt --minContigLength 1000 --minContigDepth 2"
	echo $cmd
	eval $cmd

	cmd="${jgi_summarize_bam_contig_depths} ${nres70[*]} --outputDepth ${assemdepth_dir}/nre70_depth.txt --pairedContigs ${assemdepth_dir}/nre70_paired.txt --minContigLength 1000 --minContigDepth 2"
        echo $cmd
        eval $cmd

	cmd="${jgi_summarize_bam_contig_depths} ${nres100[*]} --outputDepth ${assemdepth_dir}/nre100_depth.txt --pairedContigs ${assemdepth_dir}/nre100_paired.txt --minContigLength 1000 --minContigDepth 2"
        echo $cmd
        eval $cmd

	cmd="${jgi_summarize_bam_contig_depths} ${nres180[*]} --outputDepth ${assemdepth_dir}/nre180_depth.txt --pairedContigs ${assemdepth_dir}/nre180_paired.txt --minContigLength 1000 --minContigDepth 2"
        echo $cmd
        eval $cmd
	
	# 2) Bin contigs
        # When assembly is poor quality or from highly complex community,
        # this one has greatest number of large contigs.
        # It appears that there are significant contaminations in terms of both strain level or above, so it is not advised to use lower minContig cutoff.
        # (reference : https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices)

	cd ${assem_dir}
	nres=("nre30" "nre70" "nre100" "nre180")
	#nre="nre30"
	for nre in ${nres[@]}; do
	    echo $nre
	    cmd="${metabat2} -i ./${nre}/${nre}_assembled.contigs.fa -a ${assemdepth_dir}/${nre}_depth.txt -o ${bin_dir}/${nre} -v"
	    echo $cmd
	    eval $cmd
	done
```

After the bins were successfully generated checkm was utilized to assess the quality and completness of those bins. 

Quality check of bins with Checkm

```shell
###################################################################################
#	Binn QC using CheckM
##!!!!! Run after "conda activate checkm"
###################################################################################

        output="output2"

        # Inputs
        root="root"
        dat_root="${root}/data"
        out_root="${root}/${output}"

	clean_dir="${out_root}/clean_fastq"
	unzip_dir="${clean_dir}/unzip"
        assem_dir="${out_root}/assembly"
        assemmap_dir="${assem_dir}/map_reads2assembly"
	assembam_dir="${assemmap_dir}/bam"
	
	# Output
	bin_dir="${out_root}/bin/nre100_bins"
        assemdepth_dir="/home4/sjeong6/Paerl/output2/checkm_nre100/depth"
	binqc_dir="/home4/sjeong6/Paerl/output2/checkm_nre100/QC"
	mkdir -p $assemdepth_dir $binqc_dir
	
        # Tool
	checkm data setRoot /home4/sjeong6/tools/checkm_DB
	
	# Bin QC using checkM
	cmd="checkm lineage_wf -x fa ${bin_dir} ${binqc_dir} --tab_table -f ${binqc_dir}/MAGs_checkm.tab --reduced_tree "

	echo $cmd
	eval $cmd
```

Once the bins were assessed for quality they were taxonomically assigned using GTBD-TK, which output taxonomic composition for all of the bins. The bin quality metrics of completeness and contamiantion were resolved into barplots which can be seen here:
![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/bin_completness.png)
![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/bin_contamination.png)

```shell
###################################################################################
#	MAG taxonomic classification using GTDB-Tk
###################################################################################

        output="output2"

        # Inputs
        root="root"
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

	cd ${bin_dir}
	cmd="${gtdbtk} classify_wf --extension fa --genome_dir ${bin_dir} --out_dir ${binclass_dir}"
	echo $cmd
	eval $cmd
```

This data was then resolved into relative abundance graphs using the R code. The section is broken up into two parts; one which maganges and extracts the data and the other which outputs the desired graph. 

Data Managment and Extraction in R
```R
###########################################################################################################
#       Stacked bar chart to show taxanomic diversity of bins
###########################################################################################################

	library(readr)
	library(tidyverse)
	library(ggplot2)
	library(RColorBrewer)
	library(colorRamps)
	library(seqinr)
	library(rlist)
	library(gridExtra)
	
	# I/O
	output <- "output2"
	root <- "root"
	dat_root <- file.path(root, "data")
	out_root <- file.path(root, output)

	assem_dir <- file.path(out_root, "assembly")
	assemdepth_dir <- file.path(assem_dir, "depth")
	bin_dir <- file.path(out_root, "bin")
	binclass_dir <- file.path(out_root, "bin_class")
	binclassify_dir <- file.path(binclass_dir, "classify")
	
	bintaxa_dir <- file.path(out_root, "bin_taxa")
	taxaplot_dir <- file.path(bintaxa_dir, "taxa_plot")
	if (!dir.exists(taxaplot_dir)) dir.create(taxaplot_dir)	

        ##### Sample Info #####
	# Import meta data
	meta <- read.table(file.path(dat_root, "meta/NGS_WorkOrder_NRE_MetaG_12_2_21_JS_Sample_Names.csv"), sep=",", header=T)
	colnames(meta) <- c("sampleID", "ng.uL", "ID", "NRE", "Date")

        meta_grp <- meta %>% select("sampleID", "NRE", "Date")
        meta_grp$NRE <- factor(meta_grp$NRE, levels=c("NRE30", "NRE70", "NRE100", "NRE180"))
        meta_grp$Date <- as.Date(meta_grp$Date, "%m/%d/%Y") %>% factor(., ordered=T)

	##### Contig depths from jgi_summarize_bam_contig_depths
	nre30_depth <- read.delim(file.path(assemdepth_dir, "nre30_depth.txt"))
	nre70_depth <- read.delim(file.path(assemdepth_dir, "nre70_depth.txt"))
	nre100_depth <- read.delim(file.path(assemdepth_dir, "nre100_depth.txt"))
	nre180_depth <- read.delim(file.path(assemdepth_dir, "nre180_depth.txt"))

	depth_red <- function(x) {
		  col1 <- c("contigName", "contigLen", "totalAvgDepth")
		  col2 <- colnames(x)[grepl(".bam$",colnames(x))]
		  x <- x[, c(col1, col2)]

		  col2 <- col2 %>% gsub("_sorted_contigs.bam", "", .) %>% gsub("_S.*$", "", .)
		  colnames(x) <- c(col1, col2)
		  return(x)
	}

	nre30_depth_red <- depth_red(nre30_depth)
	nre70_depth_red <- depth_red(nre70_depth)
	nre100_depth_red <- depth_red(nre100_depth)
	nre180_depth_red <- depth_red(nre180_depth)
	

	##### Bin taxa from GTDB-Tk
	ar53_summary <- read.table(file.path(binclass_dir, "gtdbtk.ar53.summary.tsv"), sep = "\t", header = TRUE, stringsAsFactor = F)   # 12
	bac120_summary <- read.table(file.path(binclass_dir, "gtdbtk.bac120.summary.tsv"), sep = "\t", header = TRUE, stringsAsFactor = F)   # 2453

	bin_taxa <- rbind(bac120_summary, ar53_summary)   # 2465
	
	taxa.2use <- bin_taxa[,c("user_genome","classification", "closest_placement_taxonomy", "other_related_references.genome_id.species_name.radius.ANI.AF.")]
	colnames(taxa.2use) <- c("bin", "classification", "closest_placement_taxonomy", "other_related_ref_species")
	
	levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")
	taxa <- matrix(NA, nrow(taxa.2use), length(levels))

	for (i in 1:nrow(taxa.2use)) {
	    x <- strsplit(taxa.2use[i,"classification"], ";")[[1]]
	    if (length(x)==length(levels)) {
	       taxa[i,] <- gsub("^.__", "", x) %>% gsub(" ", "_", .)
	       taxa[i,][taxa[i,] == ""] <-"NA"
	    } else {
	      taxa[i,] <- rep(gsub(" ", "_", taxa.2use[i,"classification"]), length(levels))
	    }	    
	}
	
	colnames(taxa) <- levels
	taxa.2use <- cbind(taxa.2use, taxa)

	write.table(taxa.2use, file.path(bintaxa_dir, "bin_taxa.csv"), row.names = F, sep = ",", quote = F, col.names=T)

	##### contigs in a bin
	bin_fas <- list.files(bin_dir, pattern=".fa$", full.names=F)

	contigs_in_bin <- NULL
	for (fa in bin_fas) {
	    bin <- sub(".fa$", "", fa)
	    x <- read.fasta(file.path(bin_dir, fa))
	    contigs <- as.array(names(x))
	    df <- data.frame(bin=rep(bin, length(contigs)), contigName=contigs)
	    contigs_in_bin <- rbind(contigs_in_bin, df)
	}

	outfile <- "contigs_in_bin.tab"
	write.table(contigs_in_bin, file.path(bintaxa_dir, outfile), row.names = F, sep = "\t", quote = F, col.names=T)
	
	contigs_in_bin <- read.table(file.path(bintaxa_dir, "contigs_in_bin.tab"), header=T, sep="\t", stringsAsFactors=F, check.names=F)

	nre30_contigs_in_bin <- contigs_in_bin[grep("^nre30", contigs_in_bin$bin),]
	nre70_contigs_in_bin <- contigs_in_bin[grep("^nre70", contigs_in_bin$bin),]
	nre100_contigs_in_bin <- contigs_in_bin[grep("^nre100", contigs_in_bin$bin),]
	nre180_contigs_in_bin <- contigs_in_bin[grep("^nre180", contigs_in_bin$bin),]

	# Merge the contigs_in_bin and depth_red
	nre30_merged <- left_join(nre30_contigs_in_bin, nre30_depth_red, by="contigName")
	nre70_merged <- left_join(nre70_contigs_in_bin, nre70_depth_red, by="contigName")
	nre100_merged <- left_join(nre100_contigs_in_bin, nre100_depth_red, by="contigName")
	nre180_merged <- left_join(nre180_contigs_in_bin, nre180_depth_red, by="contigName")

	bins <- function(x,nre) factor(x$bin, levels=paste0(nre,".",1:length(unique(x$bin))))
	bins30 <- bins(nre30_merged,"nre30")
	bins70 <- bins(nre70_merged,"nre70")
	bins100 <- bins(nre100_merged,"nre100")
	bins180 <- bins(nre180_merged,"nre180")

	nre30_bin_lst <- split(nre30_merged, f=bins30)
	nre70_bin_lst <- split(nre70_merged, f=bins70)
	nre100_bin_lst <- split(nre100_merged, f=bins100)
	nre180_bin_lst <- split(nre180_merged, f=bins180)

	len_weighted_avg <- function(x) {
		cont.len <- as.matrix(t(x[,"contigLen"]))
		cont.depth <- as.matrix(x[, !colnames(x) %in% c("bin", "contigName", "contigLen", "totalAvgDepth") ])
		is.na(cont.depth) <- 0
		avg <- cont.len %*% cont.depth
		tot <- x[,"contigLen"] %>% sum
		result <- data.frame(contigs.tot.len=tot, avg/tot)
		return(result)
	}

	# len_weighted_avg(nre30_bin_lst[[1]])

	nre30_bin_depth <- lapply(nre30_bin_lst, len_weighted_avg) %>% list.rbind %>% rownames_to_column("bin")
	nre70_bin_depth <- lapply(nre70_bin_lst, len_weighted_avg) %>% list.rbind %>% rownames_to_column("bin")
	nre100_bin_depth <- lapply(nre100_bin_lst, len_weighted_avg) %>% list.rbind %>% rownames_to_column("bin")
	nre180_bin_depth <- lapply(nre180_bin_lst, len_weighted_avg) %>% list.rbind %>% rownames_to_column("bin")

	nre30_summary <- merge(taxa.2use[grep("^nre30",taxa.2use$bin), c("bin", levels)], nre30_bin_depth, by="bin")
	nre70_summary <- merge(taxa.2use[grep("^nre70",taxa.2use$bin), c("bin", levels)], nre70_bin_depth, by="bin")
	nre100_summary <- merge(taxa.2use[grep("^nre100",taxa.2use$bin), c("bin", levels)], nre100_bin_depth, by="bin")
	nre180_summary <- merge(taxa.2use[grep("^nre180",taxa.2use$bin), c("bin", levels)], nre180_bin_depth, by="bin")

	nres_summary <- list(nre30=nre30_summary, nre70=nre70_summary, nre100=nre100_summary, nre180=nre180_summary)

	write.table(nre30_summary, file.path(bintaxa_dir, "nre30_bin_taxa_abundance.csv"), row.names = F, sep = ",", quote = F, col.names=T)
	write.table(nre70_summary, file.path(bintaxa_dir, "nre70_bin_taxa_abundance.csv"), row.names = F, sep = ",", quote = F, col.names=T)
	write.table(nre100_summary, file.path(bintaxa_dir, "nre100_bin_taxa_abundance.csv"), row.names = F, sep = ",", quote = F, col.names=T)
	write.table(nre180_summary, file.path(bintaxa_dir, "nre180_bin_taxa_abundance.csv"), row.names = F, sep = ",", quote = F, col.names=T)
```

Graphing Relative Abundances
```R
###########################################################################################################
#       Stacked bar chart to show taxanomic diversity of bins
###########################################################################################################

	source("/tools/scripts2/4_5_1_MAG_taxa.R")
        # set level among ("phylum", "class", "order", "family", "genus", "species")
        # set num as the number of taxa that shown in plots
#	level="domain" ; num=4    X
#       level="phylum" ; num=15
        level="class" ; num=15
#       level="order" ; num=15
#       level="family" ; num=15
#       level="genus" ; num=20
#       level="species" ; num=20

	nres=c("nre30","nre70","nre100","nre180")

	####################################################
	# Export the stacked bar plot
        filename <- paste0("bin_bar_", level, ".pdf")
        pdf(file.path(taxaplot_dir, filename), width = 18, height = 12)
	
        stack_plot <- function(nre) {
            nre_summary	<- nres_summary[[nre]]
	    # 2) Excluding "unclassified" and cannot_be_assigned_to_a_(non-viral)_species, NA
            z <- nre_summary %>% group_by_at(all_of(level)) %>% summarize_at(names(.)[!names(.) %in% c("bin", "contigs.tot.len", levels)], sum) %>%
              arrange(-rowMeans(.[, names(.) != level]))
	     
            z1 <- z[1:(num),]
            z2 <- z[(num+1):nrow(z),-1] %>% colSums(.)
            z2 <- c("others", z2)
            zz <- rbind(z1, z2)
            zz[, names(zz)[names(zz)!=level]] <- sapply(zz[, names(zz)[names(zz)!=level]], as.numeric)

            ww <- zz %>% as.data.frame
            ww <- ww[!grepl('unclassified|cannot|NA', ww[,level]), ]

            rel_w <- ww %>% mutate_at(names(.)[names(.) != level], function(x) 100 * x / sum(x)) # The column sum per sample = 100
            colnames(rel_w)[1] <- "taxa"

	    col.cnt <- length(rel_w$taxa)
	    getPalette <- colorRampPalette(brewer.pal(n=12, "Set3"))
            cols <- getPalette(col.cnt); names(cols) <- rel_w$taxa

            rel_w %>% mutate_at("taxa", function(x) factor(x, levels = .[,"taxa"] %>% unlist %>% rev)) %>%
                      pivot_longer(cols = names(.)[names(.) != "taxa"], names_to = "sampleID", values_to = "bin") %>%
                      left_join(meta_grp[grep(toupper(nre), meta_grp$NRE),c("sampleID","Date")], by = "sampleID") %>%
                      arrange(taxa, Date, desc(bin)) %>% mutate(sampleID = factor(sampleID, levels = .[.[,"taxa"] == "others", "sampleID"] %>% unlist %>% as.vector)) %>%
                      ggplot(aes( x = sampleID, y = bin, fill = taxa)) + geom_col(position = "stack") +
                                 ggtitle(paste0("Taxonomic distribution of bins of ", toupper(nre), " at ", level, " level_", "w/o NAs"))    +
                                 xlab("Sample ID") + ylab("Proportion of bins (%)") + theme(legend.position = "right") +
                                 scale_y_continuous(limits = c(-1, 101), expand = c(0,0)) + scale_fill_manual(values = cols)
        }

        p1 <- stack_plot("nre30");  p2 <- stack_plot("nre70");  p3 <- stack_plot("nre100");   p4 <- stack_plot("nre180")
        grid.arrange(p1,p2,p3,p4, nrow = 2)
	
        dev.off()
```

The relative abundance chart output from the graphing R scripts can be seen here:
![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/bin_taxonomic_abundances.png)

The relative abundances depticated in the graph above are at the class level as it provided the greatest resolution consitency and distiction between samples. For lower levels of classification please see the data located on the harddrive. 

## Wrap Up

Now that you have completed the workflow you should have beatiful graphs displaying the realtive abundances of the bins generated from the Metagenomic experiment of Hurricane Florence's affect on the North Carolinian coast. For questions regarding the workflow or troublehshooting please contact Jorden Rabasco (jrabasc@ncsu.edu) or Sangmi Jeong (sjeong6@ncsu.edu).

