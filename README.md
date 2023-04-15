# florence_metagenomic_analysis
In the analysis documented below we analyzed shotgun sequencing data sampled from various locations before and after the landfall of hurricane florence. The sequencing data was quality controlled and taxonomically anlyzed to detect taxonomic shifts in the local estuary produced by Florence's landfall and the subsequence ecological perturbation that followed. A detailed summary of the goals of this analysis can be found in the included Plan of Work document. 

This document is intented as an analysis summary and sample workflow from which our procedure can be replicated. The aforementioned procedure is broken up into various steps seperated out in different bash scripts to emphasize each step in the process. The scripts are intended to be run sequentially in the order in which they are presented. 

## Setup

This script sets up the work folders and directories needed in the downstream analysis. The installation of all of the required tools is also included in this section. The tools used in this analysis are included below;

1) Fastqc - Read quality; https://anaconda.org/bioconda/fastqc
2) Trimmomatic - Read trimming; https://anaconda.org/bioconda/trimmomatic
3) Samtools - File conversions and managment; https://anaconda.org/bioconda/samtools
4) Bowtie - Sequence alignment and sequence analysis; https://anaconda.org/bioconda/bowtie
5) Kaiju - taxonomic classification;https://anaconda.org/bioconda/kaiju
6) Megahit - NGS de novo assembly; https://anaconda.org/bioconda/megahit
7) Quasit - Quality Assessment for Genome Assemblies; https://anaconda.org/bioconda/quast
8) Metabat2 - Binning and genome reconstruction; https://anaconda.org/bioconda/metabat2

The tools listed above weree installed locally outside of a conda enviorment and that code is included in the document. For ease of reproduction our recommendation is to install the tools listed via a conda enviorment, instead of what is presented here. 

```shell
# change root directory to the folder which will house all of your other analysis folder
  root="/home4/sjeong6/Paerl" 
# change data directory to the folder which will house your data (you will just need to change "data" in this instance) 
  dat_root="${root}/data"
# change out_root directory to the folder which will house your output files (you will just need to change "output2" in this instance)
  out_root="${root}/output2" 
  
# set up output sub-folders, in the order of output
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
		
# Make all subfolders described above in your local system; these lines are commented out to avoid overwritting folders in subsequent runs of the analysis   
	mkdir -p $out_root && cd $out_root
	mkdir -p $log_dir $logstrim_dir $merge_dir $trim_dir $filter_dir $clean_dir $pool_dir $qc_dir $qc2_dir $kdb_dir $pretaxa_dir $assem_dir $assem2_dir $bin_dir


	################# INSTALLATION #################
	tool_dir="/home4/sjeong6/tools"
	cd ${tool_dir}

	# metabat2
    	conda install -c bioconda metabat2
    	conda install -c bioconda/label/cf201901 metabat2
```

## FastQC,

Initally all reads were processed thorugh fastqc before trimming to assess if the trimming and intial decontamination works and has maintained the integrity of the data.

FastQC Before Trimming
```shell
# Output folder specified
	qc_dir="${out_root}/QC"

# Specify the fastqc location in your system if installed not through a conda enviorment
	fastqc="/opt/fastqc/0.11.9/fastqc"

# enter data folder for a specific data collection and get list of all files in folder
	cd ${dat_root}/NVS139B_fastq
	fqs=(`ls | grep '.fastq.gz'`)

# run fastqc on all files specified above
	for fq in ${fqs[@]}; do
	    cmd="$fastqc $fq -o ${qc_dir}"
	    echo $cmd
	    eval $cmd
	done
```
The fastqc 
Trimming via trimmomatic
```shell
# Output Directories
	trim_dir="${out_root}/trim"
	logstrim_dir="${out_root}/logs_trim"

# Specify the trimmomatic location in your system 
	trimmomatic="/opt/Trimmomatic/0.39/trimmomatic-0.39.jar"

# enter data folder for a specific data collection and get list of all files in folder
	cd ${dat_root}/NVS139B_fastq
	r1s=(`ls | grep '_R1_001.fastq.gz$'`)

	# run fastqc with print out which file is being trimmed
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

Contaminant Analysis
```shell
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
#	wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
#	unzip GRCh38_noalt_as.zip
	
	cd ${trim_dir}
       	r1s=(`ls | grep '_R1_trimmed.fastq.gz$'`)

	#r1="BF10_S16_R1_trimmed.fastq.gz"
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
