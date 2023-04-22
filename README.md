# Florence_Metagenomic_Analysis
In the analysis documented below we analyzed shotgun sequencing data sampled from various locations before and after the landfall of hurricane florence. The sequencing data was quality controlled and taxonomically anlyzed to detect taxonomic shifts in the local estuary produced by Florence's landfall and the subsequence ecological perturbation that followed. A detailed summary of the goals of this analysis can be found in the included Plan of Work document. 

This document is intented as an analysis summary and sample workflow from which our procedure can be replicated. The aforementioned procedure is broken up into various steps seperated out in different bash scripts to emphasize each step in the process. The scripts are intended to be run sequentially in the order in which they are presented.

![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/project_overveiw.png)

## Setup

This script sets up the work folders and directories needed in the downstream analysis. The installation of all of the required tools is also included in this section. The tools used in this analysis are included below;

1) Fastqc, 0.11.9 - Read quality
2) Trimmomatic, 0.39 - Read trimming
3) Samtools, 1.13 - File conversions and managment
4) Bowtie, 2.4.4 - Sequence alignment and sequence analysis 
5) Kaiju, 1.9.2 - taxonomic classification
6) Megahit, 1.2.9 - NGS de novo assembly 
7) Quast, 5.2.0 - Quality Assessment for Genome Assemblies
8) Metabat2, 2.15 - Binning and genome reconstruction
9) CheckM, 1.0.7 - Binning Quality Control

All programs used are included in the `tools` folder along with their accompanying databsed which is included here in `tar` format with the data. 

```shell
# change root directory to the folder which will house all of your other analysis folder
  root="root_dir" 
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
	tool_dir="tools"
	cd ${tool_dir}
```

## Seqeuence QC

Initally all reads were processed thorugh fastqc before trimming to assess if the trimming and intial decontamination maintained the integrity of the data.

FastQC Before Trimming
```shell
# Output folder specified
	qc_dir="${out_root}/QC"

# Specify the fastqc location in your system if installed not through a conda enviorment
	fastqc="/tools/fastqc/0.11.9/fastqc"

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

An example quailty score graph from the archetype sampel `BF3_S9` is included;

![alt text](https://github.com/jordenrabasco/florence_metagenomic_analysis/blob/main/analysis_support_docs/analysis_images/fastqc_example_before_trimming.png)

As you can see there is a decreace in sequence quality torwards the end of the sequence, but overall the quality of the reads were good. 

Trimmomatic was used to trim any residual Illumina specific sequences such as residual primers, seqeuncing bar codes, etc. found in the IlluminaClip sequence library. 

Trimming via trimmomatic
```shell
# Output Directories
	trim_dir="${out_root}/trim"
	logstrim_dir="${out_root}/logs_trim"

# Specify the trimmomatic location in your system 
	trimmomatic="/tools/Trimmomatic/0.39/trimmomatic-0.39.jar"

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


Contamination reads, sourced from the human genome, were identified and removed utlizing the GRCh38_noalt_as database.

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
	bowtie2="/tools/bowtie2/2.4.4/bowtie2"
	bowtie2_build="/tools/bowtie2/2.4.4/bowtie2-build"
	samtools="/tools/samtools/1.13/bin/samtools"
	
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

        output="output"

        # Inputs
        root="root_dir"
        dat_root="${root}/data"
        out_root="${root}/${output}"
	clean_dir="${out_root}/clean_fastq"

        # Outputs
	qc2_dir="${out_root}/QC2"

        # Tool
	fastqc="/tools/fastqc/0.11.9/fastqc"

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


## Assembly

```shell
       output="output"

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
	${kaiju_makedb} -s refseq
#	mv ./*.dmp ./refseq

	# nr_euk : Subset of NCBI BLAST nr database containing all proteins belonging to Archaea, Bacteria and Viruses + proteins from fungi and microbial eukaryotes
	${kaiju_makedb} -s nr_euk
#	mv ./*.dmp ./nr_euk

```

```shell
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
#	r1s=(`ls | grep '_cont_removed_R1.fastq.gz$'`)
	r1s="BF11_S17_cont_removed_R1.fastq.gz"
	
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


```

```shell
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
#	levels=("phylum")

	for level in ${levels[@]}; do
	    cmd="${kaiju2table} -t ${kdb_dir}/${db}/nodes.dmp -n ${kdb_dir}/${db}/names.dmp -r ${level} ${outs[*]} \
	    			-o ${pretaxa_dir}/${db}/taxa_summary/kaiju_${db}_${level}_summary.tsv -l superkingdom,phylum,class,order,family,genus,species"
	    echo $cmd
	    eval $cmd
	done

```

```shell
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
#	levels=("phylum")

	for level in ${levels[@]}; do
	    cmd="${kaiju2table} -t ${kdb_dir}/${db}/nodes.dmp -n ${kdb_dir}/${db}/names.dmp -r ${level} ${outs[*]} \
	    			-o ${pretaxa_dir}/${db}/taxa_summary/kaiju_${db}_${level}_summary.tsv -l superkingdom,phylum,class,order,family,genus,species"
	    echo $cmd
	    eval $cmd
	done

```

```R
###########################################################################################################
#       Stacked bar chart to show taxanomic diversity of reads
#       In the 2_1_3_kaiju_summary.sh, A. ALL taxa levels 
###########################################################################################################

#        install.packages("devtools")
#        devtools::install_github("jaredhuling/jcolors")
	
	library(readr)
	library(tidyverse)
	library(ggplot2)
#	library(jcolors)
	library(RColorBrewer)
	library(colorRamps)
	
	# I/O
	output <- "output2"
	root <- "/home4/sjeong6/Paerl"
	dat_root <- file.path(root, "data")
	out_root <- file.path(root, output)
	pretaxa_dir <- file.path(out_root, "pre_taxa")

	# !!! assign the DB you want to analyze
#       db="refseq"
	db="nr_euk"
	pretaxa_dir <- file.path(pretaxa_dir, db) # update pretaxa_dir with db

	taxasum_dir <- file.path(pretaxa_dir, "taxa_summary")
	taxaplot_dir <- file.path(pretaxa_dir, "taxa_plot")

        ##### Sample Info #####
	# Import meta data
	meta <- read.table(file.path(dat_root, "meta/NGS_WorkOrder_NRE_MetaG_12_2_21_JS_Sample_Names.csv"), sep=",", header=T)
	colnames(meta) <- c("sampleID", "ng.uL", "ID", "NRE", "Date")

        meta_grp <- meta %>% select("sampleID", "NRE", "Date")
        meta_grp$NRE <- factor(meta_grp$NRE, levels=c("NRE30", "NRE70", "NRE100", "NRE180"))
        meta_grp$Date <- as.Date(meta_grp$Date, "%m/%d/%Y") %>% factor(., ordered=T)

	##### Taxonomic assignment from Kaiju #####
	# Import the kaiju summary table
	taxa_tab <- readr::read_tsv(file.path(taxasum_dir, paste0("kaiju_",db,"_summary.tsv"))) %>% as.data.frame

	x <- taxa_tab
	x$file <- gsub("_S.*$", "", x$file)
	colnames(x) <- c("sampleID", "percent", "read_counts", "taxon_id", "taxon_name")

	samp <- factor(x$sampleID)
	y <- split(x, samp)

	names(y) <- levels(samp)
	y <- lapply(y, function(z) z %>% select(read_counts, taxon_id, taxon_name))

	taxa_cnts <- merge(y[[1]], y[[2]], by=c("taxon_id","taxon_name"), all=T)
        colnames(taxa_cnts)[(ncol(taxa_cnts)-1):ncol(taxa_cnts)] <- levels(samp)[1:2]
        for (i in 3:length(y)) {
            colnames(y[[i]])[1] <- levels(samp)[i]
            taxa_cnts <- merge(taxa_cnts, y[[i]], by=c("taxon_id","taxon_name"), all=T)
        }

	levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")
	taxa <- matrix(NA, nrow(taxa_cnts), length(levels))

	for (i in 1:nrow(taxa_cnts)) {
	    x <- strsplit(taxa_cnts[i,"taxon_name"], ";")[[1]]
	    if (length(x)==length(levels)) {
	       taxa[i,] <- gsub(" ", "_", x)
	    } else {
	      taxa[i,] <- rep(gsub(" ", "_", taxa_cnts[i,"taxon_name"]), length(levels))
	    }	    
	}

	colnames(taxa) <- levels
	z <- taxa_cnts[,! colnames(taxa_cnts) %in% c("taxon_id", "taxon_name")]
	z[is.na(z)] <- 0
	taxa <- cbind(taxa, z)

	filename <- paste0("taxa_reads_",db,".csv")
	write.table(taxa, file.path(pretaxa_dir, filename), row.names = F, sep = ",", quote = F, col.names=T)

	######################################## 
        ##### Stack bar plots
	########################################
	# set level among ("domain", "phylum", "class", "order", "family", "genus", "species")
	# set num as the number of taxa that shown in plots
#	level="domain" ; num=5   # max=5
#	level="phylum" ; num=15
	level="class" ; num=15
#	level="order" ; num=15
#	level="family" ; num=15
#	level="genus" ; num=20
#	level="species" ; num=20

	# Export the stacked bar plot
	filename <- paste0(level, "_bar_", db, ".pdf")
        pdf(file.path(taxaplot_dir, filename), width = 15, height = 6)
	op <- par(mfrow=c(2,1), pty="s")

        z <- taxa %>% group_by_at(all_of(level)) %>% summarize_at(names(.)[!names(.) %in% levels], sum) %>%
        arrange(-rowMeans(.[, names(.) != level]))
        z1 <- z[1:(num),]
        z2 <- z[(num+1):nrow(z),-1] %>% colSums(.)
        z2 <- c("others", z2)
        zz <- rbind(z1, z2)
        zz[, names(zz)[names(zz)!=level]] <- sapply(zz[, names(zz)[names(zz)!=level]], as.numeric)

	# 1) Including "unclassified" and "cannot_be_assigned_to_a_(non−viral)_species"
	rel_z <- zz %>% mutate_at(names(.)[names(.) != level], function(x) 100 * x / sum(x)) # The column sum per sample = 100
	colnames(rel_z)[1] <- "taxa"
		
	rel_z %>% mutate_at("taxa", function(x) factor(x, levels = .[,"taxa"] %>% unlist %>% rev)) %>%
	  	      pivot_longer(cols = names(.)[names(.) != "taxa"], names_to = "sampleID", values_to = "read") %>%
	      	      left_join(meta_grp, by = "sampleID") %>% 
	      	      arrange(taxa, Date, desc(read)) %>% mutate(sampleID = factor(sampleID, levels = .[.[,"taxa"] == "others", "sampleID"] %>% unlist %>% as.vector)) %>%
	      	      ggplot(aes( x = sampleID, y = read, fill = taxa)) + geom_col(position = "stack") +
	      			 facet_grid(.~NRE, space = "free", scales = "free") +
				 ggtitle(paste0("Taxonomic distribution of reads at ", level, " level")) +
				 xlab("Sample ID") + ylab("Proportion of reads (%)") + theme(legend.position = "right") +
				 scale_y_continuous(limits = c(-1, 101), expand = c(0,0)) #+ theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
	
	# 2) Excluding "unclassified" and cannot_be_assigned_to_a_(non-viral)_species, NA
	ww <- zz %>% as.data.frame
	ww <- ww[!grepl('unclassified|cannot|NA', ww[,level]), ]

	rel_w <- ww %>% mutate_at(names(.)[names(.) != level], function(x) 100 * x / sum(x)) # The column sum per sample = 100
        colnames(rel_w)[1] <- "taxa"
		
        rel_w %>% mutate_at("taxa", function(x) factor(x, levels = .[,"taxa"] %>% unlist %>% rev)) %>%
              	      pivot_longer(cols = names(.)[names(.) != "taxa"], names_to = "sampleID", values_to = "read") %>%
              	      left_join(meta_grp, by = "sampleID") %>%
              	      arrange(taxa, Date, desc(read)) %>% mutate(sampleID = factor(sampleID, levels = .[.[,"taxa"] == "others", "sampleID"] %>% unlist %>% as.vector)) %>%
              	      ggplot(aes( x = sampleID, y = read, fill = taxa)) + geom_col(position = "stack") +
                                 facet_grid(.~NRE, space = "free", scales = "free") +
				 ggtitle(paste0("Taxonomic distribution of reads at ", level, " level_", "w/o NAs"))	+
                                 xlab("Sample ID") + ylab("Proportion of reads (%)") + theme(legend.position = "right") +
                                 scale_y_continuous(limits = c(-1, 101), expand = c(0,0)) #+ theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
	
	par(op)
        dev.off()

	# "NA"
	show(taxa_cnts[grep(";NA;",taxa_cnts$taxon_name), ! colnames(taxa_cnts) %in% levels(samp)])
	taxa_cnts[grep(";NA;",taxa_cnts$taxon_name), "taxon_name"] %>% length    
```

```R
###########################################################################################################
#       Stacked bar chart to show taxanomic diversity of reads
###########################################################################################################

#        install.packages("devtools")
#        devtools::install_github("jaredhuling/jcolors")
	
	library(readr)
	library(tidyverse)
	library(ggplot2)
#	library(jcolors)
	library(RColorBrewer)
	library(colorRamps)
	
	# I/O
	output <- "output2"
	root <- "/home4/sjeong6/Paerl"
	dat_root <- file.path(root, "data")
	out_root <- file.path(root, output)
	pretaxa_dir <- file.path(out_root, "pre_taxa")

	# !!!!!!!!!!!!!!!! assign the DB you want to analyze
#       db="refseq"
	db="nr_euk"
	
        # set level among ("domain", "phylum", "class", "order", "family", "genus", "species")
        # set num as the number of taxa that shown in plots
#       level="domain" ; num=5   # max=5
#       level="phylum" ; num=15
        level="class" ; num=15
#       level="order" ; num=15
#       level="family" ; num=15
#       level="genus" ; num=20
#       level="species" ; num=20
	# !!!!!!!!!!!!!!!!
	
	pretaxa_dir <- file.path(pretaxa_dir, db) # update pretaxa_dir with db
        taxasum_dir <- file.path(pretaxa_dir, "taxa_summary")
        taxaplot_dir <- file.path(pretaxa_dir, "taxa_plot")

        ##### Sample Info #####
	# Import meta data
	meta <- read.table(file.path(dat_root, "meta/NGS_WorkOrder_NRE_MetaG_12_2_21_JS_Sample_Names.csv"), sep=",", header=T)
	colnames(meta) <- c("sampleID", "ng.uL", "ID", "NRE", "Date")

        meta_grp <- meta %>% select("sampleID", "NRE", "Date")
        meta_grp$NRE <- factor(meta_grp$NRE, levels=c("NRE30", "NRE70", "NRE100", "NRE180"))
        meta_grp$Date <- as.Date(meta_grp$Date, "%m/%d/%Y") %>% factor(., ordered=T)

	##### Taxonomic assignment from Kaiju #####
	# Import the kaiju summary table
	taxa_tab <- readr::read_tsv(file.path(taxasum_dir, paste0("kaiju_",db,"_",level,"_summary.tsv"))) %>% as.data.frame

	x <- taxa_tab
	x$file <- gsub("_S.*$", "", x$file)
	colnames(x) <- c("sampleID", "percent", "read_counts", "taxon_id", "taxon_name")

	samp <- factor(x$sampleID)
	y <- split(x, samp)

	names(y) <- levels(samp)
	y <- lapply(y, function(z) z %>% select(read_counts, taxon_id, taxon_name))

	taxa_cnts <- merge(y[[1]], y[[2]], by=c("taxon_id","taxon_name"), all=T)
        colnames(taxa_cnts)[(ncol(taxa_cnts)-1):ncol(taxa_cnts)] <- levels(samp)[1:2]
        for (i in 3:length(y)) {
            colnames(y[[i]])[1] <- levels(samp)[i]
            taxa_cnts <- merge(taxa_cnts, y[[i]], by=c("taxon_id","taxon_name"), all=T)
        }

	levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")
	taxa <- matrix(NA, nrow(taxa_cnts), length(levels))

	for (i in 1:nrow(taxa_cnts)) {
	    x <- strsplit(taxa_cnts[i,"taxon_name"], ";")[[1]]
	    if (length(x)==length(levels)) {
	       taxa[i,] <- gsub(" ", "_", x)
	    } else {
	      taxa[i,] <- rep(gsub(" ", "_", taxa_cnts[i,"taxon_name"]), length(levels))
	    }	    
	}

	colnames(taxa) <- levels
	z <- taxa_cnts[,! colnames(taxa_cnts) %in% c("taxon_id", "taxon_name")]
	z[is.na(z)] <- 0
	taxa <- cbind(taxa, z)

#	filename <- paste0("separate_", level, "_taxa_reads_",db,".csv")
#	write.table(taxa, file.path(pretaxa_dir, filename), row.names = F, sep = ",", quote = F, col.names=T)

	######################################## 
        ##### Stack bar plots
	########################################
	# Export the stacked bar plot
        filename <- paste0("tab_separate_", level, "_bar_", db, ".pdf")
        pdf(file.path(taxaplot_dir, filename), width = 15, height = 6)
        op <- par(mfrow=c(2,1), pty="s")

        z <- taxa %>% group_by_at(all_of(level)) %>% summarize_at(names(.)[!names(.) %in% levels], sum) %>%
             arrange(-rowMeans(.[, names(.) != level]))

        z1 <- z[1:(num),]
        z2 <- z[(num+1):nrow(z),-1] %>% colSums(.)
        z2 <- c("others", z2)
        zz <- rbind(z1, z2)
        zz[, names(zz)[names(zz)!=level]] <- sapply(zz[, names(zz)[names(zz)!=level]], as.numeric)

	# 1) Including "unclassified" and "cannot_be_assigned_to_a_(non−viral)_species"
	rel_z <- zz %>% mutate_at(names(.)[names(.) != level], function(x) 100 * x / sum(x)) # The column sum per sample = 100
	colnames(rel_z)[1] <- "taxa"

	rel_z %>% mutate_at("taxa", function(x) factor(x, levels = .[,"taxa"] %>% unlist %>% rev)) %>%
	      pivot_longer(cols = names(.)[names(.) != "taxa"], names_to = "sampleID", values_to = "read") %>%
	      left_join(meta_grp, by = "sampleID") %>% 
	      arrange(taxa, Date, desc(read)) %>% mutate(sampleID = factor(sampleID, levels = .[.[,"taxa"] == "others", "sampleID"] %>% unlist %>% as.vector)) %>%
	      ggplot(aes( x = sampleID, y = read, fill = taxa)) + geom_col(position = "stack") +
	      			 facet_grid(.~NRE, space = "free", scales = "free") +
				 ggtitle(paste0("Taxonomic distribution of reads at ", level, " level")) +
				 xlab("Sample ID") + ylab("Proportion of reads (%)") + theme(legend.position = "right") +
				 scale_y_continuous(limits = c(-1, 101), expand = c(0,0)) #+ theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
		
	# 2) Excluding "unclassified" and "cannot_be_assigned_to_a_(non−viral)_species"
        ww <- zz %>% as.data.frame
        ww <- ww[!grepl('unclassified|cannot|NA', ww[,level]), ]
	
	rel_w <- ww %>% mutate_at(names(.)[names(.) != level], function(x) 100 * x / sum(x)) # The column sum per sample = 100
        colnames(rel_w)[1] <- "taxa"

        rel_w %>% mutate_at("taxa", function(x) factor(x, levels = .[,"taxa"] %>% unlist %>% rev)) %>%
              pivot_longer(cols = names(.)[names(.) != "taxa"], names_to = "sampleID", values_to = "read") %>%
              left_join(meta_grp, by = "sampleID") %>%
              arrange(taxa, Date, desc(read)) %>% mutate(sampleID = factor(sampleID, levels = .[.[,"taxa"] == "others", "sampleID"] %>% unlist %>% as.vector)) %>%
              ggplot(aes( x = sampleID, y = read, fill = taxa)) + geom_col(position = "stack") +
                                 facet_grid(.~NRE, space = "free", scales = "free") +
				 ggtitle(paste0("Taxonomic distribution of reads at ", level, " level", "w/o NAs"))	+
                                 xlab("Sample ID") + ylab("Proportion of reads (%)") + theme(legend.position = "right") +
                                 scale_y_continuous(limits = c(-1, 101), expand = c(0,0)) #+ theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        par(op)
	dev.off()


	# "NA"
	show(taxa_cnts[grep(";NA;",taxa_cnts$taxon_name), ! colnames(taxa_cnts) %in% levels(samp)])
	taxa_cnts[grep(";NA;",taxa_cnts$taxon_name), "taxon_name"] %>% length    

```

```shell
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

```

```shell
###################################################################################
#	unzip fastqc.gz files
###################################################################################

        output="output2"

        # Inputs
        root="/home4/sjeong6/Paerl"
        dat_root="${root}/data"
        out_root="${root}/${output}"
	clean_dir="${out_root}/clean_fastq"

        # Outputs
	unzip_dir="${clean_dir}/unzip"
	mkdir -p $unzip_dir

	cd ${clean_dir}
	fqs=(`ls | grep '.fastq.gz$'`)

	# unzip files
	#fq="BF10_S16_cont_removed_R1.fastq.gz"
	for fq in ${fqs[@]}; do
	    cmd="gunzip -k $fq"
	    echo $cmd
	    eval $cmd
	done

```

```shell
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
```

```shell
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
        assem_dir="${out_root}/assembly"
	
        # Outputs
	assemqc_dir="${assem_dir}/assembly_QC"
	mkdir -p $assemqc_dir
	
        # Tool
	module add python
	quast="/home4/sjeong6/tools/quast-5.2.0/quast.py"
#	quast="/home4/sjeong6/miniconda3/pkgs/quast-5.0.2-py27pl5262h8eb80aa_5/bin/metaquast.py"

	nre=("nre30" "nre70" "nre100" "nre180")
	
	# Run metaquast for assembly files together
	cd ${assem_dir}
	fa="_assembled.contigs.fa"

	cmd="${quast} -m 50 -t 10 ./${nre[0]}/${nre[0]}${fa} ./${nre[1]}/${nre[1]}${fa} ./${nre[2]}/${nre[2]}${fa} ./${nre[3]}/${nre[3]}${fa} -o $assemqc_dir"
	echo $cmd
	eval $cmd
```

```shell
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
#	    echo $cmd
#	    eval $cmd
	done
	
	# 2) Map input reads to assembly
	cd ${unzip_dir}
	r1="_cont_removed_R1.fastq"
	r2="_cont_removed_R2.fastq"

        # NRE30
        for samp in ${nre30[@]}; do
            cmd="${bowtie2} --sensitive-local -x ${assemdb_dir}/nre30 \
                            -1 ${samp}${r1} -2 ${samp}${r2} --no-unal -p 4 -S ${assemmap_dir}/${samp}.sam 2> ${assemmap_dir}/${samp}_aligned2nre30.log"
#            echo $cmd
#            eval $cmd
        done
    
        # NRE70
        for samp in ${nre70[@]}; do
            cmd="${bowtie2} --sensitive-local -x ${assemdb_dir}/nre70 \
                            -1 ${samp}${r1} -2 ${samp}${r2} --no-unal -p 4 -S ${assemmap_dir}/${samp}.sam 2> ${assemmap_dir}/${samp}_aligned2nre70.log"
#            echo $cmd
#            eval $cmd
        done

        # NRE100
        for samp in ${nre100[@]}; do
            cmd="${bowtie2} --sensitive-local -x ${assemdb_dir}/nre100 \
                            -1 ${samp}${r1} -2 ${samp}${r2} --no-unal -p 4 -S ${assemmap_dir}/${samp}.sam 2> ${assemmap_dir}/${samp}_aligned2nre100.log"
#            echo $cmd
#            eval $cmd
        done

        # NRE180
        for samp in ${nre180[@]}; do
            cmd="${bowtie2} --sensitive-local -x ${assemdb_dir}/nre180 \
                            -1 ${samp}${r1} -2 ${samp}${r2} --no-unal -p 4 -S ${assemmap_dir}/${samp}.sam 2> ${assemmap_dir}/${samp}_aligned2nre180.log"
#            echo $cmd
#            eval $cmd
        done

	# SAM to BAM index
        cd ${assemmap_dir}
	sams=(`ls | grep '.sam'`)
	
	#sam="BF10_S16.sam"
        for sam in ${sams[@]}; do
	    pre=`echo $sam | sed 's/.sam//'`
#	    echo ${sam}
#	    echo ${pre}
	    
            # 3) Convert SAM to BAM file
	    cmd1="${samtools} view ${sam} -b -o ${assembam_dir}/${pre}_assembled_contigs.bam"
#            echo $cmd1
#            eval $cmd1
	    
	    # 4) Sort BAM file
	    cmd2="${samtools} sort -@ 10 ${assembam_dir}/${pre}_assembled_contigs.bam -o ${assembam_dir}/${pre}_sorted_contigs.bam"
#	    echo $cmd2
#	    eval $cmd2
	    
	    # 5) Index BAM file
	    cmd3="${samtools} index ${assembam_dir}/${pre}_sorted_contigs.bam"
#	    echo $cmd3
#	    eval $cmd3
        done

	# 6) Generate assembly stats
	cd ${assembam_dir}
	bams=(`ls | grep '_sorted_contigs.bam$'`)
	

	for bam in ${bams[@]}; do
	    pre=`echo $bam | sed 's/_sorted_contigs.bam//'`
	    cmd="${samtools} idxstats ${bam} > ${assemstat_dir}/${pre}_idxstats.txt"
#	    echo $cmd
#	    eval $cmd
	done
```

```shell
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
	assemmap2_dir="${assem_dir}/map_reads2assembly2"
	assembam_dir="${assemmap2_dir}/bam"
	assemstat_dir="${assemmap2_dir}/stat"

	mkdir -p $assemdb_dir $assemmap2_dir $assembam_dir $assemstat_dir
	
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

	cmd="${bowtie2_build} ${assem_dir}/nre30/nre30${fa},${assem_dir}/nre70/nre70${fa},${assem_dir}/nre100/nre100${fa},${assem_dir}/nre180/nre180${fa} ${assemdb_dir}/nres"
#	echo $cmd
#	eval $cmd
	
	# 2) Map input reads to assembly
	cd ${unzip_dir}
	r1s=(`ls | grep '_cont_removed_R1.fastq'`)

	for r1 in ${r1s[@]}; do
	    pre=`echo $r1 | sed 's/_cont_removed_R1.fastq//'`
	    r2=`echo $r1 | sed 's/_R1./_R2./'`

	    cmd="${bowtie2} --sensitive-local -x ${assemdb_dir}/nres -1 ${r1} -2 ${r2} --no-unal -p 4 -S ${assemmap2_dir}/${pre}_aligned2assembly.sam 2> ${assemmap2_dir}/${pre}_aligned2assemly.log"
	    echo $cmd
	    eval $cmd
	done

	# SAM to BAM index
        cd ${assemmap2_dir}
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

	    # 6) Generate assembly stats
	    cmd4="${samtools} idxstats ${assembam_dir}/${pre}_sorted_contigs.bam > ${assemstat_dir}/${pre}_idxstats.txt"
	    echo $cmd4
	    eval $cmd4
        done
```

```shell
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
```

