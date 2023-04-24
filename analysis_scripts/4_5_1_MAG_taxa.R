#!/opt/R/4.1.0/bin/R

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
	root <- "/home4/sjeong6/Paerl"
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

#	write.table(taxa.2use, file.path(bintaxa_dir, "bin_taxa.csv"), row.names = F, sep = ",", quote = F, col.names=T)

	##### contigs in a bin
#	bin_fas <- list.files(bin_dir, pattern=".fa$", full.names=F)

#	contigs_in_bin <- NULL
#	for (fa in bin_fas) {
#	    bin <- sub(".fa$", "", fa)
#	    x <- read.fasta(file.path(bin_dir, fa))
#	    contigs <- as.array(names(x))
#	    df <- data.frame(bin=rep(bin, length(contigs)), contigName=contigs)
#	    contigs_in_bin <- rbind(contigs_in_bin, df)
#	}

#	outfile <- "contigs_in_bin.tab"
#	write.table(contigs_in_bin, file.path(bintaxa_dir, outfile), row.names = F, sep = "\t", quote = F, col.names=T)
	
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

#	write.table(nre30_summary, file.path(bintaxa_dir, "nre30_bin_taxa_abundance.csv"), row.names = F, sep = ",", quote = F, col.names=T)
#	write.table(nre70_summary, file.path(bintaxa_dir, "nre70_bin_taxa_abundance.csv"), row.names = F, sep = ",", quote = F, col.names=T)
#	write.table(nre100_summary, file.path(bintaxa_dir, "nre100_bin_taxa_abundance.csv"), row.names = F, sep = ",", quote = F, col.names=T)
#	write.table(nre180_summary, file.path(bintaxa_dir, "nre180_bin_taxa_abundance.csv"), row.names = F, sep = ",", quote = F, col.names=T)