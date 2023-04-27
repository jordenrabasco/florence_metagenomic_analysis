#!/opt/R/4.1.0/bin/R

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

	# !!!!!!!!!!!!! assign the DB between refseq and nr_euk which you want to analyze
        db="refseq"
#	db="nr_euk"
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	pretaxa_dir <- file.path(pretaxa_dir, db) # update pretaxa_dir with db

	taxasum_dir <- file.path(pretaxa_dir, "taxa_summary")
	taxaplot_dir <- file.path(pretaxa_dir, "taxa_plot")
	
	if (!dir.exists(taxasum_dir)) dir.create(taxasum_dir)	
	if (!dir.exists(taxaplot_dir)) dir.create(taxaplot_dir)
	
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
#	level="class" ; num=15
#	level="order" ; num=15
	level="family" ; num=15
#	level="genus" ; num=20
#	level="species" ; num=20
	########################################

	# Export the stacked bar plot
	filename <- paste0(level, "_bar_", db, ".pdf")
        pdf(file.path(taxaplot_dir, filename), width = 15, height = 6)
	op <- par(mfrow=c(2,1), pty="s")

	# 1) Including "unclassified" and "cannot_be_assigned_to_a_(non−viral)_species"
        z <- taxa %>% group_by_at(all_of(level)) %>% summarize_at(names(.)[!names(.) %in% levels], sum) %>%
        arrange(-rowMeans(.[, names(.) != level]))
        z1 <- z[1:(3),]
        z2 <- z[(3+1):nrow(z),-1] %>% colSums(.)
        z2 <- c("others", z2)
        zz <- rbind(z1, z2)
        zz[, names(zz)[names(zz)!=level]] <- sapply(zz[, names(zz)[names(zz)!=level]], as.numeric)

	rel_z <- zz %>% mutate_at(names(.)[names(.) != level], function(x) 100 * x / sum(x)) # The column sum per sample = 100
	colnames(rel_z)[1] <- "taxa"

        col.cnt <- length(rel_z$taxa)
        getPalette <- colorRampPalette(brewer.pal(n=12, "Set3"))
        cols <- getPalette(col.cnt); names(cols) <- rel_z$taxa

	rel_z %>% mutate_at("taxa", function(x) factor(x, levels = .[,"taxa"] %>% unlist %>% rev)) %>%
	  	      pivot_longer(cols = names(.)[names(.) != "taxa"], names_to = "sampleID", values_to = "read") %>%
	      	      left_join(meta_grp, by = "sampleID") %>% 
	      	      arrange(taxa, Date, desc(read)) %>% mutate(sampleID = factor(sampleID, levels = .[.[,"taxa"] == "others", "sampleID"] %>% unlist %>% as.vector)) %>%
	      	      ggplot(aes( x = sampleID, y = read, fill = taxa)) + geom_col(position = "stack") +
	      			 facet_grid(.~NRE, space = "free", scales = "free") +
				 ggtitle(paste0("Taxonomic distribution of reads at ", level, " level")) +
				 xlab("Sample ID") + ylab("Proportion of reads (%)") + theme(legend.position = "right") +
				 scale_y_continuous(limits = c(-1, 101), expand = c(0,0)) + scale_fill_manual(values = cols)
	
	# 2) Excluding "unclassified" and cannot_be_assigned_to_a_(non-viral)_species, NA
        z <- taxa %>% group_by_at(all_of(level)) %>% summarize_at(names(.)[!names(.) %in% levels], sum) %>%
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
              	      pivot_longer(cols = names(.)[names(.) != "taxa"], names_to = "sampleID", values_to = "read") %>%
              	      left_join(meta_grp, by = "sampleID") %>%
              	      arrange(taxa, Date, desc(read)) %>% mutate(sampleID = factor(sampleID, levels = .[.[,"taxa"] == "others", "sampleID"] %>% unlist %>% as.vector)) %>%
              	      ggplot(aes( x = sampleID, y = read, fill = taxa)) + geom_col(position = "stack") +
                                 facet_grid(.~NRE, space = "free", scales = "free") +
				 ggtitle(paste0("Taxonomic distribution of reads at ", level, " level_", "w/o NAs"))	+
                                 xlab("Sample ID") + ylab("Proportion of reads (%)") + theme(legend.position = "right") +
                                 scale_y_continuous(limits = c(-1, 101), expand = c(0,0)) + scale_fill_manual(values = cols)
	
	par(op)
        dev.off()

	# "NA"
	show(taxa_cnts[grep(";NA;",taxa_cnts$taxon_name), ! colnames(taxa_cnts) %in% levels(samp)])
	taxa_cnts[grep(";NA;",taxa_cnts$taxon_name), "taxon_name"] %>% length    