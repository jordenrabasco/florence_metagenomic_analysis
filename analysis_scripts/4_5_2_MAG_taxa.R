#!/opt/R/4.1.0/bin/R

###########################################################################################################
#       Stacked bar chart to show taxanomic diversity of bins
###########################################################################################################

	source("/home4/sjeong6/Paerl/scripts2/4_5_1_MAG_taxa.R")
	
	######################################## 
        ##### Stack bar plots
	########################################
        ####################################################
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
	
	z <- lapply(nres_summary, function(x) x %>% group_by_at(all_of(level)) %>% summarize_at(names(.)[!names(.) %in% c("bin", "contigs.tot.len", levels)], sum)) %>%
               reduce(full_join, by = level) %>% as.data.frame
	z[is.na(z)] <- 0
	classified <- z[!grepl("Unclassified|NA",z[,level]),]
	unclassified <- c("Unclassified", z[grep("Unclassified|NA",z[,level]), -1] %>% colSums(.))
	z <- rbind(classified, unclassified)
	z[ names(z)[names(z)!=level]] <- sapply(z[, names(z)[names(z)!=level]], as.numeric)
	z <- z %>% arrange(-rowMeans(.[, names(.) != level]))
	

	# 1) Including "Unclassified"
        z1 <- z[grep("Unclassified",z[,level]),]
        z2 <- z[!grepl("Unclassified",z[,level]),-1] %>% colSums(.)
        z2 <- c("Classified", z2)
        zz <- rbind(z1, z2)
        zz[, names(zz)[names(zz)!=level]] <- sapply(zz[, names(zz)[names(zz)!=level]], as.numeric)

	rel_z <- zz %>% mutate_at(names(.)[names(.) != level], function(x) 100 * x / sum(x)) # The column sum per sample = 100
        colnames(rel_z)[1] <- "taxa"

        col.cnt <- length(rel_z$taxa)
        getPalette <- colorRampPalette(brewer.pal(n=12, "Set3"))
        cols <- getPalette(col.cnt); names(cols) <- rel_z$taxa

        p1 <- rel_z %>% mutate_at("taxa", function(x) factor(x, levels = .[,"taxa"] %>% unlist %>% rev)) %>%
                      pivot_longer(cols = names(.)[names(.) != "taxa"], names_to = "sampleID", values_to = "read") %>%
                      left_join(meta_grp, by = "sampleID") %>% 
                      arrange(taxa, Date, desc(read)) %>% mutate(sampleID = factor(sampleID, levels = .[.[,"taxa"] == "Classified", "sampleID"] %>% unlist %>% as.vector)) %>%
                      ggplot(aes( x = sampleID, y = read, fill = taxa)) + geom_col(position = "stack") +
                                 facet_grid(.~NRE, space = "free", scales = "free") +
                                 ggtitle(paste0("Taxonomic distribution of bins at ", level, " level")) +
                                 xlab("Sample ID") + ylab("Proportion of reads (%)") + theme(legend.position = "right") +
                                 scale_y_continuous(limits = c(-1, 101), expand = c(0,0)) + scale_fill_manual(values = cols)

        # 2) Excluding "Unclassified"
	w <- z[!grepl("Unclassified",z[,level]), ]	
	w1 <- w[1:(num),]
        w2 <- w[(num+1):nrow(w),-1] %>% colSums(.)
        w2 <- c("others", w2)
        ww <- rbind(w1, w2)
        ww[, names(ww)[names(ww)!=level]] <- sapply(ww[, names(ww)[names(ww)!=level]], as.numeric)

        rel_w <- ww %>% mutate_at(names(.)[names(.) != level], function(x) 100 * x / sum(x)) # The column sum per sample = 100
        colnames(rel_w)[1] <- "taxa"

	col.cnt <- length(rel_w$taxa)
	getPalette <- colorRampPalette(brewer.pal(n=12, "Set3"))
        cols <- getPalette(col.cnt); names(cols) <- rel_w$taxa

        p2 <- rel_w %>% mutate_at("taxa", function(x) factor(x, levels = .[,"taxa"] %>% unlist %>% rev)) %>%
                      pivot_longer(cols = names(.)[names(.) != "taxa"], names_to = "sampleID", values_to = "read") %>%
                      left_join(meta_grp, by = "sampleID") %>%
                      arrange(taxa, Date, desc(read)) %>% mutate(sampleID = factor(sampleID, levels = .[.[,"taxa"] == "others", "sampleID"] %>% unlist %>% as.vector)) %>%
                      ggplot(aes( x = sampleID, y = read, fill = taxa)) + geom_col(position = "stack") +
                                 facet_grid(.~NRE, space = "free", scales = "free") +
                                 ggtitle(paste0("Taxonomic distribution of bins at ", level, " level w/o unclassified bins"))    +
                                 xlab("Sample ID") + ylab("Proportion of reads (%)") + theme(legend.position = "right") +
                                 scale_y_continuous(limits = c(-1, 101), expand = c(0,0)) + scale_fill_manual(values = cols)
        
	grid.arrange(p1,p2, nrow = 2)
	
        dev.off()
