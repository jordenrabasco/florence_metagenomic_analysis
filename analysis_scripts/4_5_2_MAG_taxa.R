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
