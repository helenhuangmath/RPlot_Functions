## Create count and percentage barplot using ggplot 
## HH 01/2022

##--------------------------------------
## Two ways to create barcode:
## 1. Create summary table
## 2. Directly get barplot using ggplot 
##--------------------------------------

## 1st way: 
## Group by samples 
meta_summ <- meta %>% group_by(SampleType, celltype) %>%
	summarise( counts = n() ) %>% mutate( perc = round( counts / sum(counts) * 100, 3) ) %>% data.frame()
meta_summ$celltype <- factor(meta_summ$celltype, levels=celltype_levels)

## Plot proportion of cells per group
p_cN <- ggplot( meta_summ, aes(x = celltype, y = counts, fill = SampleType ) ) +
	geom_bar( position="stack", stat="identity", width=0.8 ) +
  scale_fill_manual( values = mycolors ) +
	xlab("") + ylab("Cells #") + ggtitle("Cell type level1+") +
	theme_bw() + 
	theme( axis.title=element_text(size=16), 
		   axis.text=element_text(size=16), axis.text.x = element_text(size=18, angle=45, hjust=1, vjust=1),
		   plot.title=element_text(size=16),
		   plot.margin = margin(1, 1, 1, 2, "cm") )
ggsave(paste0("p_N_celltype_N.pdf"), p_cN, height = 6, width = 14 )

p_cpct <- ggplot( meta_summ, aes(x = celltype, y=perc, fill = SampleType) ) +
	geom_bar( position="fill", stat="identity" , width=0.8) +
	scale_fill_manual( values = mycolors ) 		 
ggsave(paste0("p_N_celltype_perc.pdf"), p_cpct, height = 6, width = 14 )

##--------------------------------------
## 2nd way: 
meta$celltype <- factor(meta$celltype, levels=celltype_levels)

p_cN <- ggplot(meta, aes(x=celltype, fill=SampleType) ) +
	geom_bar(stat="count", width=0.8 )
  
p_cpct <- ggplot( meta, aes(x=celltype, fill=SampleType) ) +
	geom_bar(stat="count", position="fill", aes(y = (..count..)/sum(..count..)), width=0.8 ) 


##--------------------------------------


