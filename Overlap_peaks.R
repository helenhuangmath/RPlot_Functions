## Compare and plot Overlap of different set of peaks using ChIPpeakAnno
## HH 01/2022

##-------------------------------------------

library(ChIPpeakAnno) # finding overlaps
library(GenomicRanges)

## Create comparison list, maximum could compare 5 sets 
ls_peak <- list( P1=GRanges(PeakFile1), P2=GRanges(PeakFile2) )

## Use findOverlapsOfPeaks() function to decide overlaps with gap size 
Overlap <- findOverlapsOfPeaks( ls_peak, maxgap = 0  ) 
	
## 	Plot venn diagram of compared peaks
pdf( paste0("p_Venn_ATAC_AmyvsHua_DiffPeak_", compHua, "_", updown, ".pdf"), height = 4, width = 4 )
		makeVennDiagram(Overlap, 
          cex = 1.5, cex.lable = 2, cat.cex =1.5, cex.title = 6 , margin = 0.2,
		 			main = paste0( "Venn Diagram" ) )
dev.off()

## Output Unique and Overlapped peaks 
uniqP <- data.frame(Overlap$uniquePeaks)
olP <- data.frame(Overlap$overlappingPeaks$`P1///P2`)
fwrite(uniqP, paste0("Venn_NotOLPeaks.csv"), sep=",", col.names=T )
fwrite(olP, paste0("Venn_OverlapPeaks.csv"), sep=",", col.names=T )



