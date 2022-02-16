## Heatmap plot using pheatmap with annotation bars 
## HH 02/2022

##----------------------------------------------------------------------------------------------------------

library(pheatmap)

##----------------------------------------------------------------------------------------------------------

## Add Venn groups into heatmap annotation 

# Clst_annot <- readRDS( paste0("DAPs_All3Comparison_w_KMean5_cluster_annot.rds"))
Clst_annot <- Clst_annot %>% left_join( GeneOL, by=c("peakN"="Gene.Name"))
peak_clusters <- data.frame(Cluster=Clst_annot$p.kmeans.cluster, OL=Clst_annot$OL_n); rownames(peak_clusters) <- Clst_annot$Row.names
peak_clusters <- peak_clusters %>% arrange(Cluster, OL)
peak_clusters$Cluster <- factor(peak_clusters$Cluster)
peak_clusters$OL <- factor(peak_clusters$OL)

## Reorder heatmap orders 
peak_cluster_ord=peak_clusters
peak_cluster_ord$label=rownames(peak_cluster_ord)
peak_cluster_ord=peak_cluster_ord[order(peak_cluster_ord$Cluster),]
peak_cluster_ord$label=NULL


CellTypes <- data.frame( CellTypes=c(rep("Memory", 3), rep("NT", 3), rep("T",3)))
rownames(CellTypes) <- colnames(Scaled_Matrix)

mycol1 <- list( Cluster= brewer.pal(5, "Set1"), 
				OL = brewer.pal(7, "Set3"), 
				CellTypes =  mycols[-1] )
names(mycol1$Cluster) <- c(1:5)
names(mycol1$OL) <- levels(peak_clusters$OL)


pdf( paste0("pheatmap_DAPs_All3Comparison_KMean5_cluster_wVenn_", Sys.Date(), ".pdf"),width=4,height=6)
pheatmap(Scaled_Matrix[rownames(peak_cluster_ord),],
	scale="none", 
	show_rownames = F, 
	annotation_row = peak_clusters, annotation_colors = mycol1,
	annotation_col = CellTypes, 
	treeheight_row=1, cluster_rows = F, cluster_cols = F)
dev.off()


