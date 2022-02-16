
# PCA using log-transformed counts ---------------------------------------------------


pca = prcomp(t( log10(cm+1)) )
print(summary(pca))

pcaData = as.data.frame(pca$x)
pcaData$sample = rownames(pcaData)
pcaData = merge(pcaData, cdata,by.x=0,by.y=0)
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))

p_12 = ggplot(data=pcaData, aes(x = PC1, y = PC2, fill=condition)) + 
	geom_point(pch=21,size=8,stroke=1,color="grey20") + 
	theme_bw() + scale_fill_manual(values = mycols ) +
	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
	ggpubr::theme_pubr(base_size = 12, border=T, legend="right")
ggsave(filename = "PC12_wNaive.pdf", p_12, height=3, width=4.5)

p_12_wlab = p_12 + ggrepel::geom_text_repel(label=pcaData$Row.names, size=3)
ggsave(filename = "PC12_wNaive_wlab.pdf", p_12_wlab, height=3, width=4.5 )

p_34 = ggplot(data=pcaData, aes(x = PC3, y = PC4, fill=condition)) + 
	geom_point(pch=21,size=8,stroke=1,color="grey20") + 
	theme_bw() + scale_fill_manual(values = mycols ) +
	xlab(paste0("PC3: ", percentVar[3], "% variance")) +
	ylab(paste0("PC4: ", percentVar[4], "% variance")) + 
	ggpubr::theme_pubr(base_size = 12, border=T, legend="right")
ggsave(filename = "PC34_wNaive.pdf", p_34, height=3, width=4.5 )


p_13 = ggplot(data=pcaData, aes(x = PC1, y = PC3, fill=condition)) + 
	geom_point(pch=21,size=8,stroke=1,color="grey20") + 
	theme_bw() + scale_fill_manual(values = mycols ) +
	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC3: ", percentVar[3], "% variance")) + 
	ggpubr::theme_pubr(base_size = 12, border=T, legend="right")
ggsave(filename = "PC13_wNaive.pdf", p_13, height=3, width=4.5 )

p_13_wlab = p_13 + ggrepel::geom_text_repel(label=pcaData$Row.names, size=3)
ggsave(filename = "PC13_wNaive_wlab.pdf", p_13_wlab, height=3, width=4.5 )


p_23 = ggplot(data=pcaData, aes(x = PC2, y = PC3, fill=condition)) + 
	geom_point(pch=21,size=8,stroke=1,color="grey20") + 
	theme_bw() + scale_fill_manual(values = mycols ) + 
	xlab(paste0("PC2: ", percentVar[2], "% variance")) + 
	ylab(paste0("PC3: ", percentVar[3], "% variance")) + 
	ggpubr::theme_pubr(base_size = 12, border=T, legend="right")
ggsave(filename = "PC23_wNaive.pdf", p_23, height=3, width=4.5 )

## Check rotations & contributions 

varexp = data.frame(x=1:length(percentVar), y=percentVar)
varexp$x = factor(varexp$x)
p_varp = ggplot(data=varexp, aes(x=x, y=y)) + geom_bar(stat="identity") + xlab("Principal Component") + ylab("Proportion of variation (%)")+theme_bw()
ggsave(filename = "PC_VariationProportion.pdf", p_varp, height=3, width=4.5 )

loadings = abs(pca$rotation)
contribution = as.data.frame(sweep(loadings, 2, colSums(loadings), "/"))
contribution = contribution[with(contribution, order(-PC1)),]
head(contribution)



# PCA using rlog counts ---------------------------------------------------

res.pca = PCA(t(rlog_counts[unique(rbindlist(DAPsLS)$peakid),]), scale.unit=TRUE, graph=F)
pca_coords_df <- data.frame(res.pca$ind$coord)

pcaData=merge(pca_coords_df, cdata,by.x=0,by.y=0)

q_12 = ggplot(data=pcaData, aes(x = Dim.1, y = Dim.2, fill=condition)) + 
	geom_point(pch=21,size=8,stroke=1,color="grey20") + 
	theme_bw() + scale_fill_manual(values = mycols ) +
	xlab(paste0("PC1: ", round(res.pca$eig[1,2],2), "% variance")) +
	ylab(paste0("PC2: ", round(res.pca$eig[2,2],2), "% variance")) +
	ggpubr::theme_pubr(base_size = 12, border=T, legend="right")
ggsave(filename = "UseAll3CompDEPs_PC12_wNaive_rlog.pdf", q_12, height=3, width=4.5 )

#		panel.border = element_rect(colour = "black", fill=NA, size=1))

q_12_wlab = q_12 + ggrepel::geom_text_repel(label=pcaData$Row.names, size=3)
ggsave(filename = "UseAll3CompDEPs_PC12_wNaive_rlog_wlab.pdf", q_12_wlab, height=3, width=4.5 )


q_23 = ggplot(data=pcaData, aes(x = Dim.2, y = Dim.3, fill=condition)) + 
	geom_point(pch=21,size=8,stroke=1,color="grey20") + 
	theme_bw() + scale_fill_manual(values = mycols ) +
	xlab(paste0("PC2: ", round(res.pca$eig[2,2],2), "% variance")) +
	ylab(paste0("PC3: ", round(res.pca$eig[3,2],2), "% variance")) +
	ggpubr::theme_pubr(base_size = 12, border=T, legend="right")
ggsave(filename = "UseAll3CompDEPs_PC23_wNaive_rlog.pdf", q_23, height=3, width=4.5 )


q_13 = ggplot(data=pcaData, aes(x = Dim.1, y = Dim.3, fill=condition)) + 
	geom_point(pch=21,size=8,stroke=1,color="grey20") + 
	theme_bw() + scale_fill_manual(values = mycols ) +
	xlab(paste0("PC1: ", round(res.pca$eig[1,2],2), "% variance")) + 
	ylab(paste0("PC3: ", round(res.pca$eig[3,2],2), "% variance")) +
	ggpubr::theme_pubr(base_size = 12, border=T, legend="right")
ggsave(filename = "UseAll3CompDEPs_PC13_wNaive_rlog.pdf", q_13, height=3, width=4.5 )


q_34 = ggplot(data=pcaData, aes(x = Dim.3, y = Dim.4, fill=condition)) + 
	geom_point(pch=21,size=8,stroke=1,color="grey20") + 
	theme_bw() + scale_fill_manual(values = mycols ) +
	xlab(paste0("PC3: ", round(res.pca$eig[3,2],2), "% variance")) + 
	ylab(paste0("PC4: ", round(res.pca$eig[4,2],2), "% variance")) +
	ggpubr::theme_pubr(base_size = 12, border=T, legend="right")
	
ggsave(filename = "UseAll3CompDEPs_PC34_wNaive_rlog.pdf", q_34, height=3, width=4.5 )

