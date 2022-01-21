## Plot PCA and ranks of PCs
## 11/11/2021

##-------------------------------------------------------------------------------------------

# R
library(data.table)
library(dplyr)
library(stringr)

plotpca <- function( dt, name, condition, colorset ) {

  ## Chose either log-tranformed or Z-normalized RC for performing PCA to avoid any outliers driving PCs
	#prcomp_normRC <- prcomp( t( log( dt +1) ))  # using log transformed data could help reduce effects of extreme outliers
	prcomp_normRC <- prcomp( t( dt )) ## should input Z-normalized RC 

	PCs_normRC <- prcomp_normRC$x %>% data.frame
	PCA_summary <- summary(prcomp_normRC)


	## Important variables of each PC 
	pca_latentsem <- prcomp_normRC$rotation

  ## Get significant leading components
	## PC1 
	SigVarPC1_pos <- signif(sort(pca_latentsem[,1],decreasing=TRUE)[1:30],2)
	SigVarPC1_neg <- signif(sort(pca_latentsem[,1],decreasing=FALSE)[1:30],2)
	print(head(SigVarPC1_pos))
  
  ## PC2
	SigVarPC2_pos <- signif(sort(pca_latentsem[,2],decreasing=TRUE)[1:30],2)
	SigVarPC2_neg <- signif(sort(pca_latentsem[,2],decreasing=FALSE)[1:30],2)


	p_PCA_normRC <- ggplot(data = PCs_normRC, aes(x = PC1, y = PC2 , color = condition  ) ) +
    		geom_point( size = 8, alpha=0.99 ) +
	    	geom_hline(yintercept = 0, colour = "gray80", size=0.1) + geom_vline(xintercept = 0, colour = "gray80", size=0.1) +
    		xlab( paste0("PC1: ", round(PCA_summary$importance[2,1],3)*100 , "%") ) + 
    		ylab( paste0("PC2: ", round(PCA_summary$importance[2,2],3)*100 , "%") ) +
	    	ggtitle( paste0( name,"" ) ) +
	        theme_bw() +
	    	theme(plot.title = element_text(size=22, hjust=0.5), 
				axis.text=element_text(size=rel(1.3)), axis.title=element_text(size=18),
				legend.title = element_text(size=8), legend.key.size = unit(0.5, "cm") , 
				legend.position="right", legend.text = element_text(size = 8), 
  				panel.spacing=grid::unit(2, "cm"), 
				panel.grid.major = element_blank(), panel.grid.minor = element_blank()  ) +
  			scale_color_manual(values = colorset ) #+		
#    	geom_text_repel(data = PCs_normRC, aes(label = rownames(PCs_normRC) ), color = "grey50", size = 4 ) 

	ggsave(paste0("plot_PCA_", name, "_PC12_", Sys.Date(), ".pdf"), p_PCA_normRC, height = 3.2, width = 4 )


p_PCA23_normRC <- ggplot(data = PCs_normRC, aes(x = PC2, y = PC3 , color = condition  ) ) +
    		geom_point( size = 8, alpha=0.99 ) +
	    	geom_hline(yintercept = 0, colour = "gray80", size=0.1) + geom_vline(xintercept = 0, colour = "gray80", size=0.1) +
    		xlab( paste0("PC2: ", round(PCA_summary$importance[2,2],3)*100 , "%") ) + 
    		ylab( paste0("PC3: ", round(PCA_summary$importance[2,3],3)*100 , "%") ) +
	    	ggtitle( paste0( name,"" ) ) +
	        theme_bw() +
	    	theme(plot.title = element_text(size=22, hjust=0.5), 
				axis.text=element_text(size=rel(1.3)), axis.title=element_text(size=18),
				legend.title = element_text(size=8), legend.key.size = unit(0.5, "cm") , 
				legend.position="right", legend.text = element_text(size = 8), 
  				panel.spacing=grid::unit(2, "cm"), 
				panel.grid.major = element_blank(), panel.grid.minor = element_blank()  ) +
  			scale_color_manual(values = colorset ) #+		

	ggsave(paste0("plot_PCA_", name, "_PC23_", Sys.Date(), ".pdf"), p_PCA23_normRC, height = 3.2, width = 4 )


  ## Also plot ranks of variation of PCs 
	
	PC_Var <- t( PCA_summary$importance ) %>% data.frame() %>% 
		mutate(label=paste0( round(Proportion.of.Variance,3) * 100,"%")) %>% 
		mutate(PCs=rownames(.), PC_N=gsub("PC","",rownames(.)) ) %>% mutate(PCs=factor(PCs,levels=unique(PCs)), PC_N=factor(PC_N,levels=unique(PC_N)))
		
	p_rankPCs <- ggplot(PC_Var[1:10,], aes( x = PC_N, y = Proportion.of.Variance) ) +
    	geom_bar( stat="identity", fill="grey70", width=0.8 ) + 
    	xlab("PCs") + ylim(c(0,0.78)) + 
	    ggtitle( paste0( name,"" ) ) +
  		geom_text(aes(label=label), vjust=1.2, color="grey20", size=4.5)+
		theme_classic() + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), plot.title=element_text(size=20))
		    	
	ggsave(paste0("plot_PCs_", name, "_rank_", Sys.Date(), ".pdf"), p_rankPCs, height = 4, width = 5.5 )

}


mark <- "H3K27me3"

for (mark in AllMarks){

	dt <- fread(paste0("../NormRC_",mark,"_wID.txt"))
	dt_Z <- t(scale(t(dt[,-c(1:6)]),center=T,scale=T))
	rownames(dt_Z) <- dt$peakN

	condition <- str_split(colnames(dt)[-c(1:6)], "_",simplify=T)[,1]
	colorset <-  list( Naive = "#c7c7c7", Arm = "#80b1d3", Cl13 = "#fb8072" )

	plotpca(dt_Z, name=mark, condition, colorset )

}






