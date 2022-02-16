## Venn diagram 
## HH 11/2021 

##---------------------------------------------------------------------------------------------------------

library(data.table)
library(dplyr)

library(VennDiagram)
library(RColorBrewer)

library(ggVennDiagram)
library(ggplot2)
library(eulerr)

##-----------------------------------------------------------------

## Prepare input data 

merge <- fread("Genelist.txt", header=T)

list1 <- merge %>% filter( CellType == "JC_Int" )
list2 <- merge %>% filter( CellType == "JC_Progenitor1" )
list3 <- merge %>% filter( CellType == "JC_Progenitor2" )
list4 <- merge %>% filter( CellType == "JC_Term" )

mycol1 <- brewer.pal(4, "Pastel2")
mycol1 <- c( "#fddbc7","#fc8d59","#3288bd","#d1e5f0")

LS_Comps <- list( Cond1 = list1$Gene.Name, Cond2=list2$Gene.Name, Cond3=list3$Gene.Name, Cond4=list4$Gene.Name  )
LS_Names <- c("Cond1", "Cond2", "Cond3", "Cond4")

##-----------------------------------------------------------------

## Calculate overlap

OL_Comps <- calculate.overlap(LS_Comps)
OL_Comps
    
## Save OL results 
GeneOL <- lapply(OL_Comps, data.frame) %>% rbindlist( idcol=T) %>% unique()
colnames(GeneOL) <- c("OL_Groups", "Gene.Name")
    
OL_Summ <- data.frame(table(GeneOL$OL_Groups))
    
GeneOL <- GeneOL %>% inner_join(OL_Summ, by=c("OL_Groups"="Var1")) %>% arrange(OL_Groups, Gene.Name)
fwrite(GeneOL, paste0("OL_Venn_Groups.txt"), sep="\t" )


##---------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------

## 1. use VennDiagram::venn.diagram
## General used way of plotting venn diagram 

p_venn1 <- venn.diagram(
		x = LS_Comps,
		category.names = LS_Names,
		#filename = paste0("p_Venn_venn.diagram_comps.png"), output=F,
		#imagetype="png" , height = 1200 ,  width = 1200 ,  resolution = 300, 
		filename = NULL, 
		main = paste0(" Venn "), 
		lwd = 3, lty = 1, col = mycol1, 
		cex = 1.5, fontfamily = "sans", cat.fontfamily = "sans", main.fontfamily = "sans", 
		cat.default.pos="outer", # cat.pos=c(-20, 120, 80, 60), 
		main.cex = 1.3, cat.cex=1.4, margin = 0.15, cat.dist=0.1 )

pdf(paste0("p_Venn_venn.diagram_comps.pdf"), height=3.2, width=5 )
grid::grid.draw(p_venn1)
dev.off()


##-----------------------------------------------------------------

## 2. use ggVennDiagram::ggVennDiagram
## Similar as venn.diagram 

p_venn2 <- ggVennDiagram(LS_Comps, lwd = 0.2 , set_size=5, label_size=5, edge_size = 1, 
		label="count", label_txtWidth = 60,  label_alpha = 0, 
		category.names = LS_Names, 
		title = paste0(" Venn ") ) +
		scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
		scale_color_manual(values = mycol1  ) +
		theme(legend.position = "none" ) 
		
ggsave(paste0("p_Venn_ggVennDiagram_comps.pdf"), p_venn2, width=5, height=2.5)


##-----------------------------------------------------------------

## 3. use eulerr::venn
## Similar as venn.diagram 

p_venn3 <- eulerr::venn(LS_Comps)

pdf(paste0("p_Venn_eulerrvenn_comps.pdf"), width=4, height=3)
plot(p_venn3, quantities = list(type = c("counts"), cex=1.4),
		fill = mycol1, alpha=0.5, 
		main = paste0(" Venn "), 
		labels = list( labels=LS_Names, fontsize=15, font=1), adjust_labels=F ) 
dev.off()


##-----------------------------------------------------------------

## 4. use eulerr::euler
## Venn plot by scale 

p_venn4 <- eulerr::euler( LS_Comps )

## Note: if more than 4 groups, be better to use "ellipse" instead of circle 
p_venn4 <- eulerr::euler( LS_Comps, shape = "ellipse" )

pdf(paste0("p_Venn_euler_comps.pdf"), width=4.5, height=2.5)
print(plot(p_venn4, quantities = list(type = c("counts"), cex=1.4),
		fill = mycol1, alpha=0.5, 
		main = paste0("Venn"), 
		labels = list( labels=LS_Names, fontsize=16, font=1), adjust_labels=F ) )
dev.off()



##-----------------------------------------------------------------

## 5. Use ggvenn::ggvenn
   
library(ggvenn)

pdf("p_Venn_ggvenn.pdf", height=6, width=6)
ggvenn(LS_Comps, fill_color = c("red","cyan","green","orange"), text_size=4.5, set_name_size=6 ) 
dev.off()


   
   
   
