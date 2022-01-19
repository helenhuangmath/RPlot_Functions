## Create user-defined colors for ggplot 
## HH 01/2022


# R

library(RColorBrewer)

## Create multiple continous colors based on R brewer palette, use colorRampPalette() under RColorBrewer package 
colorRampPalette(brewer.pal(9, "Set1"))

## Usage example in bar plot: 
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

N=10
ggplot(meta, aes(x=celltype, fill=getPalette(N) )) +
	geom_bar(stat="count", width=0.8 )
  
ggplot(meta, aes(x=celltype, fill=sampleID ) ) +
  geom_bar(stat="count", width=0.8 ) +
	scale_fill_manual( values = getPalette(10) ) 





