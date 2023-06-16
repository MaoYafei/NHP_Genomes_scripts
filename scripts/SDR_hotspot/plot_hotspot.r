library(karyoploteR)
library("RColorBrewer")


custom.genome <- toGRanges("chm13_v1.0.len.bed")
kp <- plotKaryotype(genome = custom.genome,plot.type=2)

regions<-toGRanges('cenSat.merge.bed')
kpPlotRegions(kp, data=regions, col="grey5", border=NA, r0=-0.1, r1=-0.4,data.panel=1)

#regions<-toGRanges('chm13.draft_v1.0.SD.merge_noMY.bed')
#kpPlotRegions(kp, data=regions, col="red", border=NA, r0=-0.1, r1=-0.4,data.panel=1)



regions<-toGRanges("all_nonsyn_final.bed")
kpPlotDensity(kp, data=regions,col="black",window.size = 500000,data.panel=1,r0=0,r1=1.2,border=NA)

regions<-read.csv("hotspot.bed",sep='\t',header=T)
head(as.character(regions$col))

#kpRect(kp,chr=regions$seqnames,x0=regions$starts,x1=regions$ends,y0=0,y1=0.4,data.panel=2,col=as.character(regions$col),border=NA)
kpHeatmap(kp,chr=regions$seqnames,x0=regions$starts,x1=regions$ends,r0=1.2,r1=1.5,data.panel=1,y=regions$pvalue,col=rev(brewer.pal(9,"Reds")))
#brewer.pal(9,"Reds")
display.brewer.pal(n=9,name="Reds")
#kpPlotRegions(kp, data=regions,col=elementMetadata(regions)$col,r0=0,r1=0.4,data.panel=2)
