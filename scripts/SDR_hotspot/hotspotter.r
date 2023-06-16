library(GenomicRanges)
library(rtracklayer)

gr_obj = import("all_nonsyn_final.bed" )
fai = read.delim('/net/eichler/vol26/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa.fai', header=FALSE, sep='\t')

seqlen <- c()

for (index in seq_len(nrow(fai))) {
	seqlen[fai[index, 'V1']] = fai[index, 'V2']
}


hotspotter <- function(gr, bw, pval=1e-3, num.trial=100) {
	set.seed(123) # fix seed for random permutations of bootstrapping
	## Iterate over chromosomes and calculate p-values
	pranges.list <- GenomicRanges::GRangesList()
		for (chrom in seqlevels(gr)) {
			grc <- gr[seqnames(gr)==chrom]
			if (length(grc)>1) {
				midpoints <- (start(grc)+end(grc))/2
				kde <- stats::density(midpoints,bw=bw,kernel='gaussian')
				# Random distribution of genomic events
				kde.densities <- numeric()
				for (i1 in seq_len(num.trial)) {
					midpoints.r <- round(stats::runif(length(midpoints),1,seqlen[chrom]))
					kde.r <- stats::density(midpoints.r,bw=bw,kernel='gaussian')
					kde.densities <- c(kde.densities, kde.r$y)
				}
				# Use ecdf to calculate p-values 
				p <- 1-stats::ecdf(kde.densities)(kde$y)
				pvalues <- data.frame(chromosome=chrom,start=kde$x,pvalue=p)
				# Make GRanges
				pvalues$end <- pvalues$start
				pvalues$chromosome <- factor(pvalues$chromosome, levels=seqlevels(gr))
				pvalues <- as(pvalues,'GRanges')
				seqlevels(pvalues) <- seqlevels(gr)
				suppressWarnings(
				seqlengths(pvalues) <- seqlengths(gr)[names(seqlengths(pvalues))]
				)
				# Resize from pointsize to bandwidth
				suppressWarnings(
				pvalues <- GenomicRanges::resize(pvalues, width=bw, fix='center')
				)
				pvalues <- trim(pvalues)
				## Find regions where p-value is below specification
				mask <- pvalues$pvalue <= pval
				rle.pvals <- rle(mask)
				rle.pvals$values <- cumsum(rle.pvals$values+1)
				pvalues$group <- inverse.rle(rle.pvals)
				if (length(which(mask))>0) {
					pvalues.split <- split(pvalues[mask],pvalues$group[mask])
					pranges <- unlist(endoapply(pvalues.split, function(x) { y <- x[1]; end(y) <- end(x)[length(x)]; y$pvalue <- min(x$pvalue); return(y) }))
					pranges$group <- NULL
					pranges$num.events <- GenomicRanges::countOverlaps(pranges, grc)
					## Make sure only non-zero counts are reported
					pranges <- pranges[pranges$num.events > 0]
					pranges.list[[chrom]] <- pranges
				}
			}
		}
	pranges <- unlist(pranges.list, use.names=FALSE)
	names(pranges) <- NULL
	return(pranges)
}

hotspots=hotspotter(gr_obj, 500000)

gr=hotspots
head(gr)

df <- data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1,
  ends=end(gr),
  names=c(rep(".", length(gr))),
  events=elementMetadata(gr)$num.events,
  pvalue=elementMetadata(gr)$pvalue,
  strands=strand(gr))

head(df)

write.table(df, "hotspot.bed",sep="\t",quote=FALSE,row.names=FALSE)
#rtracklayer::export.bed(hotspots, "hotspot.bed")
