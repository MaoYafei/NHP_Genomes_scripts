rm(list = ls())
library(tidyverse)
library(PopGenome)


sparrows <- readData("./chr2/", format = "VCF", include.unknown = TRUE, FAST = TRUE)

sparrow_info <- read_delim("./pop.txt", delim = "\t")

populations <- split(sparrow_info$ind, sparrow_info$pop)
sparrows <- set.populations(sparrows, populations, diploid = T)
chr8 <- 242696138
window_size <- 20000
window_jump <- 10000

window_start <- seq(from = 1, to = chr8, by = window_jump)
window_stop <- window_start + window_size
sum(window_start > chr8)
sum(window_stop > chr8)

window_start <- window_start[which(window_stop < chr8)]
window_stop <- window_stop[which(window_stop < chr8)]
chr8 - window_stop[length(window_stop)]

sum(window_start > chr8)
sum(window_stop > chr8)

windows <- data.frame(start = window_start, stop = window_stop,mid = window_start + (window_stop-window_start)/2)

sparrows_sw <- sliding.window.transform(sparrows, width = 20000, jump = 10000, type = 2)
sparrows_sw <- diversity.stats(sparrows_sw, pi = TRUE)
#sparrows_sw <- F_ST.stats(sparrows_sw, mode = "nucleotide")
sparrows_sw <- sweeps.stats(sparrows_sw)

nd <- sparrows_sw@nuc.diversity.within/20000
pops <- c("AFR")
colnames(nd) <- paste0(pops, "_pi")
#fst <- t(sparrows_sw@nuc.F_ST.pairwise)
#dxy <- get.diversity(sparrows_sw, between = T)[[2]]/20000
#fst[fst < 0] <- 0 
#fst_dxy <- fst*100/dxy
CL <- sparrows_sw@CL

#x <- colnames(fst)
#x <- sub("pop1", pops[1], x)
#x <- sub("pop2", pops[2], x)
#x <- sub("/", "_", x)
#colnames(fst) <- paste0(x, "_fst")
#colnames(dxy) <- paste0(x, "_dxy")
#colnames(fst_dxy) <- paste0(x, "_fst/dxy")

sparrow_data <- as_tibble(data.frame(windows, nd,CL[,1]))
sparrow_data %>% select(contains("pi")) %>% summarise_all(mean)
#sparrow_data$Bonobo_Chimp_fst[sparrow_data$Bonobo_Chimp_fst < 0] <- 0 
sparrow_data$CL...1.[sparrow_data$CL...1. < -50] <- -50
#sparrow_data$CL...2.[sparrow_data$CL...2. < -50] <- -50

write.csv(sparrow_data,'chr2.csv')

pi_g <- sparrow_data %>% select(contains("pi")) %>% gather(key = "species", value = "pi")
a <- ggplot(pi_g, aes(species, pi)) + geom_boxplot() + theme_light() + xlab(NULL)
a
#a <- ggplot(sparrow_data, aes(mid/10^6, Bonobo_Chimp_fst)) + geom_line(colour = "red")
#a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
#a + theme_light()
hs <- sparrow_data %>% select(mid, AFR_pi,CL...1.,)
hs_g <- gather(hs, -mid, key = "stat", value = "value")

a <- ggplot(hs_g, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")


#library(tidyverse)
#library(PopGenome)
my_data <- hs_g %>% filter(mid > 86700000, mid < 87300000) #102000000, mid < 112000000 #mid > 85876743, mid < 88034905
a <- ggplot(my_data, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")
