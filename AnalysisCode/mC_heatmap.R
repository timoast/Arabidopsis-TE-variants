library(dplyr)
library(gplots)
library(RColorBrewer)
library(readr)

mc_ins_2kb <- read_tsv("../ProcessedData/allC_ins_2kb.tsv", col_names = F)
mc_del_2kb <- read_tsv("../ProcessedData/allC_true_del_2kb.tsv", col_names = F)


### old code ###
### re-write for updated analysis ###

# sum rows for sort order
mc_all_ins_2kb_dist.sums <- mutate(mc_all_ins_2kb_dist, sums=rowSums(mc_all_ins_2kb_dist[,2:41]))


# categorical heatmap
cat_hmap <- function(d) {
  heatmap.2(d,Colv=NA,
            #col=brewer.pal(9, "Dark2"),
            col=brewer.pal(4, "Spectral"),
            labRow = "",
            labCol = "",
            dendrogram = "none",
            trace = "none",
            par(cex.main=0.7,
                cex.lab=0.7,
                cex.axis=0.7
            ),
            cex.lab=0.7)
}

# heatmap
hmap <- function(d) {
  heatmap.2(as.matrix(d),Colv=NA,col=brewer.pal(9, "Reds"),
            labRow = "",
            labCol = "",
            dendrogram = "none",
            margins = c(5,2),
            trace = "none",
            par(cex.main=0.7,
                cex.lab=0.7,
                cex.axis=0.7
            ),
            cex.lab=0.7)
}

# scale values to max 0.5
eu.floor <- eu > 0.5
eu[eu.floor] <- 0.5
het.floor <- het > 0.5
het[het.floor] <- 0.5

pdf("../Plots/heatmap_mc.pdf", height = 6, width = 3)
hmap(data)
dev.off()
