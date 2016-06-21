library(dplyr)
library(RColorBrewer)
library(readr)
library(ggplot2)

mc_no_insertion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_allC_no_insertion.tsv.gz", col_names = T)
mc_true_deletion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_allC_true_deletions.tsv.gz", col_names = T)

# sum rows for sort order
mc_no_insertion_2kb.sums <- mutate(mc_no_insertion_2kb, sums=rowSums(mc_no_insertion_2kb[,2:41]))
mc_true_deletion_2kb.sums <- mutate(mc_true_deletion_2kb, sums=rowSums(mc_true_deletion_2kb[,2:41]))

# Sort by total mC level
mc_no_insertion_2kb.sort <- arrange(mc_no_insertion_2kb.sums, desc(sums))
mc_true_deletion_2kb.sort <- arrange(mc_true_deletion_2kb.sums, desc(sums))

# categorical heatmap
# cat_hmap <- function(d) {
#   heatmap.2(d,Colv=NA,
#             #col=brewer.pal(9, "Dark2"),
#             col=brewer.pal(4, "Spectral"),
#             labRow = "",
#             labCol = "",
#             dendrogram = "none",
#             trace = "none",
#             par(cex.main=0.7,
#                 cex.lab=0.7,
#                 cex.axis=0.7
#             ),
#             cex.lab=0.7)
# }

scale_max <- function(m, l){
  m.floor <- m > l
  m[m.floor] <- l
  return(t(data.matrix(m)))
}
both <- rbind(mc_no_insertion_2kb.sort, mc_true_deletion_2kb.sort)
ins <- scale_max(mc_no_insertion_2kb.sort[,2:40], 0.5)
del <- scale_max(mc_true_deletion_2kb.sort[,2:40], 0.5)
both_scaled <- scale_max(both[,2:40], 0.5)
color <- colorRampPalette(brewer.pal(9,"Reds"))(100)

pdf("../Plots/heatmap_mc.pdf", height = 6, width = 3)
image(ins, col = color, main = "TE insertions")
image(del, col = color, main = "TE deletions")
image(both_scaled, col=color)
dev.off()

png("../Plots/heatmap_insertion.png", height = 6, width = 3, units = "in", res = 1200)
image(ins, col = color, main = "TE insertions")
dev.off()

png("../Plots/heatmap_deletion.png", height = 6, width = 3, units = "in", res = 1200)
image(del, col = color, main = "TE deletions")
dev.off()

png("../Plots/heatmap_ins_del.png", height = 6, width = 3, units = "in", res = 1200)
image(both_scaled, col=color)
dev.off()
