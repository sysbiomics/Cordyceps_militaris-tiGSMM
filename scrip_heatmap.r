## Load in data
data <- read.csv("LS.csv", header=T, row.names=1)
data_matrix <- as.matrix(data)

library(pheatmap)
library(RColorBrewer)


display.brewer.all()
cols <- brewer.pal(11, "RdBu")

pheatmap(data_matrix,
         color = colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100),
         show_rownames = T,
         cluster_cols = F,
         cluster_rows = F,
)
res <- pheatmap(data_matrix, scale = "row",
                color = colorRampPalette(rev(brewer.pal(n=9,name="RdBu")))(50),
                show_rownames = T,
                cluster_cols = F,
                cluster_rows = F,
)

a <-data.frame(cutree(res$tree_row, k = 2))
write.csv(as.data.frame(a),file = "Clusters_DEG2_25_15.csv")

