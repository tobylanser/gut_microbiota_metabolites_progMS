# Load necessary packages
library(readr)
library(pheatmap)
library(ggplot2)

setwd("~/Dropbox (Partners HealthCare)/24-mmVEC/STOOL-trim")
# Set the path to your TSV file
file_path <- "conditionals-FORPAPER-TOBY.tsv"


data <- t(as.matrix(read.csv("conditionals-FORPAPER-TOBY.tsv", row.names = 1, sep = "\t")))
# Read the TSV file into a data frame
data <- read_tsv(file_path, col_types = cols(),)

# Optionally, you might want to set row names (if your data frame doesn't have any yet)
rownames(data) <- data$Compound_Source
data2 <- t(data)
colnames(data2) <- rownames(data)
data2 <- data2[-1, ]
data2 <- as.data.frame(data2)
colnames(data2) <- NULL

data2 <- as.numeric(data2)
# Create the heatmap
heat <- pheatmap(data, show_colnames = T, angle_col = "45")
?pheatmap
?heatmap
ggsave(paste0("toby_conditionals_FORPAPE_heatmap.pdf"), device = "pdf", 
       width = 16, height = 12, units = "in", dpi = 300, plot = heat)

data2 <- as.matrix(data2)
heatmap(data)

# Save the heatmap to a file (optional)
# pheatmap(data, filename = "heatmap.png")
