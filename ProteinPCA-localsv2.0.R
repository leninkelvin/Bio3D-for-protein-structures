
# Load the library
library(bio3d)
library(plotly)
library(ggplot2)
library(factoextra)
library(httr)
library(dplyr)
library(openxlsx)
library(DT)

# Get all PDB files in a directory
pdb_files <- list.files(path = "C:/Users/acher/Documents/PCA/ERS2/", 
                        pattern = "\\.pdb$", 
                        full.names = TRUE)

pdb_list <- lapply(pdb_files, read.pdb)

# Use pdbaln for alignment and analysis
pdbs <- pdbaln(pdb_files)


# Calculate sequence identity
pdbs$id <- substr(basename(pdbs$id),1,6)
seqidentity(pdbs)

heatmap(seqidentity(pdbs))

datatable(seqidentity(pdbs))

## Calculate RMSD
rmsd(pdbs, fit=TRUE)

core <- core.find(pdbs)

# See Figure 3.
col=rep("black", length(core$volume))
col[core$volume<2]="pink"; col[core$volume<1]="red"
plot(core, col=col)

core.inds <- print(core, vol=1.0)

write.pdb(xyz=pdbs$xyz[1,core.inds$xyz], file="quick_core.pdb")

xyz <- pdbfit( pdbs, core.inds )

rd <- rmsd(xyz)
heatmap(rd, scale = "column")


# RMSD clustering
hc.rd <- hclust(as.dist(rd))

pdbs$id <- substr(basename(pdbs$id), 1, 6)
hclustplot(hc.rd, k = 4, labels=pdbs$id, cex=1.2,
           ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)

# Perform PCA
#pc <- pca(pdbs)
# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)
pc <- pca.xyz(xyz[, gaps.pos$f.inds], use.svd = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pc$z)
colnames(pc_scores) <- paste0("PC", 1:ncol(pc_scores))

# Add labels (modify as needed)
pc_scores$label <- basename(pdb_files)

# Create interactive 2D plot
fig_2d <- plot_ly(pc_scores, 
                  x = ~PC1, 
                  y = ~PC2,
                  text = ~label,
                  type = 'scatter',
                  mode = 'markers',
                  marker = list(size = 10,
                                color = 'blue',
                                line = list(color = 'white', width = 1))) %>%
  layout(title = "PCA of Protein Structures",
         xaxis = list(title = paste0("PC1 (", round(pc$L[1]/sum(pc$L)*100, 1), "%)")),
         yaxis = list(title = paste0("PC2 (", round(pc$L[2]/sum(pc$L)*100, 1), "%)")),
         hovermode = 'closest')

fig_2d

# Create interactive 3D plot
fig_3d <- plot_ly(pc_scores,
                  x = ~PC1,
                  y = ~PC2,
                  z = ~PC3,
                  text = ~label,
                  type = 'scatter3d',
                  mode = 'markers',
                  marker = list(size = 5,
                                color = ~PC1,
                                colorscale = 'Viridis',
                                showscale = TRUE)) %>%
  layout(title = "3D PCA of Protein Structures",
         scene = list(
           xaxis = list(title = paste0("PC1 (", round(pc$L[1]/sum(pc$L)*100, 1), "%)")),
           yaxis = list(title = paste0("PC2 (", round(pc$L[2]/sum(pc$L)*100, 1), "%)")),
           zaxis = list(title = paste0("PC3 (", round(pc$L[3]/sum(pc$L)*100, 1), "%)"))
         ))

fig_3d