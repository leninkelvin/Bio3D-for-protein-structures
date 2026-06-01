library(bio3d)
library(httr)
library(dplyr)
library(openxlsx)
library(DT)
library(heatmaply)
library(plotly)

# Download some example PDB files

{ids <- c('1l2j',
          '1nde',
          '1qkm',
          '1u3q',
          '1u3r',
          '1u3s',
          '1u9e',
          '1x76',
          '1x78',
          '1x7b',
          '1x7j',
          '1yy4',
          '1yye',
          '1zaf',
          '2fsz',
          '2giu',
          '2i0g',
          '2jj3',
          '2nv7',
          '2qtu',
          '2yjd',
          '2yly',
          '2z4b',
          '3oll',
          '3ols',
          '3omo',
          '3omp',
          '3omq',
          '4j24',
          '4j26',
          '4zi1',
          '5toa',
          '7xvy',
          '7xvz',
          '7xwp',
          '7xwq',
          '7xwr',
          '9ebl',
          '9ecf'
)} # All PDB

{ids <- c('1l2j_A',
          '1l2j_B',
          '1nde_A',
          '1qkm_A',
          '1u3q_A',
          '1u3q_B',
          '1u3q_C',
          '1u3q_D',
          '1u3r_A',
          '1u3r_B',
          '1u3s_A',
          '1u3s_B',
          '1u9e_A',
          '1u9e_B',
          '1x76_A',
          '1x76_B',
          '1x78_A',
          '1x78_B',
          '1x7b_A',
          '1x7b_B',
          '1x7j_A',
          '1x7j_B',
          '1yy4_A',
          '1yy4_B',
          '1yye_A',
          '1yye_B',
          '1zaf_A',
          '1zaf_B',
          '2fsz_A',
          '2fsz_B',
          '2giu_A',
          '2i0g_A',
          '2i0g_B',
          '2jj3_A',
          '2jj3_B',
          '2nv7_A',
          '2nv7_B',
          '2qtu_A',
          '2qtu_B',
          '2yjd_A',
          '2yjd_B',
          '2yly_A',
          '2yly_B',
          '2z4b_A',
          '2z4b_B',
          '3oll_A',
          '3oll_B',
          '3ols_A',
          '3ols_B',
          '3omo_A',
          '3omo_B',
          '3omp_A',
          '3omp_B',
          '3omq_A',
          '3omq_B',
          '4j24_A',
          '4j24_B',
          '4j24_C',
          '4j24_D',
          '4j26_A',
          '4j26_B',
          '4zi1_A',
          '5toa_A',
          '5toa_B',
          '7xvy_A',
          '7xvy_B',
          '7xvz_A',
          '7xvz_B',
          '7xwp_A',
          '7xwp_B',
          '7xwq_A',
          '7xwq_B',
          '7xwr_A',
          '7xwr_B',
          '9ebl_A',
          '9ebl_B',
          '9ebl_E',
          '9ebl_F',
          '9ecf_A',
          '9ecf_B',
          '9ecf_E',
          '9ecf_F',
          '9ecf_J',
          '9ecf_I',
          '9ecf_M',
          '9ecf_N'
)} # Curated LBD oligomers

{ids <- c('1u3r_C',
          '1u3r_D',
          '1u3s_C',
          '1u3s_D',
          '1u9e_C',
          '1u9e_D',
          '1x76_C',
          '1x76_D',
          '1x78_C',
          '1x78_D',
          '1x7b_C',
          '1x7b_D',
          '1x7j_C',
          '1x7j_D',
          '1yy4_C',
          '1yy4_D',
          '1yye_C',
          '1yye_D',
          '1zaf_C',
          '1zaf_D',
          '2nv7_C',
          '2nv7_D',
          '2yjd_C',
          '2yjd_D',
          '3oll_C',
          '3oll_D',
          '3ols_C',
          '3ols_D',
          '3omo_C',
          '3omo_D',
          '3omp_C',
          '3omp_D',
          '3omq_C',
          '3omq_D',
          '4j24_I',
          '4j24_J',
          '4j24_K',
          '4j24_L',
          '4j26_I',
          '4j26_J',
          '4zi1_B',
          '7xvy_C',
          '7xvy_D',
          '7xvz_C',
          '7xvz_D',
          '7xwp_C',
          '7xwp_D',
          '7xwq_C',
          '7xwq_D',
          '7xwr_C',
          '7xwr_D',
          '9ebl_C',
          '9ebl_D',
          '9ebl_G',
          '9ebl_H',
          '9ecf_K',
          '9ecf_P',
          '9ecf_O',
          '9ecf_L',
          '9ecf_D',
          '9ecf_C',
          '9ecf_H',
          '9ecf_G'
)} # Curated attached peptides
# NCOA1 black
# NCOA2 red
# NCOA5 blue
# Synthetic green
#  

raw.files <- get.pdb(ids)

anno <- pdb.annotate(ids)
datatable(anno)

write.xlsx(anno, 'annotaciones.xlsx')
annotation <- cbind(anno, color)
datatable(annotation)

# Extract and align the chains we are interested in
files <- pdbsplit(raw.files, ids)
pdbs <- pdbaln(files)

# Calculate sequence identity
pdbs$id <- substr(basename(pdbs$id),1,6)
seqidentity(pdbs)

heatmaply(seqidentity(pdbs))

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
write.table(rd, file = "RMSD_matrix.tab", sep = "\t", row.names = T)
hist(rd, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")

		# Save as PNG
			png("RMSDheatmap.png", width=1200, height=800)
			heatmap(rd, scale = "column")
		dev.off()

# RMSD clustering
hc.rd <- hclust(as.dist(rd))

pdbs$id <- substr(basename(pdbs$id), 1, 6)
hclustplot(hc.rd, k = 4, colors=annotation[, "color"], labels=pdbs$id, cex=1.2,
           ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)

	# Save as PNG
		png("RMSD_cluster.png", width=1200, height=800)
		hclustplot(hc.rd, k = 4, colors=annotation[, "color"], labels=pdbs$id, cex=1.2, 
		ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)
    dev.off()

# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

# Do PCA
pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds], use.svd = TRUE)
pc.xray

plot(pc.xray, col=annotation[, "color"])

    # Save as PNG
    png("PCAx4.png", width=1200, height=800)
    plot(pc.xray, col=annotation[, "color"])
    dev.off()

hc <- hclust(dist(pc.xray$z[,1:2]))
grps <- cutree(hc, h=30)
cols <- c("red", "green", "blue")[grps]
plot(pc.xray, pc.axes=1:2, col=annotation[, "color"])

      # Save as PNG
      png("PCAcluster.png", width=1200, height=800)
      plot(pc.xray, pc.axes=1:2, col=annotation[, "color"])
      dev.off()

# Dendrogram plot
names(cols) <- pdbs$id
hclustplot(hc, col=annotation[, "color"], ylab="Distance in PC Space", main="PC1-2", cex=0.9)

      # Save as PNG
      png("PCA_dendogram.png", width=1200, height=800)
      hclustplot(hc, col=annotation[, "color"], ylab="Distance in PC Space", main="PC1-2", cex=0.9)
      dev.off()

