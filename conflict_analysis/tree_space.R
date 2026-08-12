# =============================================================================
# tree_space_and_ boxplot.R
# -----------------------------------------------------------------------------
# Compares gene-tree support and topological congruence between the nuclear
# (Angiosperms353 + ITS) and plastome gene-tree datasets for Sapotaceae.
#
# The script has two main parts:
#
#   1. BOOTSTRAP SUPPORT BOXPLOTS
#      For each dataset (nuclear, plastome), find the 3 genes with the
#      highest and 3 genes with the lowest mean bootstrap support, and plot
#      their node-support distributions alongside the full concatenated tree.
#
#   2. TREE SPACE (topological distance ordination)
#      Using the `treespace` package, compute Robinson-Foulds distances
#      between all individual gene trees (after pruning them to a common set
#      of taxa) and visualize the resulting ordination, colored by data
#      partition (nuclear vs. plastome), with convex hulls per group.
#      (353 nuclear loci + ITS, vs. plastid)
#
# Expected input files (relative to the working directory set below):
#   raw_treefiles/indiv_gene_trees_sapo_ufbt_353plusITS.treefile  - individual
#       nuclear (353 loci + ITS) gene trees with UFBoot support values
#   raw_treefiles/indv_gene_trees_sapo_cp.treefile                - individual
#       plastid gene trees with bootstrap support values
#   raw_treefiles/sapo_concat_ML_353plusITS.treefile              - concatenated
#       nuclear ML tree (353+ITS), used for the "full" support distribution
#   raw_treefiles/sapo_concat_cp_ML.treefile                      - concatenated
#       plastid ML tree, used for the "full" support distribution
#   sapo_a353_ITS_gene_list.txt  - one gene/locus name per line, in the same
#       order as the trees in the 353+ITS multiPhylo object (used to name them)
#   cp_genelist.txt              - one gene/locus name per line, in the same
#       order as the trees in the plastid multiPhylo object
#
#
# Winter 2025
# RC Pizzardo
# =============================================================================

library(treespace)   # tree-space ordination based on topological distances
library(ape)          # phylogenetic tree handling
library(ggplot2)       # plotting (used for the convex-hull tree-space plot)
library(dplyr)

# -----------------------------------------------------------------------------
# 0. Working directory and input files
# -----------------------------------------------------------------------------
wd <- "your working directory"
setwd(wd)

# Individual gene trees (multiPhylo objects, one tree per locus)
nuclear.353ITS <- read.tree("raw_treefiles/indiv_gene_trees_sapo_ufbt_353plusITS.treefile")
cp.79          <- read.tree("raw_treefiles/indv_gene_trees_sapo_cp.treefile")

# Concatenated ML trees
tr <- read.tree("raw_treefiles/sapo_concat_ML_353plusITS.treefile")
cp <- read.tree("raw_treefiles/sapo_concat_cp_ML.treefile")

# Gene name lists, used to label each tree in the multiPhylo objects above
names.n   <- read.table("sapo_a353_ITS_gene_list.txt")
names.cp  <- read.table("cp_genelist.txt")

names(nuclear.353ITS) <- names.n$V1
names(cp.79)          <- names.cp$V1

# =============================================================================
# 1. Bootstrap support boxplots
# =============================================================================

# -----------------------------------------------------------------------------
# 1a. Nuclear (353 + ITS): identify highest- and lowest-support genes
# -----------------------------------------------------------------------------
# Mean node support (across all internal nodes) for each individual gene tree
m <- numeric(length = 0)
for (i in 1:length(nuclear.353ITS)) {
  x <- nuclear.353ITS[[i]]$node.label
  m[i] <- mean(as.numeric(x)[!is.na(as.numeric(x))])
}

# 3 genes with the highest mean support, and 3 with the lowest
top    <- names(nuclear.353ITS)[match(sort(m, decreasing = TRUE)[1:3], m)]
bottom <- names(nuclear.353ITS)[which(m %in% sort(m)[1:3])]

# Clean up the concatenated tree's node labels: IQ-TREE-style labels can look
# like "95/87" (SH-aLRT/UFBoot) - keep only the value before the "/"
tr$node.label <- gsub("/", " ", tr$node.label)
tr$node.label <- sub(" .*", "", tr$node.label)

# Support-value vectors: full concatenated tree, plus each of the 6 selected genes
full  <- as.numeric(tr$node.label)
data1 <- as.numeric(nuclear.353ITS[[top[1]]]$node.label)[!is.na(as.numeric(nuclear.353ITS[[top[1]]]$node.label))]
data2 <- as.numeric(nuclear.353ITS[[top[2]]]$node.label)[!is.na(as.numeric(nuclear.353ITS[[top[2]]]$node.label))]
data3 <- as.numeric(nuclear.353ITS[[top[3]]]$node.label)[!is.na(as.numeric(nuclear.353ITS[[top[3]]]$node.label))]
data4 <- as.numeric(nuclear.353ITS[[bottom[3]]]$node.label)[!is.na(as.numeric(nuclear.353ITS[[bottom[3]]]$node.label))]
data5 <- as.numeric(nuclear.353ITS[[bottom[2]]]$node.label)[!is.na(as.numeric(nuclear.353ITS[[bottom[2]]]$node.label))]
data6 <- as.numeric(nuclear.353ITS[[bottom[1]]]$node.label)[!is.na(as.numeric(nuclear.353ITS[[bottom[1]]]$node.label))]

# Quick look: full tree + the 3 highest- and 3 lowest-support nuclear genes
boxplot(full, data1, data2, data3, data4, data5, data6,
        names = c("full", top[1], top[2], top[3], bottom[3], bottom[2], bottom[1]),
        col = rainbow(7))

# -----------------------------------------------------------------------------
# 1b. Plastid (cp): identify highest- and lowest-support genes
# -----------------------------------------------------------------------------
# Clean up node labels the same way as for the nuclear trees above
for (i in 1:length(cp.79)) {
  t <- cp.79[[i]]
  t$node.label <- gsub("/", " ", t$node.label)
  t$node.label <- sub(" .*", "", t$node.label)
  cp.79[[i]]$node.label <- t$node.label
}

m2 <- numeric(length = 0)
for (i in 1:length(cp.79)) {
  x <- cp.79[[i]]$node.label
  m2[i] <- mean(as.numeric(x)[!is.na(as.numeric(x))])
}

top2    <- names(cp.79)[match(sort(m2, decreasing = TRUE)[1:3], m2)]
bottom2 <- names(cp.79)[which(m2 %in% sort(m2)[1:3])]

cp$node.label <- gsub("/", " ", cp$node.label)
cp$node.label <- sub(" .*", "", cp$node.label)

full_cp <- as.numeric(cp$node.label)
data1cp <- as.numeric(cp.79[[top2[1]]]$node.label)[!is.na(as.numeric(cp.79[[top2[1]]]$node.label))]
data2cp <- as.numeric(cp.79[[top2[2]]]$node.label)[!is.na(as.numeric(cp.79[[top2[2]]]$node.label))]
data3cp <- as.numeric(cp.79[[top2[3]]]$node.label)[!is.na(as.numeric(cp.79[[top2[3]]]$node.label))]
data4cp <- as.numeric(cp.79[[bottom2[3]]]$node.label)[!is.na(as.numeric(cp.79[[bottom2[3]]]$node.label))]
data5cp <- as.numeric(cp.79[[bottom2[2]]]$node.label)[!is.na(as.numeric(cp.79[[bottom2[2]]]$node.label))]
data6cp <- as.numeric(cp.79[[bottom2[1]]]$node.label)[!is.na(as.numeric(cp.79[[bottom2[1]]]$node.label))]

# Quick look: full plastid tree + the 3 highest- and 3 lowest-support plastid genes
boxplot(full_cp, data1cp, data2cp, data3cp, data4cp, data5cp, data6cp,
        names = c("cp", top2[1], top2[2], top2[3], bottom2[3], bottom2[2], bottom2[1]),
        col = rainbow(7))

# -----------------------------------------------------------------------------
# 1c. Combined publication figure: full trees + top-3 genes, nuclear vs. plastid
# -----------------------------------------------------------------------------
pdf(file = "plots/boxplots_markers.pdf", width = 6, height = 10)
boxplot(full, full_cp, data1, data2, data3, data1cp, data2cp, data3cp,
        names = c("Nuclear", "Plastome", top[1], top[2], top[3], top2[1], top2[2], top2[3]),
        col = "gray95", xlab = "Bootstraps values", ylab = "Markers", pch = 20,
        las = 1, notch = FALSE, boxwex = 0.6, horizontal = TRUE)
dev.off()


# =============================================================================
# 2. Tree space: nuclear (353 + ITS) vs. plastid
# =============================================================================
# Robinson-Foulds tree-space analysis requires every input tree to share the
# exact same set of tips, so genes/trees with too few taxa are dropped, and
# the remaining trees are pruned down to their common taxon set.

# -----------------------------------------------------------------------------
# 2a. Drop gene trees with too few tips
# -----------------------------------------------------------------------------
# Nuclear: flag any gene tree with fewer than 250 tips
genes.to.drop.n <- numeric()
for (i in 1:length(nuclear.353ITS)) {
  x <- length(nuclear.353ITS[[i]]$tip.label)
  if (x < 250) {
    print(paste("tree", i, "has less than 250 tips.", "NTips = ", x))
    genes.to.drop.n <- c(genes.to.drop.n, i)
  }
}
# Tree 245 has only 54 samples (gene 6514) and is excluded.
# See `genes.to.drop.n` for the full list of excluded trees.

# Plastid: flag any gene tree with fewer than 200 tips
genes.to.drop.c <- numeric()
for (i in 1:length(cp.79)) {
  x <- length(cp.79[[i]]$tip.label)
  if (x < 200) {
    print(paste("tree", i, "has less than 250 tips.", "NTips = ", x))
    genes.to.drop.c <- c(genes.to.drop.c, i)
  }
}
# Trees 11, 27, and 33 have only 79, 80, and 93 samples (genes infA, petG, psaI)
# and are excluded. See `genes.to.drop.c` for the full list of excluded trees.

nuclear.353.d <- nuclear.353ITS[-genes.to.drop.n]  # 289 gene trees remain
cp.79.d       <- cp.79[-genes.to.drop.c]           # 12 gene trees remain

# -----------------------------------------------------------------------------
# 2b. Find the taxa shared by every remaining tree, and prune to that set
# -----------------------------------------------------------------------------
tips.list.n <- list()
for (i in 1:length(nuclear.353.d)) {
  tp <- nuclear.353.d[[i]]$tip.label
  tips.list.n[[i]] <- tp
}
length(tips.list.n)  # 290

tips.list.c <- list()
for (i in 1:length(cp.79.d)) {
  tp <- cp.79.d[[i]]$tip.label
  tips.list.c[[i]] <- tp
}
length(tips.list.c)  # 69

tips.list <- append(tips.list.n, tips.list.c)
length(tips.list)  # 359

to.keep <- Reduce(intersect, tips.list)  # 65 tips shared by every tree

for (i in 1:length(nuclear.353.d)) {  # prune nuclear trees down to the shared taxa
  to.drop <- nuclear.353.d[[i]]$tip.label[nuclear.353.d[[i]]$tip.label %in% to.keep == FALSE]
  nuclear.353.d[[i]] <- drop.tip(nuclear.353.d[[i]], to.drop)
}

for (i in 1:length(cp.79.d)) {  # prune plastid trees down to the shared taxa
  to.drop <- cp.79.d[[i]]$tip.label[cp.79.d[[i]]$tip.label %in% to.keep == FALSE]
  cp.79.d[[i]] <- drop.tip(cp.79.d[[i]], to.drop)
}

all <- append(nuclear.353.d, cp.79.d)

# -----------------------------------------------------------------------------
# 2c. Compute and plot the tree-space ordination
# -----------------------------------------------------------------------------
tspace.all <- treespace(all, method = "RF", nf = 10)  # Robinson-Foulds tree-space, 10 axes
head(tspace.all$pco$li)  # ordination coordinates (principal coordinates) per tree
plot(tspace.all$pco$li[, 1], tspace.all$pco$li[, 2])  # quick look at axes 1 vs. 2
round(tspace.all$pco$eig / sum(tspace.n$pco$eig), 3)

# Color palette: one color per data partition (nuclear vs. plastid)
partitions <- c("nuclear", "plastome")
palettes_hcl <- hcl.pals()
colors <- hcl.colors(12, palette = "SunsetDark", alpha = 1)[c(1, 6)]
names(colors) <- partitions

# Group-membership factor: first 290 trees are nuclear, next 69 are plastid
y  <- gl(1, 290, labels = "nuclear")
yy <- gl(1, 69,  labels = "plastome")

# Scatter plot of the first two ordination axes, colored by partition
s.class(tspace.all$pco$li, fac = as.factor(c(y, yy)), cellipse = 0, cstar = 1,
        label = NULL, pch = 19, cgrid = 1, cpoint = 0.5,
        col = colors[levels(as.factor(c(y, yy)))])
# Convex hull around each partition's points
s.chull(tspace.all$pco$li, fac = as.factor(c(y, yy)),
        add.plot = TRUE, clabel = 0.7, col = colors[levels(as.factor(c(y, yy)))],
        optchull = 1)

# Save the same plot to PDF for publication
pdf(file = "plots/tree_space_353plusITS.pdf", width = 10, height = 10)
s.class(tspace.all$pco$li, fac = as.factor(c(y, yy)), cellipse = 0, cstar = 1,
        label = NULL, pch = 19, cgrid = 1, cpoint = 0.5,
        col = colors[levels(as.factor(c(y, yy)))])
s.chull(tspace.all$pco$li, fac = as.factor(c(y, yy)),
        add.plot = TRUE, clabel = 0.7, col = colors[levels(as.factor(c(y, yy)))],
        optchull = 1)
dev.off()

# -----------------------------------------------------------------------------
# 2d. ggplot2 version, with convex hulls drawn via dplyr
# -----------------------------------------------------------------------------
x <- cbind(tspace.all$pco$li, c(y, yy))
colnames(x)[11] <- "markers"
hull <- x %>% group_by(markers) %>% slice(chull(A1, A2))  # convex-hull vertices per group

p <- ggplot(x, aes(x = A1, y = A2)) +
  geom_point(size = 2, shape = 1, aes(color = markers)) +  # grey edges
  theme_bw()
p + geom_polygon(data = hull, alpha = 0, aes(fill = markers, color = markers)) +
  scale_color_manual(values = c("#7D1D67FF", "#F34A70FF"))


### THE END ###