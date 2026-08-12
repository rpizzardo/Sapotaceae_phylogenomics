#### Getting bootstrap values per tree part - gene trees ####
#
# This script roots, ladderizes, and time-calibrates the
# individual nuclear (353 loci + ITS) and plastome gene trees for
# Sapotaceae, then summarizes how bootstrap values are distributed across each
# tree's part and averaging within each quartile, per gene.
#
# Expected input files:
#   raw_treefiles/indiv_gene_trees_sapo_ufbt_353plusITS.treefile - individual
#       nuclear (353 loci + ITS) gene trees
#   raw_treefiles/indv_gene_trees_sapo_cp.treefile                - individual
#       chloroplast gene trees
#   sapo_a353_ITS_gene_list.txt  - one gene/locus name per line, matching the
#       order of trees in the 353+ITS multiPhylo object
#   cp_genelist.txt              - one gene/locus name per line, matching the
#       order of trees in the chloroplast multiPhylo object
#
# Winter 2025
# RC Pizzardo
# =============================================================================

# Set working directory
wd <- ("your directory")
setwd(wd)

# Load necessary libraries
library(ape)         # For phylogenetic analysis
library(parallel)    # For potential parallel processing
library(phytools)
# -----------------------------------------------------------------------------
# Helper function: reroot a tree on a given node, then ladderize it
# -----------------------------------------------------------------------------
# Function to reroot and ladderize a tree based on a given root node
reroot_and_ladderize_tree <- function(inputTree, rootNode) {
  if (is.rooted(inputTree) == TRUE) {               # If tree is already rooted, unroot it
    inputTree <- unroot(inputTree)
  }
  rootedTree <- reroot(inputTree, rootNode, 
                       position=0.5 * inputTree$edge.length[which(inputTree$edge[,2] == rootNode)])  # Re-root at midpoint
  ladderizedTree <- ladderize(rootedTree, right=TRUE)  # Ladderize tree for consistent visualization
  return(ladderizedTree)
}

# =============================================================================
# 1. Load gene trees and their name lists
# =============================================================================

# Load raw nuclear and chloroplast gene trees
nuclear.353ITS <- read.tree("raw_treefiles/indiv_gene_trees_sapo_ufbt_353plusITS.treefile")
cp.79 <- read.tree("raw_treefiles/indv_gene_trees_sapo_cp.treefile")

# Load gene name lists for labeling
names.n <- read.table("sapo_a353_ITS_gene_list.txt")
names.cp <- read.table("cp_genelist.txt")

# Assign gene names to nuclear and cp trees
names(nuclear.353ITS) <- names.n$V1
names(cp.79) <- names.cp$V1

# =============================================================================
# 2. Root each gene tree on its outgroup, and ladderize
# =============================================================================
# Outgroup taxa belong to the families Marcgraviaceae or Tetrameristaceae.
# If a gene has exactly one outgroup tip present, root directly on it; if it
# has several, root on their most recent common ancestor (MRCA); if it has
# none, skip that gene and report it so it can be rooted manually later.

# --- Nuclear gene trees ------------------------------------------------------
nuclear.353ITS.root <- list()

for (g in 1:length(nuclear.353ITS)) {
  gt <- nuclear.353ITS[[g]]
  p <- grep("Marcgraviaceae_|Tetrameristaceae_", gt$tip.label)         # Identify outgroup tips
  outgroup <- gt$tip.label[p]
  foundOGs <- gt$tip.label[!is.na(match(gt$tip.label,outgroup))]
  if (length(foundOGs) == 1) {
    rrNode <- which(gt$tip.label == foundOGs[1])
  } else if (length(foundOGs) == 0){
    print(paste("Gene", names.n[g,1], "is missing all outgroups. Try rerooting manually.", sep =" "))
    next
  } else {
    rrNode <- getMRCA(gt, tip=foundOGs)   # If multiple outgroups, root at their MRCA
  }
  gt.rerooted <- reroot_and_ladderize_tree(gt, rrNode)
  nuclear.353ITS.root[[g]] <- gt.rerooted
}

class(nuclear.353ITS.root) <- "multiPhylo"
names(nuclear.353ITS.root) <- names(nuclear.353ITS)

# --- Chloroplast gene trees --------------------------------------------------
# Root all chloroplast gene trees similarly
cp.79.root <- list()

for (g in 1:length(cp.79)) {
  gt <- cp.79[[g]]
  p <- grep("Marcgraviaceae_|Tetrameristaceae_", gt$tip.label)
  outgroup <- gt$tip.label[p]
  foundOGs <- gt$tip.label[!is.na(match(gt$tip.label,outgroup))]
  if (length(foundOGs) == 1) {
    rrNode <- which(gt$tip.label == foundOGs[1])
  } else if (length(foundOGs) == 0){
    # NOTE: reports the nuclear gene-name list (names.n) here rather than
    # names.cp - likely a copy-paste carryover from the nuclear loop above,
    # so the gene name printed for a missing cp outgroup may be wrong.
    print(paste("Gene", names.n[g,1], "is missing all outgroups. Try rerooting manually.", sep =" "))
    next
  } else {
    rrNode <- getMRCA(gt, tip=foundOGs)
  }
  gt.rerooted <- reroot_and_ladderize_tree(gt, rrNode)
  cp.79.root[[g]] <- gt.rerooted
}

class(cp.79.root) <- "multiPhylo"
names(cp.79.root) <- names(cp.79)

# =============================================================================
# 3. Time-calibrate each rooted gene tree with chronos
# =============================================================================
# `chronos()` fits a penalized-likelihood relaxed-clock model to convert each
# rooted tree into an ultrametric (time-calibrated) tree, from which node
# ages ("branching times") can be extracted.

# --- Nuclear trees ------------------------------------------------------------
nuclear.353ITS.dated <- list()
bt <- list()
for (i in 1:length(nuclear.353ITS.root)) {
  tr <- nuclear.353ITS.root[[i]]
  if(length(tr) == 0){
    print(paste("Tree", names.n[i,1], "is missing.", sep =" "))
  }
  c <- chronos(tr)                                 # Estimate ultrametric tree
  nuclear.353ITS.dated[[i]] <- c
  bt[[i]] <- branching.times(c)                    # Extract branching times
  print(paste("Tree", names.n[i,1], "done!"))
}

# Store and annotate dated nuclear trees and their branching times
dated_trees <- nuclear.353ITS.dated
bt.nuclear <- bt
class(dated_trees) <- "multiPhylo"
names(bt.nuclear) <- names(nuclear.353ITS.root)

# --- Plot: nuclear branching-time distribution in 4 relative-depth quartiles -
# For each nuclear gene tree, split its branching times into four bins based
# on their relative position between the root and the tips, and
# boxplot the values within each bin. One boxplot panel is produced per gene.
pdf(file = "test_nuclear.pdf")
par(mfrow=c(2,2))
for (i in 1:length(bt.nuclear)) {
  b.now <- bt.nuclear[i]
  b.now <- b.now[[1]][]
  if (is.null(b.now) == TRUE){
    next
  }
  value1 <- as.numeric(names(b.now[b.now <= 0.25]))
  value.no.na1 <- value1[!is.na(value1)]
  value2 <- as.numeric(names(b.now[b.now >= 0.25 & b.now <= 0.5]))
  value.no.na2 <- value2[!is.na(value2)]
  value3 <- as.numeric(names(b.now[b.now >= 0.5 & b.now <= 0.75]))
  value.no.na3 <- value3[!is.na(value3)]
  value4 <- as.numeric(names(b.now[b.now >= 0.75 & b.now <= 1]))
  value.no.na4 <- value4[!is.na(value4)]
  
  boxplot(value.no.na1, value.no.na2, value.no.na3, value.no.na4, main = names(bt.nuclear[i]), 
          names = c("1", "2", "3", "4"))
}
dev.off()

# --- Chloroplast trees --------------------------------------------------------
# Date all rooted chloroplast gene trees
cp.79.dated <- list()
bt <- list()
for (i in 1:length(cp.79.root)) {
  tr <- cp.79.root[[i]]
  if(length(tr) == 0){
    print(paste("Tree", names.cp[i,1], "is missing.", sep =" "))
    next
  }
  c <- chronos(tr)
  cp.79.dated[[i]] <- c
  bt[[i]] <- branching.times(c)
  print(paste("Tree", names.cp[i,1], "done!"))
}

# Store and annotate dated cp trees and branching times
dated_trees_cp <- cp.79.dated
bt.cp <- bt
class(dated_trees_cp) <- "multiPhylo"
names(bt.cp) <- names(cp.79.root)

# --- Plot: chloroplast branching-time distribution in 4 relative-depth quartiles
# Same idea as the nuclear plot above.
pdf(file = "test_cp.pdf")
par(mfrow=c(2,2))
for (i in 1:length(bt.cp)) {
  b.now <- bt.cp[i]
  b.now <- b.now[[1]][]
  if (is.null(b.now) == TRUE){
    next
  }
  value1 <- as.numeric(sapply(strsplit(names(b.now[b.now <= 0.25]), "/"), `[`, 1))
  value.no.na1 <- value1[!is.na(value1)]
  value2 <- as.numeric(sapply(strsplit(names(b.now[b.now >= 0.25 & b.now <= 0.5]), "/"), `[`, 1))
  value.no.na2 <- value2[!is.na(value2)]
  value3 <- as.numeric(sapply(strsplit(names(b.now[b.now >= 0.5 & b.now <= 0.75]), "/"), `[`, 1))
  value.no.na3 <- value3[!is.na(value3)]
  value4 <- as.numeric(sapply(strsplit(names(b.now[b.now >= 0.75 & b.now <= 1]), "/"), `[`, 1))
  value.no.na4 <- value4[!is.na(value4)]
  boxplot(value.no.na1, value.no.na2, value.no.na3, value.no.na4, main = names(bt.cp[i]), 
          names = c("1", "2", "3", "4"))
}
dev.off()

# =============================================================================
# 4. Summarize mean bootstrap value per quartile, per gene
# =============================================================================

# --- Chloroplast: mean bootstrap value per quartile ---------------------------------
# Calculate mean node ages per tree quartile for cp trees
m <- matrix(nrow = 4, ncol = length(dated_trees_cp))
for (i in 1:length(bt.cp)) {
  b.now <- bt.cp[i]
  b.now <- b.now[[1]][]
  if (is.null(b.now) == TRUE){
    next
  }
  value1 <- as.numeric(sapply(strsplit(names(b.now[b.now <= 0.25]), "/"), `[`, 1))
  value.no.na1 <- value1[!is.na(value1)]
  m[1,i] <- mean(value.no.na1)
  value2 <- as.numeric(sapply(strsplit(names(b.now[b.now >= 0.25 & b.now <= 0.5]), "/"), `[`, 1))
  value.no.na2 <- value2[!is.na(value2)]
  m[2,i] <- mean(value.no.na2)
  value3 <- as.numeric(sapply(strsplit(names(b.now[b.now >= 0.5 & b.now <= 0.75]), "/"), `[`, 1))
  value.no.na3 <- value3[!is.na(value3)]
  m[3,i] <- mean(value.no.na3)
  value4 <- as.numeric(sapply(strsplit(names(b.now[b.now >= 0.75 & b.now <= 1]), "/"), `[`, 1))
  value.no.na4 <- value4[!is.na(value4)]
  m[4,i] <- mean(value.no.na4)
}
m  # rows = quartiles (1 = root-ward ... 4 = tip-ward), columns = cp genes

# --- Nuclear: mean bootstrap value per quartile --------------------------------------
# Calculate mean bootstrap value per tree quartile for nuclear trees
m2 <- matrix(nrow = 4, ncol = length(dated_trees))
for (i in 1:length(bt.nuclear)) {
  b.now <- bt.nuclear[i]
  b.now <- b.now[[1]][]
  value1 <- as.numeric(names(b.now[b.now <= 0.25]))
  value.no.na1 <- value1[!is.na(value1)]
  m2[1,i] <- mean(value.no.na1)
  value2 <- as.numeric(names(b.now[b.now >= 0.25 & b.now <= 0.5]))
  value.no.na2 <- value2[!is.na(value2)]
  m2[2,i] <- mean(value.no.na2)
  value3 <- as.numeric(names(b.now[b.now >= 0.5 & b.now <= 0.75]))
  value.no.na3 <- value3[!is.na(value3)]
  m2[3,i] <- mean(value.no.na3)
  value4 <- as.numeric(names(b.now[b.now >= 0.75 & b.now <= 1]))
  value.no.na4 <- value4[!is.na(value4)]
  m2[4,i] <- mean(value.no.na4)
}
m2  # rows = quartiles (1 = root-ward ... 4 = tip-ward), columns = nuclear genes

# =============================================================================
# 5. Combined summary plot and file outputs
# =============================================================================
pdf(file = "genes_per_tree_part.pdf")
matplot(m, type = "b",pch=1,col = 1:5, main = "CP", xlab = "Tree part", ylab = "Mean bootstrap value")
matplot(m2, type = "b",pch=1,col = 1:5, main = "Nuclear", xlab = "Tree part", ylab = "Mean bootstrap value")
dev.off()

# --- Save chloroplast summary table and dated trees --------------------------
write.csv(m, "cp_genes_per_tr_part.csv")
# Keep only genuinely dated trees (drop any failed/placeholder entries)
valid_trees <- Filter(function(tr) inherits(tr, "phylo") && !is.null(tr$Nnode), dated_trees_cp)
write.tree(valid_trees, "dated_gene_tr_cp.tre")

# --- Save nuclear summary table and dated trees -------------------------------
write.csv(m2, "nuclear_genes_per_tr_part.csv")
valid_trees <- Filter(function(tr) inherits(tr, "phylo") && !is.null(tr$Nnode), dated_trees)
write.tree(valid_trees, "dated_gene_tr_nuclear.tre")

### THE END ###
