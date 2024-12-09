install.packages("ape")
install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("phyloseq")
library(ape)
library(ggtree)
library(phyloseq)

tree <- read.tree("Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")
pheno_df <- read.csv("Phenotypes_Zoonomia_3.csv")

pheno_df_subset <- pheno_df[!is.na(pheno_df$Homeotherm), c("Homeotherm", "Species")]
rownames(pheno_df_subset) <- pheno_df_subset$Species
pheno_df_subset$Homeotherm <- factor(pheno_df_subset$Homeotherm)

remove_species <- setdiff(tree$tip.label, pheno_df_subset$Species)
pruned_tree <- drop.tip(tree, remove_species)

plot(pruned_tree)

useTipColorsV <- c("red", "blue")
tipColorsV <- useTipColorsV[pheno_df_subset[pruned_tree$tip.label, "Homeotherm"]]

plot(pruned_tree, type = "fan", cex = 0.3, tip.color = tipColorsV)



