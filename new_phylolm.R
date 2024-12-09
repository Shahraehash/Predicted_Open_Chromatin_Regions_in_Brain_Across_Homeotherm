options(repos = c(CRAN = "https://cloud.r-project.org/"))

# List of required packages
packages <- c("ape", "BiocManager", "phylolm")

# Install packages that are not already installed
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load the installed packages
library(ape)
library(BiocManager)  # Load BiocManager to use for Bioconductor packages

# Install Bioconductor packages
bioc_packages <- c("ggtree", "phyloseq")
BiocManager::install(bioc_packages)

# Load the Bioconductor packages
library(ggtree)
library(phyloseq)
library(phylolm)


#Note: did not use wget, but loaded it into ssh'd Jupyter notebook
astrocytePredictionsPath <- "Astrocyte/240_predictions_MatrixStacked.tsv"
astrocytePredictions <- read.delim(astrocytePredictionsPath, header = FALSE, sep = "\t")
astrocytePredictionsNamesPath <- "Astrocyte/240_predictions_NamesList.txt"
astrocytePredictionsNames <- read.delim(astrocytePredictionsNamesPath, header = FALSE, sep = "\t")

vipPredictionsPath <- "VIP/240_predictions_MatrixStacked.tsv"
vipPredictions <- read.delim(vipPredictionsPath, header = FALSE, sep = "\t")
vipPredictionsNamesPath <- "VIP/240_predictions_NamesList.txt"
vipPredictionsNames <- read.delim(vipPredictionsNamesPath, header = FALSE, sep = "\t")

# Extract the species names from astrocytePredictionsNames and make them column names
species_names <- astrocytePredictionsNames[[1]]
# Check if the length matches the number of columns (excluding the first)
if (length(species_names) == (ncol(astrocytePredictions) - 1)) {
  # Rename the columns in astrocytePredictions (excluding the first column)
  colnames(astrocytePredictions)[2:ncol(astrocytePredictions)] <- species_names
  colnames(vipPredictions)[2:ncol(vipPredictions)] <- species_names
  # Print the new column names to verify
  #print(colnames(astrocytePredictions))
} else {
  stop("The number of species names does not match the number of columns in astrocytePredictions.")
}


# Transpose the astrocytePredictions data frame so species are rows and peaks are columns
astrocytePredictions <- t(astrocytePredictions)
# Convert it back to a data frame
astrocytePredictions <- as.data.frame(astrocytePredictions)
#Fix it so peaks are actually the column names
colnames(astrocytePredictions) <- astrocytePredictions[1, ]
astrocytePredictions <- astrocytePredictions[-1, ]

# Transpose the astrocytePredictions data frame so species are rows and peaks are columns
vipPredictions <- t(vipPredictions)
# Convert it back to a data frame
vipPredictions <- as.data.frame(vipPredictions)
#Fix it so peaks are actually the column names
colnames(vipPredictions) <- vipPredictions[1, ]
vipPredictions <- vipPredictions[-1, ]

common_columns <- intersect(names(astrocytePredictions), names(vipPredictions))
length(common_columns)


#Keep only common columns
#astrocytePredictions <- astrocytePredictions[, common_columns, drop = FALSE]
#vipPredictions <- vipPredictions[, common_columns, drop = FALSE]

#To get random 10,000 that are not in common columns
all_astrocyte_indices <- 1:ncol(astrocytePredictions)
#non_common_astrocyte_indices <- setdiff(all_astrocyte_indices, common_columns)
#subset_astrocyte_indices <- sample(non_common_astrocyte_indices, size = 10000, replace = FALSE)
#species_astrocyte_indices <- grep("mm10", names(astrocytePredictions))
species_astrocyte_indices <- grep("hg38", names(astrocytePredictions))
subset_astrocyte_indices <- sample(species_astrocyte_indices, size = 10000, replace = FALSE)
astrocytePredictions <- astrocytePredictions[, subset_astrocyte_indices, drop = FALSE]
astrocytePredictions[astrocytePredictions == -1] <- NA

all_vip_indices <- 1:ncol(vipPredictions)
#non_common_vip_indices <- setdiff(all_vip_indices, common_columns)
#subset_vip_indices <- sample(non_common_vip_indices, size = 10000, replace = FALSE)
#species_vip_indices <- grep("mm10", names(vipPredictions))
species_vip_indices <- grep("hg38", names(vipPredictions))
subset_vip_indices <- sample(species_vip_indices, size = 10000, replace = FALSE)
vipPredictions <- vipPredictions[, subset_vip_indices, drop = FALSE]
vipPredictions[vipPredictions == -1] <- NA




#Load in original data to subset it --> added csv file to home directory of bridges ~___
originalZoonomiaData <- read.csv("Phenotypes_Zoonomia_3.csv", header = TRUE)
sortedZoonomia <- originalZoonomiaData[order(originalZoonomiaData$Name), ]
#Remove rows where Homeotherm value is NA
sortedZoonomia <- sortedZoonomia[!is.na(sortedZoonomia$Homeotherm), ]
# Step 1: Get the species names from sortedZoonomia
zoonomia_species <- sortedZoonomia$Name


astrocytePredictions <- astrocytePredictions[rownames(astrocytePredictions) %in% zoonomia_species, ]
astrocytePredictionsNames <- astrocytePredictionsNames[astrocytePredictionsNames$V1 %in% zoonomia_species, ]
#Need to convert the names back into a df
astrocytePredictionsNames <- data.frame(Species = astrocytePredictionsNames, stringsAsFactors = FALSE)
#Add homeotherm column
astrocytePredictionsNames$Homeotherm <- sortedZoonomia$Homeotherm

vipPredictions <- vipPredictions[rownames(vipPredictions) %in% zoonomia_species, ]
vipPredictionsNames <- vipPredictionsNames[vipPredictionsNames$V1 %in% zoonomia_species, ]
#Need to convert the names back into a df
vipPredictionsNames <- data.frame(Species = vipPredictionsNames, stringsAsFactors = FALSE)
#Add homeotherm column
vipPredictionsNames$Homeotherm <- sortedZoonomia$Homeotherm




# Tree info - added tree file to ~ path
# Read the tree
tree <- read.tree("Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")
total_species <- length(tree$tip.label)
cat("Total number of species in the original tree:", total_species, "\n")
# Read the phenotypes data
pheno_df <- read.csv("Phenotypes_Zoonomia_3.csv")
# Create a subset of pheno_df to include only species with Homeotherm values
pheno_df_subset <- pheno_df[!is.na(pheno_df$Homeotherm), c("Homeotherm", "Species")]
rownames(pheno_df_subset) <- pheno_df_subset$Species
pheno_df_subset$Homeotherm <- factor(pheno_df_subset$Homeotherm)
# Trim whitespace from species names
tree$tip.label <- trimws(tree$tip.label)
pheno_df_subset$Species <- trimws(pheno_df_subset$Species)
# Identify species in pheno_df_subset that are not in the tree
species_not_in_tree <- setdiff(pheno_df_subset$Species, tree$tip.label)
# Remove species not present in the tree
remove_species <- setdiff(tree$tip.label, pheno_df_subset$Species)
pruned_tree <- drop.tip(tree, remove_species)
pruned_tree$tip.label <- sort(pruned_tree$tip.label)
# Print the number of species in the pruned tree
num_species_pruned <- length(pruned_tree$tip.label)
cat("Number of species in the pruned tree:", num_species_pruned, "\n")
# Prepare tip colors based on Homeotherm factor
useTipColorsV <- c("red", "blue")
tipColorsV <- useTipColorsV[pheno_df_subset[pruned_tree$tip.label, "Homeotherm"]]
# Plot the pruned tree
# plot(pruned_tree, type = "fan", cex = 0.3, tip.color = tipColorsV)



vector_to_remove <- c('Balaenoptera_acutorostrata', 'Bison_bison', 'Cebus_capucinus', 'Colobus_angolensis', 'Dicerorhinus_sumatrensis', 'Equus_asinus',
                     'Gorilla_gorilla', 'Marmota_marmota', 'Neophocaena_asiaeorientalis', 'Oryctolagus_cuniculus', 'Ovis_canadensis', 'Panthera_tigris',
                     'Perognathus_longimembris', 'Peromyscus_maniculatus', 'Saimiri_boliviensis')
astrocytePredictionsNames = astrocytePredictionsNames[-(which(astrocytePredictionsNames$Species %in% vector_to_remove)),]
astrocytePredictions = astrocytePredictions[-(which(rownames(astrocytePredictions) %in% vector_to_remove)),]
vipPredictionsNames = vipPredictionsNames[-(which(vipPredictionsNames$Species %in% vector_to_remove)),]
vipPredictions = vipPredictions[-(which(rownames(vipPredictions) %in% vector_to_remove)),]



print("Astrocyte phylolm")
print(length(colnames(astrocytePredictions)))
testPeaksV <- colnames(astrocytePredictions)
loopPeaksV <- testPeaksV[1:length(colnames(astrocytePredictions))]
peakPhyloResultsF <- data.frame(peakId = loopPeaksV, pvalue=rep(NA,length(loopPeaksV)),correlation=rep(NA,length(loopPeaksV)),adjCorrelation=rep(NA,length(loopPeaksV)))
rownames(peakPhyloResultsF) <- peakPhyloResultsF$peakId
options(warn=-1)
for(curPeak in loopPeaksV) {
    speciesDetailedInfoTmpF <- astrocytePredictionsNames
    matching_order <- match(pruned_tree$tip.label, speciesDetailedInfoTmpF$Species)
    
    speciesDetailedInfoTmpF <- speciesDetailedInfoTmpF[matching_order, ]
    astrocytePredictions <- astrocytePredictions[matching_order, ]
    
    speciesDetailedInfoTmpF$curPeak <- astrocytePredictions[, curPeak]
    speciesDetailedInfoTmpF$curPeak <- as.numeric(speciesDetailedInfoTmpF$curPeak)
    
    
    rownames(speciesDetailedInfoTmpF) <- NULL
    rownames(speciesDetailedInfoTmpF) = speciesDetailedInfoTmpF$Species

    # Perform the phylogenetic regression
    #curLmFit = phylolm(curPeak ~ Homeotherm, data=speciesDetailedInfoTmpF, phy=pruned_tree)
    curLmFit = phyloglm(Homeotherm ~ curPeak, data=speciesDetailedInfoTmpF, phy=pruned_tree)
    curLmFitSum <- summary(curLmFit)
    #peakPhyloResultsF[curPeak, "pvalue"] <- curLmFitSum$coefficients["Homeotherm", "p.value"]
    peakPhyloResultsF[curPeak, "pvalue"] <- curLmFitSum$coefficients["curPeak", "p.value"]
   
    #peakPhyloResultsF[curPeak, "correlation"] <- curLmFitSum$r.squared
    #peakPhyloResultsF[curPeak, "adjCorrelation"] <- curLmFitSum$adj.r.squared
    #print("end of peak")
}
options(warn=0)
file_name <- sprintf("Peak_Phylolm/peakPhyloResults_astrocyte_human_new2.csv")
write.csv(peakPhyloResultsF, file = file_name, row.names = FALSE)
flush.console()
Sys.sleep(1)




print("VIP phylolm")
print(length(colnames(vipPredictions)))
testPeaksV <- colnames(vipPredictions)
loopPeaksV <- testPeaksV[1:length(colnames(vipPredictions))]
peakPhyloResultsF <- data.frame(peakId = loopPeaksV, pvalue=rep(NA,length(loopPeaksV)),correlation=rep(NA,length(loopPeaksV)),adjCorrelation=rep(NA,length(loopPeaksV)))
rownames(peakPhyloResultsF) <- peakPhyloResultsF$peakId
options(warn=-1)
for(curPeak in loopPeaksV) {
    speciesDetailedInfoTmpF <- vipPredictionsNames
    matching_order <- match(pruned_tree$tip.label, speciesDetailedInfoTmpF$Species)
    
    speciesDetailedInfoTmpF <- speciesDetailedInfoTmpF[matching_order, ]
    vipPredictions <- vipPredictions[matching_order, ]
    
    speciesDetailedInfoTmpF$curPeak <- vipPredictions[, curPeak]
    speciesDetailedInfoTmpF$curPeak <- as.numeric(speciesDetailedInfoTmpF$curPeak)
  
    rownames(speciesDetailedInfoTmpF) <- NULL
    rownames(speciesDetailedInfoTmpF) = speciesDetailedInfoTmpF$Species

    # Perform the phylogenetic regression
    curLmFit = phyloglm(Homeotherm ~ curPeak, data=speciesDetailedInfoTmpF, phy=pruned_tree)
    curLmFitSum <- summary(curLmFit)
    
    peakPhyloResultsF[curPeak, "pvalue"] <- curLmFitSum$coefficients["curPeak", "p.value"]
    #peakPhyloResultsF[curPeak, "correlation"] <- curLmFitSum$r.squared
    #peakPhyloResultsF[curPeak, "adjCorrelation"] <- curLmFitSum$adj.r.squared
}
options(warn=0)
file_name <- sprintf("Peak_Phylolm/peakPhyloResults_vip_human_new2.csv")
write.csv(peakPhyloResultsF, file = file_name, row.names = FALSE)
flush.console()
Sys.sleep(1)
