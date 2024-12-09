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


astrocytePredictionsPath <- "240_predictions_MatrixStacked.tsv"
astrocytePredictions <- read.delim(astrocytePredictionsPath, header = FALSE, sep = "\t")


astrocytePredictionsNamesPath <- "240_predictions_NamesList.txt"
astrocytePredictionsNames <- read.delim(astrocytePredictionsNamesPath, header = FALSE, sep = "\t")

# Extract the species names from astrocytePredictionsNames and make them column names
species_names <- astrocytePredictionsNames[[1]]

# Check if the length matches the number of columns (excluding the first)
if (length(species_names) == (ncol(astrocytePredictions) - 1)) {

  # Rename the columns in astrocytePredictions (excluding the first column)
  colnames(astrocytePredictions)[2:ncol(astrocytePredictions)] <- species_names

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

#Load in original data to subset it --> added csv file to home directory of bridges ~___
originalZoonomiaData <- read.csv("Phenotypes_Zoonomia_3.csv", header = TRUE)
sortedZoonomia <- originalZoonomiaData[order(originalZoonomiaData$Name), ]

#Remove rows where Homeotherm value is NA
sortedZoonomia <- sortedZoonomia[!is.na(sortedZoonomia$Homeotherm), ]


# Step 1: Get the species names from sortedZoonomia
zoonomia_species <- sortedZoonomia$Name
#zoonomia_species

astrocytePredictions <- astrocytePredictions[rownames(astrocytePredictions) %in% zoonomia_species, ]
astrocytePredictionsNames <- astrocytePredictionsNames[astrocytePredictionsNames$V1 %in% zoonomia_species, ]

#Need to convert the names back into a df
astrocytePredictionsNames <- data.frame(Species = astrocytePredictionsNames, stringsAsFactors = FALSE)


#Add homeotherm column
astrocytePredictionsNames$Homeotherm <- sortedZoonomia$Homeotherm



# Tree info - added tree file to ~ path


# Read the tree
tree <- read.tree("Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")
total_species <- length(tree$tip.label)


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

# Prepare tip colors based on Homeotherm factor
useTipColorsV <- c("red", "blue")
tipColorsV <- useTipColorsV[pheno_df_subset[pruned_tree$tip.label, "Homeotherm"]]


vector_to_remove <- c('Balaenoptera_acutorostrata', 'Bison_bison', 'Cebus_capucinus', 'Colobus_angolensis', 'Dicerorhinus_sumatrensis', 'Equus_asinus',
                     'Gorilla_gorilla', 'Marmota_marmota', 'Neophocaena_asiaeorientalis', 'Oryctolagus_cuniculus', 'Ovis_canadensis', 'Panthera_tigris',
                     'Perognathus_longimembris', 'Peromyscus_maniculatus', 'Saimiri_boliviensis')

astrocytePredictionsNames = astrocytePredictionsNames[-(which(astrocytePredictionsNames$Species %in% vector_to_remove)),]
astrocytePredictions = astrocytePredictions[-(which(rownames(astrocytePredictions) %in% vector_to_remove)),]

slice_indices <- c(10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,110000,120000,130000,140000,150000,160000,170000,180000,190000,200000,205151)
for (i in seq_along(slice_indices)){
    testPeaksV <- colnames(astrocytePredictions)
    start_index <- ifelse(i == 1, 1, slice_indices[i-1] + 1)
    loopPeaksV <- testPeaksV[start_index:slice_indices[i]]
    print(start_index)
    print(slice_indices[i])
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
        curLmFit = phylolm(curPeak ~ Homeotherm, data=speciesDetailedInfoTmpF, phy=pruned_tree, model="BM")
        curLmFitSum <- summary(curLmFit)

        peakPhyloResultsF[curPeak, "pvalue"] <- curLmFitSum$coefficients["Homeotherm", "p.value"]
        peakPhyloResultsF[curPeak, "correlation"] <- curLmFitSum$r.squared
        peakPhyloResultsF[curPeak, "adjCorrelation"] <- curLmFitSum$adj.r.squared
    }
    options(warn=0)
    file_name <- sprintf("Peak_Phylolm/peakPhyloResults_%d_%d.csv", start_index, slice_indices[i])
    write.csv(peakPhyloResultsF, file = file_name, row.names = FALSE)
    flush.console()
    Sys.sleep(1)
}



