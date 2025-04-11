# Load necessary libraries ----
library(ape)
# install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("BiocGenerics")
library(Biostrings)
library(DECIPHER)
library(flextable)
# BiocManager::install("ggtree")
library(ggtree)
library(ggiraph)
library(phangorn)
library(rentrez)
library(tidyverse)
library(tidytree)

# Build Table ----
tab1 <- data.frame(
        Code = c("μ01","μ03","μ05","μ17","μ23","μ33","μ07","μ13","μ19","μ25",
                 "μ45","μ09","μ27","μ39","μ11","μ29","μ47","μ21","μ31","μ35",
                 "μ41","μ15","μ37","μ43","μ49"),
        Culture = c("Reference","Reference","Reference","Farmhouse Ale","Farmhouse Ale",
                    "Farmhouse Ale","Sour Ale","Sour Ale","Sour Ale","Sour Ale",
                    "Sour Ale","Tibicos","Tibicos","Tibicos","Kombucha","Kombucha",
                    "Kombucha","Kefir A","Kefir A","Kefir A",
                    "Kefir A","Kefir B","Kefir B","Kefir B","Kefir B"),
        Species = c("Propionibacterium freudenreichii","Saccharomyces cerevisiae","Lactococcus lactis",
                    "Saccharomyces cerevisiae","Wickerhamomyces anomalus","Candida pseudolambica",
                    "Lanchancea fermentati","Saccharomyces cerevisiae","Levilactobacillus brevis",
                    "Saccharomyces cerevisiae","Pichia occidentalis","Pichia membranifaciens",
                    "Saccharomyces cerevisiae","Acetobacter fabarum","Saccharomyces cerevisiae",
                    "Torulaspora delbrueckii","Streptococcus salivarius*","Saccharomyces cerevisiae",
                    "Saccharomycodes ludwigii","Levilactobacillus brevis","Lacticaseibacillus paracasei", 
                    "Debaryomyces fabryi","Pichia membranifaciens","Levilactobacillus brevis",
                    "Lacticaseibacillus paracasei"),
        GenBankID = c("—","—","—","ON763764","ON763766","ON763771","ON763759",
                      "ON763762","ON758946","ON763767","ON763773","ON763760",
                      "ON763768","ON758948","ON763761","ON763769","—",
                      "ON763765","ON763770","ON758947","ON758949","ON763763",
                      "ON763772","ON758950","ON758951"),
        stringsAsFactors = FALSE
)

saveRDS(tab1, "rds_objects/tab1.rds")

# Build Tree ----
# Filter out rows with missing GenBank IDs
filt_data <-  dplyr::filter(tab1, !is.na(GenBankID))

# Extract the list of GenBank IDs
ids <- filt_data$GenBankID
ids <- ids[ids != "—"]

# Extract the codes
codes <- filt_data$Code[filt_data$GenBankID %in% ids]

# Fetch sequences or load cached
if (!file.exists("rds_objects/fasta_seqs.rds")) {
        seqs_raw <- lapply(ids, function(id){
                entrez_fetch(db = "nuccore", id = id, rettype = "fasta", retmode = "text")
        })
        saveRDS(seqs_raw, "rds_objects/fasta_seqs.rds")
} else {
        seqs_raw <- readRDS("rds_objects/fasta_seqs.rds")
}

# Split the input into lines
seq_lines <- lapply(seqs_raw, function(seq){
        strsplit(seq, "\n")[[seq_along(seq)]]
})

# Extract the sequences
seqs <- lapply(seq_lines, function(line){
        paste(line[-1], collapse = "")
})

# Create the DNAStringSet objects
dna_set <- DNAStringSet(unlist(seqs))

# Set the name of the DNAStringSet object
names(dna_set) <- codes

# Align the sequences using DECIPHER
alignment <- AlignSeqs(dna_set, verbose = FALSE)

# Transform to a PhyDat object
phang_alignment <- phyDat(as.matrix(alignment), type = "DNA")

# Calculate the distance matrix 
dm <- dist.ml(phang_alignment, model = "JC69")

# Construct the phylogenetic tree using neighbor-joining method
tree <- NJ(dm)

# Using Maximum Likelihood
ml_tree <- pml(tree, phang_alignment)

# Optimization of the phylogenetic likelihood model
ml_tree <- optim.pml(
        ml_tree,
        model = "GTR",
        optNni = TRUE,
        control=pml.control(trace=0))

# Root the tree
eubacteria <- c("μ19", "μ35", "μ39", "μ41", "μ43", "μ49")
saccharomyces <- c("μ11", "μ13", "μ17", "μ21", "μ25", "μ27")
fungi <- codes[!codes %in% eubacteria]

ml_rooted_tree <- root(ml_tree$tree, outgroup = eubacteria, resolve.root = TRUE)

eubacteria_mrca <- getMRCA(ml_rooted_tree, eubacteria)
saccharomyces_mrca <- getMRCA(ml_rooted_tree, saccharomyces)
fungi_mrca <- getMRCA(ml_rooted_tree, fungi)

# Plot the phylogenetic tree ----
ggtree(ml_rooted_tree, layout='circular', branch.length='none') + 
        geom_tiplab(offset = 0.2) +
        # geom_text(aes(label = node), hjust = -0.5, size = 3, color = "red") +
        geom_hilight(node = eubacteria_mrca, 
                     fill = "darkgreen", 
                     alpha = 0.5,
                     extend = 0.2) +
        geom_cladelab(node = eubacteria_mrca,
                      label = "Eubacteria",
                      angle = 337,
                      fontsize = 5,
                      offset = 4,
                      vjust = -0.5,
                      hjust = 0.5,
                      horizontal = FALSE,
                      barsize = 1,
                      barcolour = "darkgreen",
                      textcolour = "darkgreen") +
        geom_hilight(node = saccharomyces_mrca,
                     fill = "navy",
                     alpha = 0.75,
                     extend = 0.2) +
        geom_cladelab(node = saccharomyces_mrca,
                      label = "Saccharomyces",
                      angle = 43,
                      fontsize = 5,
                      offset = 4,
                      vjust = 1.5,
                      hjust = 0.35,
                      horizontal = FALSE,
                      barsize = 2,
                      barcolour = "navy",
                      textcolour = "navy") +
        geom_hilight(node = fungi_mrca,
                     fill = "lightgray",
                     alpha = 0.25,
                     extend = 0.2) +
        geom_cladelab(node = fungi_mrca,
                      label = "Non-Saccharomyces",
                      angle = 310,
                      fontsize = 5,
                      offset = 4,
                      vjust = 1.5,
                      hjust = 0.8,
                      horizontal = FALSE,
                      barsize = 1,
                      barcolour = "darkgray",
                      textcolour = "darkgray")
ggsave(filename = "images/cladogram.png",
       height = 9, width = 12, units = "cm", device = "png", dpi = 300)
