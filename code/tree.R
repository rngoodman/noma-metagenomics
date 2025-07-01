# Noma Metagenomics: Recovery and analysis of Treponema MAGs from noma samples - Figure 4B
# By: Angus M. O'Ferrall
# Shotgun metagenomic analysis of the oral microbiomes of children with noma reveals a novel disease-associated organism
# Michael Olaleye, Angus M O’Ferrall, Richard N Goodman, Deogracia W Kabila, Miriam Peters, Gregoire Falq, Joseph Samuel, Donal Doyle, Diana Gomez, Gbemisola Oloruntuyi, Shafiu Isah, Adeniyi S Adetunji, Elise N Farley, Nicholas J Evans, Mark Sherlock, Adam P Roberts, Mohana Amirtharajah, Stuart Ainsworth

# Load packages
library(tidyverse)
library(ggtree)
library(RColorBrewer)
library(ggtext)

# Import tree file and metadata
tree <- read.tree("../data/core_genome_snp_sites.aln.treefile")
tree_metadata <- read.csv("../data/Tree_metadata.csv")

# Rename genome column to match the node tips/labels
tree_metadata <- tree_metadata %>% 
  rename(label = ID)

# Join the metadata with the tree tips
tree_data <- full_join(as_tibble(tree), tree_metadata, by = "label")

# Preprocess metadata: Italicize species selectively
italicize_terms <- function(species_name) {
  # List of full species names to fully italicize
  species_list <- c(
    "Treponema lecithinolyticum",
    "Treponema sp905372235",
    "Treponema denticola",
    "Spirochaeta thermophila",
    "Treponema succinifaciens",
    "Treponema brennaborense",
    "Treponema vincentii",
    "Treponema maltophilum",
    "Treponema pallidum",
    "Treponema pedis",
    "Treponema phagedenis",
    "Treponema ruminis",
    "Treponema peruense",
    "Treponema medium",
    "Treponema parvum",
    "Treponema socranskii",
    "Treponema putidum"
  )
  
  # Italicize the full species names
  for (species in species_list) {
    species_name <- gsub(species, paste0("*", species, "*"), species_name, ignore.case = TRUE)
  }
  
  ## If "sp." is in the species name, italicize "Treponema" and bold the whole name (for novel species)
  if (grepl("sp\\.", species_name, ignore.case = TRUE)) {
    # Italicize "Treponema" and add bold formatting for the whole name
    species_name <- gsub("Treponema", "*Treponema*", species_name, ignore.case = TRUE)
    species_name <- paste0("**", species_name, "**")  # Wrap the entire name in bold
  }
  
  return(species_name)
}


# Apply the italicization function to the metadata
tree_metadata <- tree_metadata %>%
  mutate(Species = sapply(Species, italicize_terms))  # Italicize species selectively

# Create a distinct color palette using RColorBrewer
distinct_colors <- brewer.pal(12, "Set3")
if (nrow(tree_metadata) > 12) {
  # Scale colors to more species dynamically if needed
  distinct_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nrow(tree_metadata))
}

# Create the tree visualization
tree_plot <- ggtree(tree) %<+% tree_metadata + 
  geom_tippoint(aes(color = Species), size = 2.5) + 
  geom_tiplab(aes(label = Genome), size = 3, align = TRUE, offset = 0.01) +
  geom_treescale(x = 0, y = 33, width = 0.2, fontsize = 3, linesize = 0.5) + # Add scale bar
  theme_tree2() +
  scale_color_manual(values = distinct_colors, guide = guide_legend(ncol = 1)) +  # Force legend into 1 column
  theme(
    legend.text = element_markdown(),
    axis.line = element_blank(),     # Remove the axis line
    axis.ticks = element_blank(),    # Remove the axis ticks
    axis.text = element_blank(),     # Remove the axis text
    axis.title = element_blank()     # Remove axis titles
  ) +
  ggplot2::xlim(0, 2.65)  # Adjust limits

# Save the plot
ggsave("../imgs/Figure_4B.png", plot = tree_plot, width = 8, height = 6.5, dpi = 600)
