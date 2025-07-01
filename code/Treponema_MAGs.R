# Noma Metagenomics: Recovery and analysis of Treponema MAGs from noma samples - Figure 4A
# By: Angus M. O'Ferrall
# Shotgun metagenomic analysis of the oral microbiomes of children with noma reveals a novel disease-associated organism
# Michael Olaleye, Angus M O’Ferrall, Richard N Goodman, Deogracia W Kabila, Miriam Peters, Gregoire Falq, Joseph Samuel, Donal Doyle, Diana Gomez, Gbemisola Oloruntuyi, Shafiu Isah, Adeniyi S Adetunji, Elise N Farley, Nicholas J Evans, Mark Sherlock, Adam P Roberts, Mohana Amirtharajah, Stuart Ainsworth

# Load packages
library(tidyverse)
library(ggplot2)
library(ggtext)

# Import MAG data
MAG_data <- read.csv("../data/Treponema_MAG_data.csv")

# Preprocess metadata: italicize species selectively
italicize_terms <- function(species_name) {
  # Italicize specific species names
  species_name <- gsub("Spirochaeta thermophila", "*Spirochaeta thermophila*", species_name, ignore.case = TRUE)
  species_name <- gsub("Treponema maltophilum", "*Treponema maltophilum*", species_name, ignore.case = TRUE)
  species_name <- gsub("Treponema denticola", "*Treponema denticola*", species_name, ignore.case = TRUE)
  species_name <- gsub("Treponema medium", "*Treponema medium*", species_name, ignore.case = TRUE)
  species_name <- gsub("Treponema parvum", "*Treponema parvum*", species_name, ignore.case = TRUE)
  species_name <- gsub("Treponema lecithinolyticum", "*Treponema lecithinolyticum*", species_name, ignore.case = TRUE)
  species_name <- gsub("Treponema socranskii", "*Treponema socranskii*", species_name, ignore.case = TRUE)
  species_name <- gsub("Treponema sp905372235", "*Treponema sp905372235*", species_name, ignore.case = TRUE)
  species_name <- gsub("Treponema sp014334325", "*Treponema sp014334325*", species_name, ignore.case = TRUE)
  
  # If "sp." is in the species name, italicize "Treponema" and bold the whole name (for novel species)
  if (grepl("sp\\.", species_name, ignore.case = TRUE)) {
    # Italicize "Treponema" and add bold formatting for the whole name
    species_name <- gsub("Treponema", "*Treponema*", species_name, ignore.case = TRUE)
    species_name <- paste0("**", species_name, "**")  # Wrap the entire name in bold
  }
  
  return(species_name)
}

# Apply the italicization function to the metadata
MAG_data <- MAG_data %>%
  mutate(Species = sapply(Species, italicize_terms))  # Italicize species selectively

# Plot the species of medium & high quality Treponema MAGs
Treponema_MAG_plot <- ggplot(mapping = aes(x = fct_rev(fct_infreq(Species)), fill = Quality), data = MAG_data) +
  geom_bar() +
  theme_minimal() +
  labs(x = "Species", y = "Count", fill = "MAG quality") +
  coord_flip() +
  theme(axis.text.y = element_markdown()) +  # Enable Markdown rendering for y-axis labels
  scale_fill_manual(
    values = c(
      "High" = "#CCCCCC",  # Light Gray
      "Medium" = "#666666" # Dark Gray
    )
  )

# Save the plot
ggsave("../imgs/Figure_4A.png", plot = Treponema_MAG_plot, width = 5, height = 4, dpi = 600)
