# Noma Metagenomics: AMR - Figure 5
# By: Angus M. O'Ferrall
# Shotgun metagenomic analysis of the oral microbiomes of children with noma reveals a novel disease-associated organism
# Michael Olaleye, Angus M O’Ferrall, Richard N Goodman, Deogracia W Kabila, Miriam Peters, Gregoire Falq, Joseph Samuel, Donal Doyle, Diana Gomez, Gbemisola Oloruntuyi, Shafiu Isah, Adeniyi S Adetunji, Elise N Farley, Nicholas J Evans, Mark Sherlock, Adam P Roberts, Mohana Amirtharajah, Stuart Ainsworth

# Load packages
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

## Prepare the antibiotic class data Define the colour palette for antibiotic classes in plots

# Define the AMR classes and their associated genes
aminoglycoside = "Aac|Aad|Aph"
beta_lactam = "Cfx|CTX|OXA|Pen|SPU|TEM"
chloramphenicol = "Cat"
colistin = "MCR"
lincosamide = "Lnu|Lsa"
macrolide = "Mef|Msr"
MDR = "Erm"
nitroimidazole = "Nim"
quinolone = "Qep"
sulphonamide = "Sul"
tetracycline = "Tet"
trimethoprim = "Dfr"

# Define a consistent color palette for all antibiotic classes
class_colors <- c(
  "Aminoglycoside" = "#66c2a5",
  "Beta-lactam" = "#fc8d62",
  "Chloramphenicol" = "#8da0cb",
  "Colistin" = "#e78ac3",
  "Lincosamide" = "#a6d854",
  "Macrolide" = "#ffd92f",
  "MDR" = "#e5c494",
  "Nitroimidazole" = "#b3b3b3",
  "Quinolone" = "#ff9999",
  "Sulphonamide" = "#c2e699",
  "Tetracycline" = "#fdb462",
  "Trimethoprim" = "#80b1d3"
)

## Prepare the noma oral metagenome AMR data

# Import noma ARG data
ARGs <- read.csv("../data/ARGs_noma.csv", check.names = FALSE)

# Convert to present (1) or absent (0)
ARGs[, -1] <- apply(ARGs[, -1], 2, function(x) {
  ifelse(is.na(as.numeric(x)), 0, 1)
})

# Remove any columns in which there are no ARGs after transformations
ARGs <- ARGs[, colSums(ARGs != 0) > 0]

# Mutate the data to have Participant and ARG columns
ARGs_long <- ARGs %>%
  pivot_longer(
    cols = -Participant,  # All columns except 'Participant'
    names_to = "ARG",  # New column name for original column names
    values_to = "Presence"  # New column name for values (0 or 1)
  )

# filter for only relevant rows (where presence is 1)
ARGs_long <- ARGs_long %>%
  filter(Presence == 1) %>%  # Keep only rows where ARG is present (1)
  select(-Presence)          # Remove the Presence column

# print all unique AMR determinants in data set
(colnames(ARGs)[-1])

# Calculate the absolute abundance (counting only where Presence is 1) and label ARG classes
ARG_abundance <- ARGs_long %>%
  group_by(ARG) %>%                   # Group by ARG names
  summarise(Count = n(), .groups = "drop") %>%         # Count occurrences of each ARG
  # Label ARG classes for each determinant
  mutate(class = case_when(
    str_detect(ARG, aminoglycoside) ~ "Aminoglycoside",
    str_detect(ARG, beta_lactam) ~ "Beta-lactam",
    str_detect(ARG, chloramphenicol) ~ "Chloramphenicol",
    str_detect(ARG, lincosamide) ~ "Lincosamide",
    str_detect(ARG, macrolide) ~ "Macrolide",
    str_detect(ARG, nitroimidazole) ~ "Nitroimidazole",
    str_detect(ARG, quinolone) ~ "Quinolone",
    str_detect(ARG, sulphonamide) ~ "Sulphonamide",
    str_detect(ARG, tetracycline) ~ "Tetracycline",
    str_detect(ARG, trimethoprim) ~ "Trimethoprim",
    TRUE ~ NA_character_  # Set to NA for unmatched ARGs
  ))

# Add percentage (proportion)
ARG_abundance$Percentage <- (ARG_abundance$Count/19)*100

# Reorder ARGs for plotting: first by class, then by count (highest to lowest)
ARG_abundance <- ARG_abundance %>%
  arrange(class, desc(Count)) %>%  # Sort by class and then by count
  mutate(ARG = factor(ARG, levels = ARG))  # Reorder ARG factor levels

# Reverse the order of classes and ARGs for top-to-bottom plotting
ARG_abundance <- ARG_abundance %>%
  mutate(class = factor(class, levels = sort(unique(class)))) %>%
  mutate(ARG = factor(ARG, levels = rev(ARG[order(ARG_abundance$class, -ARG_abundance$Count)])))  # Reverse the order of ARGs within each class

## Prepare the global healthy oral metagenome AMR data

# Import healthy ARG data
ARGs_ctrls <- read.csv("../data/ARGs_ctrls.csv", check.names = FALSE)

# Convert to present (1) or absent (0)
ARGs_ctrls[, -1] <- apply(ARGs_ctrls[, -1], 2, function(x) {
  ifelse(is.na(as.numeric(x)), 0, 1)
})

# Remove any columns in which there are no ARGs after transformations
ARGs_ctrls <- ARGs_ctrls[, colSums(ARGs_ctrls != 0) > 0]

# Mutate the data to have Participant and ARG columns
ARGs_ctrls_long <- ARGs_ctrls %>%
  pivot_longer(
    cols = -Sample,  # All columns except 'Sample'
    names_to = "ARG",  # New column name for original column names
    values_to = "Presence"  # New column name for values (0 or 1)
  )

# filter for only relevant rows (where presence is 1)
ARGs_ctrls_long <- ARGs_ctrls_long %>%
  filter(Presence == 1) %>%  # Keep only rows where ARG is present (1)
  select(-Presence)          # Remove the Presence column

# print all unique AMR determinants in data set
(colnames(ARGs_ctrls)[-1])

# Calculate the absolute abundance (counting only where Presence is 1) and label ARG classes
ARG_ctrls_abundance <- ARGs_ctrls_long %>%
  group_by(ARG) %>%                   # Group by ARG names
  summarise(Count = n(), .groups = "drop") %>%         # Count occurrences of each ARG
  # Label ARG classes for each determinant
  mutate(class = case_when(
    str_detect(ARG, aminoglycoside) ~ "Aminoglycoside",
    str_detect(ARG, beta_lactam) ~ "Beta-lactam",
    str_detect(ARG, chloramphenicol) ~ "Chloramphenicol",
    str_detect(ARG, colistin) ~ "Colistin",
    str_detect(ARG, lincosamide) ~ "Lincosamide",
    str_detect(ARG, macrolide) ~ "Macrolide",
    str_detect(ARG, MDR) ~ "MDR",
    str_detect(ARG, sulphonamide) ~ "Sulphonamide",
    str_detect(ARG, tetracycline) ~ "Tetracycline",
    str_detect(ARG, trimethoprim) ~ "Trimethoprim",
    TRUE ~ NA_character_  # Set to NA for unmatched ARGs
  ))

# Add percentage (proportion)
ARG_ctrls_abundance$Percentage <- (ARG_ctrls_abundance$Count/20)*100

# Reorder ARGs for plotting: first by class, then by count (highest to lowest)
ARG_ctrls_abundance <- ARG_ctrls_abundance %>%
  arrange(class, desc(Count)) %>%  # Sort by class and then by count
  mutate(ARG = factor(ARG, levels = ARG))  # Reorder ARG factor levels

# Reverse the order of classes and ARGs for top-to-bottom plotting
ARG_ctrls_abundance <- ARG_ctrls_abundance %>%
  mutate(class = factor(class, levels = sort(unique(class)))) %>%
  mutate(ARG = factor(ARG, levels = rev(ARG[order(ARG_ctrls_abundance$class, -ARG_ctrls_abundance$Count)])))  # Reverse the order of ARGs within each class

## Generate metagenome plots

# Plot for Noma oral metagenome AMR (no legend)
AMR_plot <- ggplot(ARG_abundance, aes(y = ARG, x = Percentage, fill = class)) +
  geom_bar(stat = "identity") +
  labs(y = "AMR Determinant", x = "Percentage of oral microbiomes", title = "Noma dataset") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "none"  # Remove legend
  ) +
  scale_x_continuous(limits = c(0, 100)) +  # Ensure x-axis limit is 100
  scale_fill_manual(values = class_colors, drop = FALSE)

# Plot for Global healthy oral metagenome AMR (no legend)
AMR_ctrls_plot <- ggplot(ARG_ctrls_abundance, aes(y = ARG, x = Percentage, fill = class)) +
  geom_bar(stat = "identity") +
  labs(y = "AMR Determinant", x = "Percentage of oral microbiomes", title = "Global healthy dataset") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "none"  # Remove legend
  ) +
  scale_x_continuous(limits = c(0, 100)) +  # Ensure x-axis limit is 100
  scale_fill_manual(values = class_colors, drop = FALSE)

# Create a dummy plot for the legend only including all antibiotic classes in noma and healthy samples
legend_plot <- ggplot(data.frame(class = names(class_colors)), aes(x = 1, y = 1, fill = class)) +
  geom_bar(stat = "identity") +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_manual(values = class_colors, name = "Class", drop = FALSE)  # Set legend title to "Class"

# Extract only the legend from the dummy plot
legend <- cowplot::get_legend(legend_plot)

# Combine plots side by side with the legend to the right
combined_plot <- (AMR_plot | AMR_ctrls_plot) + 
  patchwork::plot_layout(guides = "collect") & theme(legend.position = "none")

# Add the legend to the combined plot
final_plot <- cowplot::plot_grid(combined_plot, legend, rel_widths = c(4, 1))

# Save the updated plot
ggsave("../imgs/Figure_5A.png", final_plot, height = 7, width = 10)

## Noma MAG AMR

# Import MAG-associated ARG data
MAG_AMR <- read.csv("../data/MAG_AMR.csv", check.names = FALSE)

# Reverse the order of AMR determinant variable for A-> Z plotting
MAG_AMR <- MAG_AMR %>%
  mutate(AMR_determinant = factor(AMR_determinant, levels = rev(sort(unique(AMR_determinant)))))

# Plot the AMR determinants in mid & high quality MAGs
MAG_AMR_plot <- ggplot(mapping = aes(x = AMR_determinant, fill = Genus), data = MAG_AMR) +
  geom_bar() +
  theme_minimal() +
  labs(x = "AMR Determinant", y = "Number of MAGs", title = "Noma MAGs") +
  coord_flip() +
  scale_fill_discrete(labels = function(x) {
    sapply(x, function(genus) {
      if (str_detect(genus, "\\(unclassified\\)$")) {
        genus_parts <- str_split_fixed(genus, " \\(unclassified\\)", 2)
        bquote(italic(.(genus_parts[1])) ~ "(unclassified)")
      } else {
        bquote(italic(.(genus)))
      }
    })
  })

# Save the plot
ggsave("../imgs/Figure_5B.png", plot = MAG_AMR_plot, width = 4.5, height = 4, dpi = 300)
