###BIOL503_Lab8_2025-03-17###

## load libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(pairwiseAdonis)
library(gridExtra)
library(dplyr)
library(ggplot2)

## set working directory
setwd("C:/Users/allyw/Desktop/UBC_Master_Zoology/Courses_2025/Winter_Term_2025/BIOL503_MicrobialEcology/BIOL503_Project/")

## Read in data
seagrassrare<-readRDS("rarefied_seagrass_default_phyloseq_2025_03_10.RDS")

###Hypothesis: Recovery of the rhizosphere post-transplantation will be different between sod & wash treatments due to different magnitudes of disturbances in microbial community composition during transplanting.

#Set.seed so that NDMS plots are reproducible
set.seed(4)

#create data frame speficying days for the for loop
days <- c("d0", "d1", "d3", "d7", "d14", "d21", "d28")

#create a list where plots, the stress values, and the permanova outputs from the for-loop can be saved
plots<- list()
stress.day<- vector()
perm.assum<-list()
permanova<- list()
p_value<- vector()

for (day in days) {
  # Subset samples for the current day
  d <- subset_samples(seagrassrare, samp_time == day)
  d
  
  # Ordinate your data to make the NMDS plot
  seagrass.ord <- ordinate(d, "NMDS", "bray")
  
  # Check the stress for the solution
  print(seagrass.ord)
  
  # Save the ordination in the list with a unique name
  stress.day[[day]] <- seagrass.ord$stress
  
  # Make the NMDS plot
  plot <- plot_ordination(d, seagrass.ord, color="treat") +
    geom_point(size=3) + theme_classic() +
    ggtitle(paste("Day", day)) +theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_line(linewidth = 1))
  
  plot_data <- as.data.frame(scores(seagrass.ord, display = "sites")) # Extract site scores
  colnames(plot_data) <- c("Axis1", "Axis2")
  
  #Add stress for that day to the plot
  plot<- plot +
    annotate("text", x = max(plot_data$Axis1)-0.4, y = max(plot_data$Axis2),
             label = paste("Stress =", round(stress.day[[day]], 4)), 
             hjust = 0, size = 4, color = "black")
  print(plot)
  
  # Save the plot in the list with a unique name
  plots[[day]] <- plot
  
  # PERMANOVA analysis
  # Calculate the distance within your data using the same ordination method as you do for your PERMANOVA
  distanceinfo <- phyloseq::distance(d, method = "bray")
  
  # Get your metadata out from phyloseq
  sample_df <- data.frame(sample_data(d))  # Use 'd' instead of 'seagrassrare'
  
  # Check lengths and unique values
  print(length(distanceinfo))
  print(length(sample_df$treat))
  print(unique(sample_df$treat))
  
  # Calculate the betadispersion within each region
  treat.dispers <- betadisper(distanceinfo, sample_df$treat)
  
  # Betadispersion test to see if all regions have the same beta dispersion
  beta.treat <- permutest(treat.dispers)
  
  # Save the permanova outputs in the list with a unique name
  perm.assum[[day]] <- beta.treat
  
  print(perm.assum)
  
  # PERMANOVA Analysis
  perm <- adonis2(distanceinfo ~ treat, data = sample_df, method = "bray")
  
  # Save PERMANOVA results
  permanova[[day]] <- perm
  
  # Extract the p-value
  p_value[[day]] <- perm$`Pr(>F)`[1]
}

#arrange the plots into one 
grid.arrange(plots[["d0"]],plots[["d1"]],plots[["d3"]],plots[["d7"]], 
             plots[["d14"]], plots[["d21"]], plots[["d28"]])


#Call up each days ordination
stress.day

#Call assumption for permanova
perm.assum

#Call up each days Permanova
permanova

# Print the p-value
print(p_value)

# Apply the Benjamini-Hochberg correction
adjusted_p_values <- p.adjust(p_value, method = "BH")

