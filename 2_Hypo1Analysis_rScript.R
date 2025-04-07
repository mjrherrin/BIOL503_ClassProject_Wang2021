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
seagrassrare<-readRDS("OUTPUT_rarefied_seagrass_default_phyloseq_2025_03_10.RDS")

###Hypothesis: Recovery of the rhizosphere post-transplantation will be different between sod & wash treatments due to different magnitudes of disturbances in microbial community composition during transplanting.

#Set.seed so that NDMS plots are reproducible
set.seed(4)

#create data frame speficying days for the for loop
days <- c("d0", "d1", "d3", "d7", "d14", "d21", "d28")

#create a list where plots, the stress values, and the permanova outputs from the for-loop can be saved
  #Following loop was used to generate Figure 2 please disregard the stats as a nested PERMANOVA (see below) is better for what we are trying to analyse
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
    geom_point(size=3) +
    theme_classic() +
    ggtitle(paste("Day", day)) +theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_line(linewidth = 1))+
    scale_color_manual(values =c("Wash" = "blue", "Sod" = "red"))+
    theme(legend.position = "none")
 
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
  
  # Calculate the betadispersion within each treatments
  treat.dispers <- betadisper(distanceinfo, sample_df$treat)
  
  # Betadispersion test to see if all rtreatments have the same beta dispersion
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


#Call up each days stress value
stress.day

#Call assumption for permanova
perm.assum

#Call up each days Permanova
permanova

# Print the p-value
print(p_value)

# Apply the Benjamini-Hochberg correction
adjusted_p_values <- p.adjust(p_value, method = "BH")


###Nested PERMANOVA --> the following stats were used for the poster
# Calculate the Bray-Curtis distances
overall.distanceinfo <- phyloseq::distance(seagrassrare, method = "bray")

# Get metadata from phyloseq
overall.sample_df <- data.frame(sample_data(seagrassrare))  # Correctly retrieve sample data

# Ensure that 'day' and 'treat' are factors
overall.sample_df$samp_time <- as.factor(overall.sample_df$samp_time)
overall.sample_df$treat <- as.factor(overall.sample_df$treat)

# Betadispersion Analysis (optional, testing assumptions)
overall.dispers <- betadisper(overall.distanceinfo, overall.sample_df$treat)

# Test for homogeneity of dispersions (assumption for PERMANOVA)
beta.overall <- permutest(overall.dispers)

# get data frames
otu = as.data.frame(t(as.matrix(seagrassrare@otu_table)))
meta = as.data.frame(as.matrix(seagrassrare@sam_data))

## merge the data
mo = merge(meta, otu, by=0)

## count colum numbers
metacols = ncol(meta)+1

## run permanova to test if seagrasses from different regions are different
perm= adonis2(mo[,-c(1:metacols)] ~ samp_time/treat,
              data=mo, method = "bray")
perm
    #For Figure 1

# Run pairwise PERMANOVA with nested factors (post-hoc test)
pairwise_results <- pairwise.adonis(
  mo[,-c(1:metacols)],                   # OTU table
  factors = interaction(mo$samp_time, mo$treat),  # Nested factors as interaction terms
  sim.method = "bray",                   # Bray-Curtis dissimilarity
  p.adjust.m = "BH"                      # Benjamini-Hochberg correction for p-values
)
  #For Figure 1 and 2

# Print pairwise PERMANOVA results
print(pairwise_results)

#save as csv
write.csv(pairwise_results, file ="BIOL503_pairwisePERMANOVA_results.csv", row.names = TRUE)


####Overall NMDS (Figure 1)

# Ordinate your data to make the NMDS plot
overall.seagrass.ord <- ordinate(seagrassrare, "NMDS", "bray")

# Save the ordination in the list with a unique name
stress <- overall.seagrass.ord$stress


# Define the gradient colors for each treatment
treatment1_colors <- c("#FFCCCC", "#FF9999", "#FF6666", "#FF3333", "#FF0000", "#CC0000", "#990000")
treatment2_colors <- c("#CCCCFF", "#9999FF", "#6666FF", "#3333FF", "#0000FF", "#0000CC", "#000099")


# Combine gradient colors for both treatments
time_colors <- c("d0_Sod" = treatment1_colors[1], "d1_Sod" = treatment1_colors[2], 
                 "d3_Sod" = treatment1_colors[3], "d7_Sod" = treatment1_colors[4], 
                 "d14_Sod" = treatment1_colors[5], "d21_Sod" = treatment1_colors[6], 
                 "d28_Sod" = treatment1_colors[7], "d0_Wash" = treatment2_colors[1], 
                 "d1_Wash" = treatment2_colors[2], "d3_Wash" = treatment2_colors[3], 
                 "d7_Wash" = treatment2_colors[4], "d14_Wash" = treatment2_colors[5], 
                 "d21_Wash" = treatment2_colors[6], "d28_Wash" = treatment2_colors[7])

# Extract sample data from the phyloseq object
sample_data <- data.frame(sample_data(seagrassrare))

# Create a combined variable for time and treatment
sample_data$samp_time_treat <- paste(sample_data$samp_time, sample_data$treat, sep="_")

# Add the combined variable back to the phyloseq object
sample_data(seagrassrare) <- sample_data

# Make the NMDS plot
overall.plot <- plot_ordination(seagrassrare, overall.seagrass.ord, color="samp_time_treat", shape="treat") +
  geom_point(size=3) + 
  scale_color_manual(values=time_colors) +  # Use custom colors for time points within treatments
  scale_shape_manual(values=c(16, 17)) +  # Different shapes for treatments (e.g., 16 for circles, 17 for triangles)
  theme_classic() +
  # not sure if we want the circle in there its a bit much stat_ellipse(aes(group=samp_time_treat), type="norm", linetype=2) +  # Add ellipses for each treatment group
  theme(panel.border = element_rect(color="black", fill=NA, linewidth=1),
        axis.line = element_line(linewidth=1))

print(overall.plot)

overall.plot_data <- as.data.frame(scores(overall.seagrass.ord, display = "sites")) # Extract site scores
colnames(overall.plot_data) <- c("Axis1", "Axis2")

#Add stress for that day to the plot
overall.plot<- overall.plot +
  annotate("text", x = max(overall.plot_data$Axis1)-0.4, y = max(overall.plot_data$Axis2),
           label = paste("Stress = 0.201"),
           hjust = 0, size = 4, color = "black")
print(overall.plot)
