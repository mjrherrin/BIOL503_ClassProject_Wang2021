###BIOL503_Project_FilterData_2025-02-10###

## load libraries
library(tidyverse)
library(phyloseq)
library(vegan)

##Set working directory --> Command setwd() tells R where to find the data file
setwd("C:/Users/allyw/Desktop/UBC_Master_Zoology/Courses_2025/Winter_Term_2025/BIOL503_MicrobialEcology/BIOL503_Lab/BIOL503_Lab5")
getwd ()

# RDS files can only be opened in R --> readRDS is the function to do this
nonrare = readRDS("~/miranda/Biol503/raw_data/Wang_seagrass_unfiltered_phyloseq.RDS")


###INVESTIGATE PHYOSEQ OBJECT
## View the 3 different parts of the phyloseq object (use @ to view specific tables)
# metadata
View(nonrare@sam_data)
# otu table
View(nonrare@otu_table)
# view taxonomy
View(nonrare@tax_table)


###FILTERING 
# remove off target taxa
nonrare = subset_taxa(nonrare,
                      domain!="Unassigned"&
                        species!="Chloroplast" &
                        species!="Mitochondria" &
                        domain!="Eukaryota")

## look at your taxonomy table to check everyhting got removed with View()
#after first round of filtering there were still mitochondrial categories in rank 5 so added additional restriction to subset function

## print the sample sums (number of reads) of each sample in the console
sample_sums(nonrare)

## make a plot of the sample sums
plot(sort(sample_sums(nonrare)))

## add total number of reads in each sample to the metadata
nonrare@sam_data$sample_sums_unfiltered = as.numeric(sample_sums(nonrare))

##Investigate min and max reads per sample
min(nonrare@sam_data$sample_sums_unfiltered)
#Output:5
max(nonrare@sam_data$sample_sums_unfiltered)
#Output:11932

## remove the samples by your threshold (equal or greater to 6000, because this is where most sample lie)
nonrare.high <- prune_samples(sample_sums(nonrare) >= 6000, nonrare)

## make another object with only the samples you lost (less than 6000 reads)
nonrare.below6000 <- prune_samples(sample_sums(nonrare) < 6000, nonrare)

## get the metadata out of phyloseq from your low reads object --> 15 samples below cutoff
nonrare.below6000 = as.matrix(nonrare.below6000@sam_data)

## write file to know which samples were lost here. This is important for the methods section.
write.csv(nonrare.below6000, "default_seagrass_samples_less_than_6000.csv")


## extract otu dataframe (asv table) from phyloseq object
# notice the t before as matrix here! This is pivoting the otu table 90 degrees.
otutab <- as.data.frame(t(as.matrix(otu_table(nonrare.high@otu_table))))

## check your asvtab dataframe visually to see if your ASVs are the ROWS an the sample names are the COLUMNS.
view(otutab)

# This is critical otherwise you will not filter to correct thing.
## calculate the sum of each row in your otu table
otutab$asv_abundance = rowSums(otutab)

## get the minimum asv_abundance value
min(otutab$asv_abundance)
#Output: 0, low abundance ASVs indicates polymerase mistake --> remove low ab.

## remove low frequency asvs, we chose 100 because bellow that seemed low but above reasonable
otu.pruned = subset(otutab, otutab$asv_abundance>=100)

## now get the minimum asv abundance from your cleaned up otu table.
# this should be at or above your filtering treshold.
min(otu.pruned$asv_abundance)
#Output: 100

## remove asv_abundance column from your OTU table since we don't want to analyse it
# our asv_abundance column gets tacked onto the end when you make it, so you just need to delete the last ## how many columns in dataframe?
widthotu = ncol(otu.pruned)
#Output: 106 columns in otu.prume

## see how widthotu appears as a Value in the environement?
## keep everything in the otu.pruned dataset except the last columns
otu.pruned = otu.pruned[,-c(widthotu)] 
#run line 89 again to confirm that a column was removed
#Output after removing last column: 105 columns in otu.prume


## create a function to count the number of occurence along rows where the number in the cell (ASV occurence)
ASVoccur = function(x){return(sum(x>0))}

## calculate the occurence of each ASV in your dataframe
otu.pruned$asv_occur_count = apply(otu.pruned, 1, ASVoccur)

## see what the ASV occurence counts look like
summary(otu.pruned$asv_occur_count)
# #Output:
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    8.00   14.00   22.53   28.00  105.00  

## lets remove ASVs found two or less times in the dataset
otu.highfreq = subset(otu.pruned, otu.pruned$asv_occur_count>2)

## see if the filtering worked.The minimum should now be above your threshold (2 in this case)
summary(otu.highfreq$asv_occur_count)
# #Output: 
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.00    8.00   15.00   23.28   29.00  105.00 

## remove the asv_occur_count column
otu.highfreq = otu.highfreq[,-c(widthotu)]


##De-noising the data
otu.clean <- mutate_all(otu.highfreq, funs(ifelse(. < 5, 0, .)))


## make the cleaned phyloseq object
cleanseagrass = phyloseq(sample_data(nonrare.high@sam_data),
                         tax_table(nonrare.high@tax_table),
                         otu_table(as.matrix(otu.clean), taxa_are_rows = TRUE))

## basic check to see that the object was created
cleanseagrass

## add the number of reads after filtering
cleanseagrass@sam_data$sample_sums_filtered = sample_sums(cleanseagrass)

## save rarefied dataframe
write_rds(cleanseagrass, "filtered_seagrass_default_phyloseq.RDS")

library(tidyr)
#this is seperating and cleaning our metadata, which had all this imforamtion together
meta <- separate(as.data.frame(as.matrix(cleanseagrass@sam_data)), col=Library.Name, 
                 into=c("ID", "treat", "samp_time", "samp_type","remove"), sep="_")

ps = phyloseq(sample_data(meta),
              tax_table(cleanseagrass@tax_table),
              otu_table(cleanseagrass@otu_table, taxa_are_rows = F))

###RAREFACTION
## use the rarecurve function in the package vegan to plot the rarefaction
## there are a lot of things going on here because you need your SampleID to be the row name and the ASV ## If you open the otu.clean dataframe, you will see that it's not in that format, so we need to fix this.
## R works like math, you read the inner most () first and work outwards
rarecurve(
  as.data.frame( ## 3. Turn the matrix back into a dataframe
    t( ## 2. turn the otu.clean matrix 90 degrees to have the correct orientation
      as.matrix( ## 1. turn the otu.clean dataframe into a matrix
        otu.clean))), ## 0. the dataframe we are doing stuff to
  step=50, cex=0.5, label=FALSE) ## now that the data are formatted, sample the reads (these are just parameters


## open the cleanseagrass metadata
View(cleanseagrass@sam_data)

## set seed to tell R which set of random numbers to use.
# This is important because it will allow you to sample randomly the same way every time
set.seed(5)
## rarefy every sample to a set number of reads here we chose 5686 because the rarification curves level off just before that meaning we capture most of the species at that sampple size
raresg <- rarefy_even_depth(cleanseagrass, sample.size = 5686)


## calculate the rarefied sample sums
raresg@sam_data$rare_sample_sums = sample_sums(raresg)
## see distribution of rare_sample_sums
summary(raresg@sam_data$rare_sample_sums)
#Output: 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2341    2341    2341    2341    2341    2341 

## save rarefied dataframe
write_rds(raresg, "rarefied_seagrass_default_phyloseq.RDS")



