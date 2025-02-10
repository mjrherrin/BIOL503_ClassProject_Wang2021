###BIOL503_Lab5_TaxaPlot_2025-02-10###

#Set working directory
setwd("C:/Users/allyw/Desktop/UBC_Master_Zoology/Courses_2025/Winter_Term_2025/BIOL503_MicrobialEcology/BIOL503_Lab/BIOL503_Lab5")
getwd ()

#read in data
seagrassdataset<-readRDS("filtered_seagrass_wMeta_phyloseq.RDS")

## load libraries
library(phyloseq)
library(tidyverse)
library(plyr)
library(qualpalr)
library(ggpubr)


## only plot root communities today
seagrass = subset_samples(seagrassdataset, samp_type=="RT")

# summarize at rank (or change this to be the taxonomic level you have in your dataset)
seagrass = tax_glom(seagrass, taxrank = "genus")

## calculate the number of reads in each sample. This is important for relative abundance calculations later
seagrass@sam_data$read_depth = sample_sums(seagrass)

## get seagrass data out of phyloseq and into a dataframe
# metadata
meta = as.data.frame(as.matrix(seagrass@sam_data)) |>
  rownames_to_column(var="sample_id")
metacols = ncol(meta) +1 # add 1 because of the full_join later to then pivot

# asv table
otu = as.data.frame(t(as.matrix(seagrass@otu_table))) |>
  rownames_to_column(var="sample_id")
# taxonomy
tax = as.data.frame(as.matrix(seagrass@tax_table)) |>
  rownames_to_column(var="placeholder_name")
## combine metadata and otu data
mo = full_join(meta, otu)

## pivot mo longer to be able to join with taxonomy
mo.long = mo |> pivot_longer(cols=-c(1:metacols),
                             names_to = "placeholder_name",
                             values_to = "taxa_abundance")

## combine with taxonomy
seagrassdf = full_join(mo.long, tax)

## caluclate relative abundance of each genus within each sample
seagrassdf$relativeabundance = as.numeric(seagrassdf$taxa_abundance)/as.numeric(seagrassdf$read_depth)

seagrassdf$plotnames = paste0(seagrassdf$order, ";", seagrassdf$genus)

## summarize data by what you want to plot
seagrass.sum = ddply(seagrassdf, c("plotnames"),
                     summarise,
                     sum = sum(relativeabundance))
## sort data by relative abundance. This is how you pick the most abundant taxa
sorted = seagrass.sum[order(-seagrass.sum$sum),]
## get top 15 genera
top = sorted[c(1:15),]
top$place = "top"

## join the top taxa and existing dataframe
## left_join removes the taxa that are not in top
alldata_tops = left_join(top, seagrassdf)

## caluclate relative abundance of other taxa
# get the per sample relative abundance of all the taxa that are left (in the top taxa)
# here I'm including region to keep it in the metadata of this all other dataset
# this is for subsetting and makign the plot later.
allothers = ddply(alldata_tops, c("sample_id", "treat"),
                  summarise,
                  top_taxa_sumra = sum(relativeabundance))
# get what the others value is supposed to be
# per-sample relative abundances sum to 1
# 1 - the total relative abundance of the top taxa = relative abundance of other taxa
allothers$relativeabundance = 1 - allothers$top_taxa_sumra
# combine datasets to have the relative abundan of top taxa and then the rest of the sample be others
alldata = full_join(alldata_tops, allothers)

## make the empty "place" cells say bottom. This workes because we used full_join
alldata$place = replace(alldata$place, is.na(alldata$place), "bottom")
## replace plot_names that have bottom taxa as their "place" with Other
alldata[alldata$place == "bottom",]$plotnames <- "Others"

# 1. find out how many colors you need
numcol <- length(unique(alldata$plotnames))
# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(15)
# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")
# 4. Extract hex colors
hex = as.data.frame(newpal$hex)
colnames(hex) <- c("taxa_color")
# 5. Get list of taxa
tops = as.data.frame(c(unique(alldata$plotnames)))
colnames(tops) <- c("plotnames")

#6. Join color list and taxa names
topcolors = cbind(tops, hex)
# 7. for the "others" plot name, replace that with grey 90 (this is just an astetic thing)
topcolors[topcolors$plotnames == "Others",]$taxa_color <- "grey90"
# 8. Make an object R can pull form for the colors
plotcolors <- topcolors$taxa_color
names(plotcolors) <- topcolors$plotnames

## order by decreasing relative abundance
alldata = alldata[order(-alldata$relativeabundance),]
## get list of factors in order
natural.genus.order = as.list(c(unique(alldata$plotnames)))
## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']
## add Others to end of list
plot.order = append(no.others, "Others")
## set plot_names levels
plot.order = unlist(plot.order)
## order dataframe by relative abundance
alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))

## these plots get pretty big, so let's only plot the goose region
wash = subset(alldata, alldata$treat=="Wash")

## make the plot
ggplot(wash, aes(x=sample_id, y=as.numeric(relativeabundance),
                  fill=as.factor(plotnames)))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=plotcolors)+
  guides(fill=guide_legend(ncol=1))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size=10, face="bold"),
        strip.text = element_text(color="black", size=10),
        legend.text=element_text(size=6),
        axis.line = element_line(colour = "black"),
        
        axis.text.x = element_blank())+
  labs(y="Relative Abundance", x="Sample", fill="Taxa")
