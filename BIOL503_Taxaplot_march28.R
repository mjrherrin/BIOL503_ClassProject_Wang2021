###BIOL503_TaxaPlot_2025-03-26###

#Set working directory
#setwd("C:/Users/allyw/Desktop/UBC_Master_Zoology/Courses_2025/Winter_Term_2025/BIOL503_MicrobialEcology/BIOL503_Lab/BIOL503_Lab5")
getwd ()

#read in data
seagrassdataset<-readRDS("~/UBC/biol503/filtered_seagrass_wMeta_phyloseq.RDS")

## load libraries
library(tidyverse)
library(phyloseq)
library(plyr)
library(qualpalr)
library(ggh4x)
library(vegan)
#library(ggpattern) #has not been loading properly but doesnt seeem to be needed
library(ggpubr)
library(dplyr)
library(ggplot2)+theme_set(theme_bw()+
                             theme(strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 12),
                                   axis.text.x = element_text(colour = "black", size = 8),
                                   legend.text = element_text(size = 8, colour ="black"),
                                   legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                                   legend.title = element_text(size = 14, colour = "black", face = "bold"),
                                   legend.key=element_blank(),
                                   # axis.ticks = element_blank(),
                                   panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank()
                             ))

## load functions
dephyloseq = function(phylo_obj){
  
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  
  ## how many metadata columns you have 
  metacols = ncol(meta)+1
  
  ## get out the otu table 
  ## if your metadta is empty after running this, you need to use 
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  #otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  
  ## get out the taxonomy file 
  tax = as.data.frame(phylo_obj@tax_table)
  
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="ASVid")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "ASVid", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}

## For plotting only sod communities today
sod = subset_samples(seagrassdataset, treat=="Sod")
## For plotting only Wash communities today
wash = subset_samples(seagrassdataset, treat=="Wash")

seagrass = merge_phyloseq(sod, wash)
seagrass

seagrass = prune_taxa(taxa_sums(seagrass)>0, seagrass)
seagrass

## calculate the number of reads in each sample. This is important for relative abundance calculations later
seagrass@sam_data$rd_depth = sample_sums(seagrass)
mean(seagrass@sam_data$rd_depth)
sum(seagrass@sam_data$rd_depth)

samplist = rownames(seagrass@sam_data)

#### SET UP GROUPS ####

## group at genus level
seagrass.gen = tax_glom(seagrass, taxrank = "genus")

## calculate read depth
seagrass.gen@sam_data$rd_depth = sample_sums(seagrass.gen)

## only sod
only.sod = subset_samples(seagrass.gen, treat == "Sod")

## only Wash
only.wash = subset_samples(seagrass.gen, treat == "Wash")

## get data out of phyloseq
seagrassdf = dephyloseq(seagrass.gen)

## make groups
grouplist = c(unique(seagrassdf$samp_time))

## calculate relative abundance
seagrassdf$ra = as.numeric(seagrassdf$asv_abundance)/as.numeric(seagrassdf$rd_depth)

## make plotnames
seagrassdf$plotnames = paste0(seagrassdf$order,"; ", seagrassdf$genus)

## summarize data by taxaplot group type. 
seagrassdf.sum = ddply(seagrassdf, c("samp_time", "plotnames", "samp_type", "treat"),
                    summarise,
                    sum = sum(ra))

## sort data by relative abundance. This is how the loop will pick the mos tabundant taxa
sorted = seagrassdf.sum[order(-seagrassdf.sum$sum),]

#### PICK TOP TAXA FOR EACH GROUP ####
## make empty dataframe to store output from the loop
top.df = NULL

## start loop
for(i in grouplist) {
  for(j in i) {
    
    ## subset dataframe by samples
    #!# Remeber to change te substrate to your group! 
    sample = subset(sorted, sorted$samp_time %in% c(j))
    
    ## get top 15 genera
    top = sample[c(1:10),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    
    ## close loop 
  }
}

##### SUBSET AND JOIN DATAFRAMES #####
toplist = c(unique(top.df$plotnames))

## calculate read depth
seagrass.gen@sam_data$rd_depth = sample_sums(seagrass.gen)

## get out of phyloseq
seadf = dephyloseq(seagrass.gen)
## make plotnames column
seadf$plotnames = paste0(seadf$order,"; ", seadf$genus)
## only keep target taxa
seagrassdf.sub = subset(seadf, seadf$plotnames %in% c(toplist))


## add info about groups from taxaplot
seasubtop = full_join(seagrassdf, top.df)
seasubtop = seasubtop[,c("Row.names", "plotnames", "samp_time", "samp_type", "treat")]

seadf.taxaplot = right_join(seasubtop, seagrassdf.sub)

## calculate relative abundance
seadf.taxaplot$ra = as.numeric(seadf.taxaplot$asv_abundance)/as.numeric(seadf.taxaplot$rd_depth)
head(seadf.taxaplot$ra)

## group the dataframe to get mean RA for each taxa
tp = ddply(seadf.taxaplot, 
           c("plotnames", "samp_time", "Row.names", "treat", "samp_type", "ra"),
           summarise,
           top_taxa_sumra = sum(ra),
           meanra = mean(ra),
           sdra = sd(ra))

## calculate others
tpothers = ddply(tp,  c("samp_time", "Row.names","treat", "samp_type", "ra"),
                 summarise,
                 sumra = sum(meanra))
tpothers$meanra = 1-tpothers$sumra
tpothers$plotnames = "Others"

alldata = full_join(tp, tpothers)

####GET COLORS FOR TAXAPLOT #####

# 1. find out how many colors you need
numcol <- length(unique(alldata$plotnames))

# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(3)

# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")

# 4. Extract hex colors
hex = as.data.frame(newpal$hex)
colnames(hex) <- c("taxa_color")

# 5. Get list of taxa
tops = as.data.frame(c(unique(alldata$plotnames)))
colnames(tops) <- c("plotnames")

# 6. Join color list and taxa names
topcolors = cbind(tops, hex)

# 7. for the "others" plot name, replace that with grey 90 (this is just an astetic thing)
topcolors[topcolors$plotnames == "Others",]$taxa_color <- "grey90"

# 8. Make an object R can pull form for the colors
plotcolors <- topcolors$taxa_color
names(plotcolors) <- topcolors$plotnames

## ORDER THE TAXA SO OTHERS ARE AT THE BOTTOM #####
## order by decreasing relative abundance
alldata = alldata[order(-alldata$meanra),]

## get list of factors in order
natural.genus.order = as.list(c(unique(alldata$plotnames)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")

## order disease
alldata$samp_time = factor(alldata$samp_time, levels=c("d0", "d1", "d3", "d7", "d14", "d21", "d28"))

## set plot_names levels
plot.order = unlist(plot.order)

##### MAKE PLOT #####
alldata = subset(alldata, alldata$samp_time!="NA")
soddata = subset(alldata, alldata$treat=="Sod")
washdata = subset(alldata, alldata$treat=="Wash")

## order dataframe by relative abundance
alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))

#alldata$samp_time = ifelse(alldata$sample_target %in% c("seawater", "airline"), "Symp.Env.", alldata$pink_health)

bp=ggplot(alldata, aes(x=samp_time, treat, samp_time, as.character(Row.names),
                       y= as.numeric(mean(ra)),
                       fill=plotnames))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=plotcolors)+
  theme(axis.text.x = element_text(angle=90))+
  facet_nested(.~samp_time+samp_type+treat, scales="free", space="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=guide_legend(ncol=1))
bp

sodplot=ggplot(soddata, aes(x=paste(samp_time, treat, samp_time, as.character(Row.names)),
                       y=as.numeric(ra),
                       fill=plotnames))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=plotcolors)+
  theme(axis.text.x = element_text(angle=90))+
  facet_nested(.~samp_time+samp_type, scales="free", space="free")+
  theme(axis.text.x = element_blank())+
  guides(fill=guide_legend(ncol=1))
sodplot

# WHY DOES IT GO ABOVE 1!

washplot=ggplot(washdata, aes(x=paste(samp_time, treat, samp_time, as.character(Row.names)),
                            y=as.numeric(ra),
                            fill=plotnames))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=plotcolors)+
  theme(axis.text.x = element_text(angle=90))+
  facet_nested(.~samp_time+samp_type, scales="free", space="free")+
  theme(axis.text.x = element_blank())+
  guides(fill=guide_legend(ncol=1))
washplot

ggarrange(sodplot,washplot, ncol=1, heights=c(0.4, 1))

#######
only.algicola.df = dephyloseq(only.algicola)
only.algicola.df$ra = as.numeric(only.algicola.df$asv_abundance)/as.numeric(only.algicola.df$rd_no_off_with_algicola)
only.algicola.df$presabs = ifelse(only.algicola.df$ra==0, "abs", "pres")

only.algicola.df$sample_target = factor(only.algicola.df$sample_target, levels=c("spool", "airline", "seawater"))
only.algicola.df$pink_health = factor(only.algicola.df$pink_health, levels=c("pink_healthy", "pink_diseased","airline", "seawater"))

ggplot(only.algicola.df, aes(x=paste(location, date_sampled, as.character(Row.names)),
                             y=paste(genus, "   spacing"),
                             size=sqrt(ra), alpha=presabs, color=genus))+
  geom_point()+
  facet_nested(.~year_sampled+pink_health, scales="free", space="free")+
  scale_color_manual(values=c("magenta"))+
  scale_size(range=c(2,8),
             breaks=c(0,0.2, 0.4, 0.6, 0.8),labels=c("0","0.2","0.4","0.6","0.8"),guide="legend")

only.algicola.df$ra = as.numeric(ifelse(only.algicola.df$ra>0, only.algicola.df$ra, "NA"))

dp=ggplot(only.algicola.df, aes(x=paste(location, date_sampled, as.character(Row.names)),
                                y=paste(genus, "   spacing"),
                                fill=sqrt(ra)))+
  geom_tile()+
  facet_nested(.~year_sampled+pink_health, scales="free", space="free")+
  scale_fill_gradient(low="#b4b4b4", high="black", na.value = "white")+
  theme(axis.text.x = element_blank())

ggarrange(dp,bp, ncol=1, heights=c(0.4, 1))

ggsave("C:/Users/siobh/OneDrive - The University Of British Columbia/Project - Outplant Saccharina/git_cascadia_outplant/Project - kelp_disease/output/taxaplot_NoAlgicola.pdf",
       width=11.5, height=8.5, units="in")



########################################################################################
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
ggplot(wash, aes(x=as.character(Row.names), y=as.numeric(ra),
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

