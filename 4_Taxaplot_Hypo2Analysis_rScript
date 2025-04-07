###BIOL503_TaxaPlot_2025-03-26###

#Set working directory
#setwd("C:/Users/Spencer/OneDrive/Documents/miranda/Biol503")
getwd ()

## load libraries
library(tidyverse)
library(phyloseq)
library(plyr)
library(qualpalr)
library(ggh4x)
library(vegan)
library(ggpubr)
library(dplyr)
library(ggplot2)+theme_set(theme_bw()+
                             theme(strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 18),
                                   axis.text.x = element_text(colour = "black", size = 18),
                                   legend.text = element_text(size = 14, colour ="black"),
                                   legend.position = "right", 
                                   axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 12, colour = "black"),
                                   legend.title = element_text(size = 14, colour = "black", face = "bold"),
                                   plot.title = element_text(size = 18, face = "bold"),
                                   # Make the legend more compact vertically
                                   legend.key.height = unit(0.5, "line"),
                                   legend.spacing.y = unit(0.1, "cm"),
                                   legend.margin = margin(0, 0, 0, 0),
                                   legend.box.margin = margin(0, 0, 0, 0),
                                   legend.key=element_blank(),
                                   # axis.ticks = element_blank(),
                                   panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank()
                             ))
#read in data
seagrassdataset<-readRDS("~/UBC/biol503/filtered_seagrass_wMeta_phyloseq.RDS")

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

## Sod subset
sod = subset_samples(seagrassdataset, treat=="Sod")
## Wash subset
wash = subset_samples(seagrassdataset, treat=="Wash")

seagrass = merge_phyloseq(sod, wash)

#### SET UP GROUPS ####

## group at genus level
seagrass = tax_glom(seagrassdataset, taxrank = "genus")
sod = tax_glom(sod, taxrank = "genus")
wash = tax_glom(wash, taxrank = "genus")

## calculate read depth
seagrass@sam_data$read_depth = sample_sums(seagrass)
sod@sam_data$read_depth = sample_sums(sod)
wash@sam_data$read_depth = sample_sums(wash)

## get data out of phyloseq
seagrassdf = dephyloseq(seagrass)
soddf = dephyloseq(sod)
washdf = dephyloseq(wash)

## make groups
grouplist = c(unique(seagrassdf$samp_time))
#warning occurs this is OK
grouplistsod = c(unique(soddf$samp_time))
grouplistwash = c(unique(washdf$samp_time))

## calculate relative abundance
seagrassdf$ra = as.numeric(seagrassdf$asv_abundance)/as.numeric(seagrassdf$read_depth)
soddf$ra = as.numeric(soddf$asv_abundance)/as.numeric(soddf$read_depth)
washdf$ra = as.numeric(washdf$asv_abundance)/as.numeric(washdf$read_depth)

## make plotnames
seagrassdf$plotnames = paste0(seagrassdf$order,"; ", seagrassdf$genus)
soddf$plotnames = paste0(soddf$order,"; ", soddf$genus)
washdf$plotnames = paste0(washdf$order,"; ", washdf$genus)

## summarize data by taxaplot group type. 
seagrassdf.sum = ddply(seagrassdf, c("samp_time", "plotnames", "treat"),
                       summarise,
                       sum = sum(ra))
soddf.sum = ddply(soddf, c("samp_time", "plotnames", "treat"),
                  summarise,
                  sum = sum(ra))
washdf.sum = ddply(washdf, c("samp_time", "plotnames", "treat"),
                   summarise,
                   sum = sum(ra))

## sort data by relative abundance. This is how the loop will pick the most abundant taxa
sorted = seagrassdf.sum[order(-seagrassdf.sum$sum),]
sortedsod = soddf.sum[order(-soddf.sum$sum),]
sortedwash= washdf.sum[order(-washdf.sum$sum),]

#### PICK TOP TAXA FOR EACH GROUP ####
## make empty dataframe to store output from the loop
top.df = NULL
top.sod.df = NULL
top.wash.df = NULL

## start loop for all
for(i in grouplist) {
  for(j in i) {
    
    ## subset dataframe by samples
    sample = subset(sorted, sorted$samp_time %in% c(j))
    
    ## get top 15 genera
    top = sample[c(1:10),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    
    ## close loop 
  }
}

## start loop for sod
for(i in grouplistsod) {
  for(j in i) {
    
    ## subset dataframe by samples
    samples = subset(sortedsod, sortedsod$samp_time %in% c(j))
    
    ## get top 15 genera
    tops = samples[c(1:10),]
    
    ## save list of top  abundance taxa
    t.tmps <- tops
    top.sod.df <- rbind.fill(top.sod.df, t.tmps)
    
    ## close loop 
  }
}

## start loop for wash
for(i in grouplistwash) {
  for(j in i) {
    
    ## subset dataframe by samples
    samplew = subset(sortedwash, sortedwash$samp_time %in% c(j))
    
    ## get top 15 genera
    topw = samplew[c(1:10),]
    
    ## save list of top  abundance taxa
    t.tmpw <- topw
    top.wash.df <- rbind.fill(top.wash.df, t.tmpw)
    
    ## close loop 
  }
}

##### SUBSET AND JOIN DATAFRAMES #####
toplist = c(unique(top.df$plotnames))
toplistsod = c(unique(top.sod.df$plotnames))
toplistwash = c(unique(top.wash.df$plotnames))

#check if different
unique(toplist)
unique(toplistsod)
unique(toplistwash)

identical(toplist, toplistsod) && identical(toplist, toplistwash)
#they are different so continue doing each seperatly

## calculate read depth
seagrass@sam_data$rd_depth = sample_sums(seagrass)
sod@sam_data$rd_depth = sample_sums(sod)
wash@sam_data$rd_depth = sample_sums(wash)

## get out of phyloseq
seagrassdf = dephyloseq(seagrass)
soddf = dephyloseq(sod)
washdf = dephyloseq(wash)

## make plotnames column
seagrassdf$plotnames = paste0(seagrassdf$order,"; ", seagrassdf$genus)
soddf$plotnames = paste0(soddf$order,"; ", soddf$genus)
washdf$plotnames = paste0(washdf$order,"; ", washdf$genus)

## only keep target taxa
seagrass.sub = subset(seagrassdf, seagrassdf$plotnames %in% c(toplist))
sod.sub = subset(soddf, soddf$plotnames %in% c(toplistsod))
wash.sub = subset(washdf, washdf$plotnames %in% c(toplistwash))

## add info about groups from taxaplot
seasubtop = full_join(seagrassdf, top.df)
seasubtop = seasubtop[,c("Row.names", "plotnames", "samp_time", "treat")]
seadf.taxaplot = right_join(seasubtop, seagrass.sub)

sodsubtop = full_join(soddf, top.sod.df)
sodsubtop = sodsubtop[,c("Row.names", "plotnames", "samp_time", "treat")]
soddf.taxaplot = right_join(sodsubtop, sod.sub)

washsubtop = full_join(washdf, top.wash.df)
washsubtop = washsubtop[,c("Row.names", "plotnames", "samp_time", "treat")]
washdf.taxaplot = right_join(washsubtop, wash.sub)

## calculate relative abundance
seadf.taxaplot$ra = as.numeric(seadf.taxaplot$asv_abundance)/as.numeric(seadf.taxaplot$rd_depth)
head(seadf.taxaplot$ra)

soddf.taxaplot$ra = as.numeric(soddf.taxaplot$asv_abundance)/as.numeric(soddf.taxaplot$rd_depth)
head(soddf.taxaplot$ra)

washdf.taxaplot$ra = as.numeric(washdf.taxaplot$asv_abundance)/as.numeric(washdf.taxaplot$rd_depth)
head(washdf.taxaplot$ra)

## group the dataframe to get mean RA for each taxa
tp = ddply(seadf.taxaplot, 
           c("plotnames", "samp_time", "Row.names", "treat"),
           summarise,
           meanra = mean(ra),
           sdra = sd(ra))

tpsod = ddply(soddf.taxaplot, 
              c("plotnames", "samp_time", "Row.names", "treat"),
              summarise,
              meanra = mean(ra),
              sdra = sd(ra))

tpwash = ddply(washdf.taxaplot, 
               c("plotnames", "samp_time", "Row.names", "treat"),
               summarise,
               meanra = mean(ra),
               sdra = sd(ra))

## calculate others
tpothers = ddply(tp,  c("samp_time", "Row.names","treat"),
                 summarise,
                 sumra = sum(meanra))
tpothers$meanra = 1-tpothers$sumra
tpothers$plotnames = "Others"

alldata = full_join(tp, tpothers)

tpotherssod = ddply(tpsod,  c("samp_time", "Row.names","treat"),
                    summarise,
                    sumra = sum(meanra))
tpotherssod$meanra = 1-tpotherssod$sumra
tpotherssod$plotnames = "Others"

soddata = full_join(tpsod, tpotherssod)


tpotherswash = ddply(tpwash,  c("samp_time", "Row.names","treat"),
                     summarise,
                     sumra = sum(meanra))
tpotherswash$meanra = 1-tpotherswash$sumra
tpotherswash$plotnames = "Others"

washdata = full_join(tpwash, tpotherswash)

####GET COLORS FOR TAXAPLOT for  #####
# 1. Combine all unique taxa from all datasets to create a master taxa list
all_unique_taxa <- unique(c(unique(alldata$plotnames), 
                            unique(soddata$plotnames), 
                            unique(washdata$plotnames)))
# 2. Find out how many colors you need (based on the combined list)
numcol <- length(all_unique_taxa)

# 3. Use a seed for reproducibility
set.seed(3)

# 4. Create one master color palette for all taxa
master_palette <- qualpal(n = numcol, colorspace = "pretty")

# 5. Create a master dataframe with taxa names and colors
master_taxa_colors <- data.frame(
  plotnames = all_unique_taxa,
  taxa_color = master_palette$hex,
  stringsAsFactors = FALSE
)

# 6. Set "Others" to grey90
master_taxa_colors$taxa_color[master_taxa_colors$plotnames == "Others"] <- "grey90"

# 7. Create the color vectors for each dataset by looking up colors from master palette
plotcolors <- master_taxa_colors$taxa_color
names(plotcolors) <- master_taxa_colors$plotnames

# Now use the master color mapping to create consistent colors for each dataset
# This ensures the same taxon gets the same color in all datasets
plotcolorssod <- plotcolors[match(unique(soddata$plotnames), names(plotcolors))]
names(plotcolorssod) <- unique(soddata$plotnames)

plotcolorswash <- plotcolors[match(unique(washdata$plotnames), names(plotcolors))]
names(plotcolorswash) <- unique(washdata$plotnames)

## ORDER THE TAXA SO OTHERS ARE AT THE BOTTOM #####
## order by decreasing relative abundance
alldata = alldata[order(-alldata$meanra),]
soddata = soddata[order(-soddata$meanra),]
washdata = washdata[order(-washdata$meanra),]

## get list of factors in order
natural.genus.order = as.list(c(unique(alldata$plotnames)))
natural.genus.order.sod = as.list(c(unique(soddata$plotnames)))
natural.genus.order.wash = as.list(c(unique(washdata$plotnames)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']
no.others.sod=natural.genus.order.sod[!natural.genus.order.sod == 'Others']
no.others.wash=natural.genus.order.wash[!natural.genus.order.wash == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")
plot.order.sod = append(no.others.sod, "Others")
plot.order.wash = append(no.others.wash, "Others")

## order time and change d to day
alldata$samp_time = factor(alldata$samp_time, levels=c("d0", "d1", "d3", "d7", "d14", "d21", "d28"),
                           labels = c("day 0", "day 1", "day 3", "day7", "day 14", "day 21", "day 28"))

soddata$samp_time = factor(soddata$samp_time, levels=c("d0", "d1", "d3", "d7", "d14", "d21", "d28"),
                           labels = c("day 0", "day 1", "day 3", "day7", "day 14", "day 21", "day 28"))

washdata$samp_time = factor(washdata$samp_time, levels=c("d0", "d1", "d3", "d7", "d14", "d21", "d28"),
                            labels = c("day 0", "day 1", "day 3", "day7", "day 14", "day 21", "day 28"))

## set plot_names levels
plot.order = unlist(plot.order)
plot.order.sod = unlist(plot.order.sod)
plot.order.wash = unlist(plot.order.wash)

#### Make Functional groups ###############
#look at unique names for taxa to look into functional groups
unique(alldata$plotnames)
unique(soddata$plotnames)
unique(washdata$plotnames)
#all data has 19, sod 20 and wash 21. All unique taxa were looked up in liturature (see our appendix) to assign functional group to each taxa

#Asign top taxa to functional groups
# create a mapping table using only unique plotnames
function_mapping <- data.frame(
  plotnames = unique(alldata$plotnames)) %>%
  mutate(functional_group = case_when(
    plotnames %in% c("Campylobacterales; Sulfurimonas", 
                     "Campylobacterales; Arcobacteraceae", 
                     "Campylobacterales; Sulfurovum", 
                     "Milano-WF1B-44; Milano-WF1B-44") ~ "Sulphur oxidizing",
    
    plotnames %in% c("Alteromonadales; Colwellia", 
                     "Lachnospirales; Lachnospiraceae", 
                     "Desulfobulbales; Desulfocapsaceae") ~ "Strict facultative anaerobes",
    
    plotnames %in% c("Chromatiales; Sedimenticolaceae", 
                     "Bacteroidales; Marinifilaceae",
                     "Steroidobacterales; Woeseia") ~ "Not strict anaerobes",
    
    plotnames %in% c("Desulfobulbales; Desulforhopalus", 
                     "Bacteroidales; Bacteroidetes_BD2-2") ~ "Sulfate-reducing",
    
    plotnames %in% c("Nitrosococcales; Methylophagaceae") ~ "Methylotrophs",
    
    plotnames %in% c("Alteromonadales; Psychromonas",
                     "Gammaproteobacteria_Incertae_Sedis; Unknown_Family") ~ "Unknown",
    
    plotnames %in% c("Vibrionales; Aliivibrio", 
                     "Vibrionales; Vibrio") ~ "Vibrio",
    
    plotnames %in% c("Others") ~ "Others"))

function_mapping.sod <- data.frame(
  plotnames = unique(soddata$plotnames)) %>%
  mutate(functional_group = case_when(
    plotnames %in% c("Campylobacterales; Sulfurimonas", 
                     "Campylobacterales; Arcobacteraceae", 
                     "Campylobacterales; Sulfurovum", 
                     "Milano-WF1B-44; Milano-WF1B-44",
                     "Thiotrichales; Thiotrichaceae") ~ "Sulphur oxidizing",
    
    plotnames %in% c( "Lachnospirales; Lachnospiraceae", 
                      "LCP-89; LCP-89") ~ "Strict anaerobes",
    
    plotnames %in% c("Alteromonadales; Psychromonas",
                     "Chromatiales; Sedimenticolaceae", 
                     "Bacteroidales; Marinifilaceae",
                     "Steroidobacterales; Woeseia",
                     "Alteromonadales; Colwellia",
                     "Bacteroidales; Draconibacterium") ~ "Not strict anaerobes",
    
    plotnames %in% c("Desulfobulbales; Desulforhopalus", 
                     "Bacteroidales; Bacteroidetes_BD2-2",
                     "Desulfobulbales; Desulfocapsaceae") ~ "Sulfate-reducing",
    
    plotnames %in% c("Nitrosococcales; Methylophagaceae") ~ "Methylotrophs",
    
    plotnames %in% c("Gammaproteobacteria_Incertae_Sedis; Unknown_Family") ~ "Unknown",
    
    plotnames %in% c("Vibrionales; Aliivibrio", 
                     "Vibrionales; Vibrio") ~ "Vibrio",
    
    plotnames %in% c("Others", 
                     "Cellvibrionales; Halioglobus",
                     "Desulfobacterales; Desulfosarcina") ~ "Others"))

function_mapping.wash <- data.frame(
  plotnames = unique(washdata$plotnames)) %>%
  mutate(functional_group = case_when(
    plotnames %in% c("Campylobacterales; Sulfurimonas", 
                     "Campylobacterales; Arcobacteraceae", 
                     "Campylobacterales; Sulfurovum") ~ "Sulphur oxidizing",
    
    plotnames %in% c("Lachnospirales; Lachnospiraceae") ~ "Strict anaerobes",
    
    plotnames %in% c("Chromatiales; Sedimenticolaceae", 
                     "Bacteroidales; Marinifilaceae",
                     "Steroidobacterales; Woeseia",
                     "Alteromonadales; Colwellia",
                     "Alteromonadales; Psychromonas") ~ "Not strict anaerobes",
    
    plotnames %in% c("Desulfobulbales; Desulforhopalus",
                     "Desulfobulbales; Desulfocapsaceae") ~ "Sulfate-reducing",
    
    plotnames %in% c("Nitrosococcales; Methylophagaceae") ~ "Methylotrophs",
    
    plotnames %in% c("Gammaproteobacteria_Incertae_Sedis; Unknown_Family",
                     "Peptostreptococcales-Tissierellales; Fusibacter") ~ "Unknown",
    
    plotnames %in% c("Vibrionales; Aliivibrio", 
                     "Vibrionales; Vibrio") ~ "Vibrio",
    
    plotnames %in% c("Alteromonadales; Algicola",
                     "Alteromonadales; Glaciecola") ~ "Strict aerobic",
    
    plotnames %in% c("Others",
                     "Campylobacterales; Halarcobacter",
                     "Bacteroidales; Bacteroidetes_BD2-2",
                     "Milano-WF1B-44; Milano-WF1B-44") ~ "Others"))

# Check that the mapping table is correct and has only 18 rows
# (one for each unique plotname)
print(dim(function_mapping))
print(function_mapping)

print(dim(function_mapping.sod))
print(function_mapping.sod)
#has 20 as mentioned before
print(dim(function_mapping.wash))
print(function_mapping.wash)
#has 21 as mentioned before

# Now add the functional group to the original data
alldata <- left_join(alldata, function_mapping, by = "plotnames")
soddata <- left_join(soddata, function_mapping.sod, by = "plotnames")
washdata <- left_join(washdata, function_mapping.wash, by = "plotnames")

#Now fix the too many top taxa problem by changing the bottom extra top taxa to others
unique(soddata$plotnames)
soddata <- soddata %>%
  mutate(plotnames = case_when(
    plotnames %in% c("Cellvibrionales; Halioglobus", "Desulfobacterales; Desulfosarcina", "Bacteroidales; Draconibacterium") ~ "Others",
    TRUE ~ plotnames))
unique(soddata$plotnames)
#17 plot names good

unique(washdata$plotnames)
washdata <- washdata %>%
  mutate(plotnames = case_when(
    plotnames %in% c("Milano-WF1B-44; Milano-WF1B-44", "Bacteroidales; Bacteroidetes_BD2-2", "Campylobacterales; Halarcobacter", "Alteromonadales; Glaciecola") ~ "Others",
    TRUE ~ plotnames))
unique(washdata$plotnames)
#17 plot names now good

##### MAKE PLOT #####

## order dataframe by relative abundance
alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))

#plot all sod and wash together
bp=ggplot(alldata, aes(x=paste(samp_time, treat, as.character(Row.names)),
                       y= as.numeric(meanra),
                       fill=plotnames))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=plotcolors)+
  theme(axis.text.x = element_text(angle=90))+
  facet_nested(.~samp_time+treat, scales="free", space="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=guide_legend(ncol=1))
bp

#sod taxa plot
#order plotnamse
soddata$plotnames = factor(soddata$plotnames, levels=c(plot.order.sod))
#plot
sodplot=ggplot(soddata, aes(x=paste(samp_time, treat, as.character(Row.names)),
                            y=as.numeric(meanra),
                            fill=plotnames))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=plotcolorssod)+
  theme(axis.text.x = element_text(angle=90))+
  facet_nested(.~samp_time, scales="free", space="free")+
  theme(axis.text.x = element_blank())+
  guides(fill=guide_legend(ncol=1))+
  ggtitle("Sod") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample",
       y = "Relative abundance",
       fill = "Taxa")
sodplot

#Wash taxa plot
#order plot names
washdata$plotnames = factor(washdata$plotnames, levels=c(plot.order.wash))
#plot
washplot=ggplot(washdata, aes(x=paste(samp_time, treat, as.character(Row.names)),
                              y=as.numeric(meanra),
                              fill=plotnames))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=plotcolors)+
  theme(axis.text.x = element_text(angle=90))+
  facet_nested(.~samp_time, scales="free", space="free")+
  theme(axis.text.x = element_blank())+
  guides(fill=guide_legend(ncol=1))+
  ggtitle("Wash") +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Sample",
       y = "Relative abundance",
       fill = "Taxa")
washplot

#plot sod and wash plots together
ggarrange(sodplot, washplot, ncol = 1, heights = c(1, 1))


#### Functional group plots
##Color and order aisgn for functional groups
# 1. Combine all unique functional groups from all datasets
all_unique_functional_groups <- unique(c(
  unique(alldata$functional_group),
  unique(soddata$functional_group),
  unique(washdata$functional_group)
))

# 2. Find out how many colors you need (based on the combined list)
num_fg_colors <- length(all_unique_functional_groups)

# 3. Generate one master color palette for all functional groups
set.seed(3)
master_fg_palette <- qualpal(n = num_fg_colors, colorspace = "pretty")

# 4. Create a master dataframe with functional group names and colors
master_fg_colors <- data.frame(
  functional_group = all_unique_functional_groups,
  functional_group_color = master_fg_palette$hex,
  stringsAsFactors = FALSE
)

# 5. Set "Other" to grey90
master_fg_colors$functional_group_color[master_fg_colors$functional_group == "Others"] <- "grey90"

# 6. Create a named vector for easy lookup
all_fg_colors <- master_fg_colors$functional_group_color
names(all_fg_colors) <- master_fg_colors$functional_group

# 7. Create dataset-specific color vectors based on the master palette
# For alldata
functional_group_colors <- data.frame(
  functional_group = unique(alldata$functional_group),
  functional_group_color = all_fg_colors[match(unique(alldata$functional_group), names(all_fg_colors))],
  stringsAsFactors = FALSE
)

# For soddata
functional_group_colors.s <- data.frame(
  functional_group = unique(soddata$functional_group),
  functional_group_color = all_fg_colors[match(unique(soddata$functional_group), names(all_fg_colors))],
  stringsAsFactors = FALSE
)

# For washdata
functional_group_colors.w <- data.frame(
  functional_group = unique(washdata$functional_group),
  functional_group_color = all_fg_colors[match(unique(washdata$functional_group), names(all_fg_colors))],
  stringsAsFactors = FALSE
)

# 8. Create plotcolors objects for ggplot
functioncolors <- functional_group_colors$functional_group_color
names(functioncolors) <- functional_group_colors$functional_group

functioncolors.s <- functional_group_colors.s$functional_group_color
names(functioncolors.s) <- functional_group_colors.s$functional_group

functioncolors.w <- functional_group_colors.w$functional_group_color
names(functioncolors.w) <- functional_group_colors.w$functional_group

# ORDER FUNCTIONAL GROUPS WITH OTHERS AT THE BOTTOM #####
# 9. Order functional groups (keeping the rest of your ordering code)
alldata <- alldata[order(-alldata$meanra),]
functional_group_order <- as.list(unique(alldata$functional_group))
no.others.F <- functional_group_order[!functional_group_order == 'Others']
functional_group_order <- append(no.others.F, "Others")
alldata$functional_group <- factor(alldata$functional_group, 
                                   levels = functional_group_order)

soddata <- soddata[order(-soddata$meanra),]
functional_group_order.s <- as.list(unique(soddata$functional_group))
no.others.F.s <- functional_group_order.s[!functional_group_order.s == 'Others']
functional_group_order.s <- append(no.others.F.s, "Others")
soddata$functional_group <- factor(soddata$functional_group, 
                                   levels = functional_group_order.s)

washdata <- washdata[order(-washdata$meanra),]
functional_group_order.w <- as.list(unique(washdata$functional_group))
no.others.F.w <- functional_group_order.w[!functional_group_order.w == 'Others']
functional_group_order.w <- append(no.others.F.w, "Others")
washdata$functional_group <- factor(washdata$functional_group, 
                                    levels = functional_group_order.w)

#sodfunction plot
sodplotfunction=ggplot(soddata, aes(x=paste(samp_time, treat, as.character(Row.names)),
                                    y=as.numeric(meanra),
                                    fill=functional_group))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=functioncolors.s)+
  theme(axis.text.x = element_text(angle=90))+
  facet_nested(.~samp_time, scales="free", space="free")+
  theme(axis.text.x = element_blank())+
  guides(fill=guide_legend(ncol=1))+
  ggtitle("Sod") +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Sample",
       y = "Relative abundance",
       fill = "Functional group")
sodplotfunction

#wash function plot
washplotfuction=ggplot(washdata, aes(x=paste(samp_time, treat, as.character(Row.names)),
                                     y=as.numeric(meanra),
                                     fill=functional_group))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=functioncolors.w)+
  theme(axis.text.x = element_text(angle=90))+
  facet_nested(.~samp_time, scales="free", space="free")+
  theme(axis.text.x = element_blank())+
  guides(fill=guide_legend(ncol=1))+
  ggtitle("Wash") +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Sample",
       y = "Relative abundance",
       fill = "Functional group")
washplotfuction

#plot together
ggarrange(sodplotfunction,washplotfuction, ncol=1, heights=c(1, 1))

