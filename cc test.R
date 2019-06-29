# load libraries
library(phyloseq)
library(ggplot2) 
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(knitr)

# save session info (packages and versions loaded)
session <- sessionInfo()
# Create function to plot bar plots with contaminants in grey scale and expected mock microbial sequences in color

expCompBarPlot <- function(physeq, exp_taxa, title){
  ## physeq - phyloseq object that will be plotted
  ## exp_taxa - taxa that are expected to be in the mock community 
  ## title - title for plot
  #set up data_table
  data_table <- as.data.frame(t(physeq@otu_table))
  
  data_table$reference = FALSE
  data_table$reference[rownames(data_table) %in% exp_taxa] = TRUE
  sample_names <- sample_names(physeq)
  data_table$id <- paste0('ASV_', 1:nrow(data_table))
  dilution_labels <- sample_data(physeq)$Dilutions
  
  set.seed(444)
  
  # define the colors to use for reference and non-reference OTUs/ASVs
  ref_colors <- brewer.pal(sum(data_table$reference), "Paired")
  other_colors <- sample(grey.colors(5, start = 0.5, end = 0.9), sum(!data_table$reference), replace = TRUE)
  
  # add a color variable to the data table
  data_table$color <- rep(NA, nrow(data_table))
  data_table$color[data_table$reference] <- ref_colors
  data_table$color[!data_table$reference] <- other_colors
  
  # reshape the data table into a ggplot-friendly format, by gathering samples into a single column called "count"
  
  color_gg <- data_table %>% select(id, sample_names, color) %>% gather("sample", "count", sample_names)
  legend_color <- c(bright = ref_colors[2], dull = other_colors[2])
  data_gg <- data_table %>% gather("sample", "count", sample_names)
  
  data_gg <- inner_join(data_gg,color_gg)
  
  # create the composition bar plot
  comp_bar <- ggplot(data_gg, aes(x = sample, y = count)) +
    geom_col(aes(fill = color, group = reference, alpha = ifelse(reference, "bright", "dull")), width = 0.7, position = position_fill()) +
    scale_fill_identity(guide = FALSE) +
    scale_alpha_manual(name = "Sequence type",
                       labels = c("expected sequences", "other"),
                       values = c(bright = 1, dull = 1),
                       guide = guide_legend(override.aes = list(fill = c(ref_colors[4], "#AEAEAE")),
                                            keywidth = NULL, keyheight = NULL)) +
    labs(title = title, x = "sample", y = "Relative Abundance") +
    theme(legend.position = "right", legend.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 16))
  comp_bar
}
## Load the dataset
load("mockDilutions.RData")
# Create profile of only expected sequences from the undiluted mock microbial communtiy sample

# subset the undiluted mock microbial sample (sample name 'D0')
mock_ps_pure<-subset_samples(mock_ps,sample_names(mock_ps)== 'D0')

# remove ASVs that are not present in the undiluted sample
mock_ps_pure<-prune_taxa(taxa_sums(mock_ps_pure)>0,mock_ps_pure)

# change the SampleType and sample_name of the pure mock microbial community sample 
sample_data(mock_ps_pure)$SampleType <-'MockCommunityProfile'
sample_names(mock_ps_pure) <-paste('mc',sample_names(mock_ps_pure),sep = '_')

# display a summary of the new phyloseq object
mock_ps_pure
# Remove the unexpected ASVs from the undiluted mock microbial community dilution series

# make a list of the top 9 abundant ASV taxa names (this is plausible for filtering since the 9 sequences we want to remove are present in low abundance)
mock_taxa = names(sort(taxa_sums(mock_ps_pure), decreasing = TRUE)[1:9])

# subset the taxa in mock_ps_pure so only the expected sequences are present
mock_ps_pure<-prune_taxa(mock_taxa,mock_ps_pure)
# display a summary of the mock microbial dilution series phyloseq object
mock_ps
# create a phyloseq object that is normalized to 100 (relative abundance)
ps_norm<-transform_sample_counts(ps,function(x) 100* x/sum(x))
mock_ps_norm <- transform_sample_counts(mock_ps,function(x) 100* x/sum(x))

# Identify the proportion of each sample that is the expected mock microbial ASVs
ps_norm_exp<-prune_taxa(mock_taxa,ps_norm)

# Create a table with the dilution, number of reads per sample, and proportion of contaminants per sample
dilutionSummary <- data.frame(DilutionSeries = sample_names(ps),NumberOfReads = sample_sums(ps), PercentContaminants = 100-sample_sums(ps_norm_exp))

# Create a variable to indicate the sample order of the plots
dilutions<-c('D0','D1','D2','D3','D4','D5','D6','D7','D8', 'Blank')

# Create plots to summarize these data
## Plot Figure 1A - number of reads per sample across dilution series
ggplot(dilutionSummary, aes(x = DilutionSeries, y = NumberOfReads)) + geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() + scale_x_discrete(limits = dilutions)

## Plot Figure 1B - Percent of contaminants across dilution series
ggplot(dilutionSummary, aes(x = DilutionSeries, y = PercentContaminants)) + geom_point() + scale_x_discrete(limits = dilutions)

## Plot Figure 1C - Stacked bar plot of Mock microbial dilution series
expCompBarPlot(ps_norm, mock_taxa, 'Initial Mock Microbial Community Dilution')  + scale_x_discrete(limits = dilutions)
# create a list of unexpected sequences (contaminants)
# create a list of all ASV taxa names
contaminant_taxa<-taxa_names(mock_ps)
# remove the expected mock microbial community ASV taxa names
contaminant_taxa <- contaminant_taxa[!(contaminant_taxa %in% mock_taxa)]

# create a phyloseq object that only contains the contaminant sequences
contaminants_ps<-prune_taxa(contaminant_taxa,mock_ps)
contaminants_ps<- prune_taxa(taxa_sums(contaminants_ps)>0,contaminants_ps)

# change the sample names to indicate that these samples only contain contmaminant ASVs
sample_names(contaminants_ps)<-paste('con',sample_names(contaminants_ps),sep = '_')
sample_data(contaminants_ps)$SampleType<-'ContaminantProfile'
# Number of ASVs in common between contaminants in mock community and the blank control
print(paste('Total number of contaminant ASVs', length(taxa_names(contaminants_ps))))
print(paste('Number of contaminant ASVs also present in blank', length(intersect(taxa_names(contaminants_ps),taxa_names(blank_ps)))))

# create a list of contaminants taxa that are not present in the blank control
contaminant_taxa_no_blank<-taxa_names(contaminants_ps)
contaminant_taxa_no_blank <- contaminant_taxa_no_blank[!(contaminant_taxa_no_blank %in% taxa_names(blank_ps))]

# Create  a binary list of contaminant ASVs indicating if the ASV is present in the blank control (1) or not (0)
contaminants_in_blank <- data.frame(matrix(1, ncol = length(taxa_names(contaminants_ps)), nrow = 1))
colnames(contaminants_in_blank) <- taxa_names(contaminants_ps)
contaminants_in_blank[,contaminant_taxa_no_blank] <- 0 
contaminants_in_blank <- t(contaminants_in_blank)
# Identify the contribution per sample of contaminants that are not present in blanks
# generate a phyloseq object with contaminants only normalized to 100 
contaminant_ps_norm <- transform_sample_counts(contaminants_ps,function(x) 100* x/sum(x))
contaminant_no_blanks<-prune_taxa(contaminant_taxa_no_blank,contaminant_ps_norm)

# Plot the proportion of contaminant ASVs per sample that were not present in the blank control
plot_bar(contaminant_no_blanks,fill='Genus',, title = ' Proportion of Contaminant ASVs Not in the Blank Control Sample') + theme(legend.position='none') + ylim(c(0,100))

# sum the amount of contaminant signal not arising from ASVs in blank control
100 - sample_sums(contaminant_no_blanks)
summary(100 - sample_sums(contaminant_no_blanks))
# Count number of contaminants present in only one sample
contaminant_bin<-as.data.frame(contaminants_ps@otu_table)
contaminant_bin[contaminant_bin>0]<-1
contaminant_bin = t(contaminant_bin)
contaminant_nsamples <- rowSums(contaminant_bin)
table(contaminant_nsamples)
# identify the maximum proportion of reads for each contaminant ASV
contaminant_max <- apply(as.data.frame(mock_ps_norm@otu_table), 2, max)
contaminant_max <- subset(contaminant_max, names(contaminant_max) %in% contaminant_taxa)

# identify the minimum proportion of reads for each contaminant ASV
contaminant_min <- apply(as.data.frame(mock_ps_norm@otu_table), 2, min)
contaminant_min <- subset(contaminant_min, names(contaminant_min) %in% contaminant_taxa)

# summarize in table (Supplemental Table 2)
contaminantSummary <- cbind( ASV = names(contaminant_max),as.data.frame(contaminants_ps@tax_table), Maximum = contaminant_max, Minimum = contaminant_min, NumberOfSamples = contaminant_nsamples, InBlank =  contaminants_in_blank)
mock_alpha <- estimate_richness(mock_ps, measures=c("Observed", "Shannon", "InvSimpson"))
mock_alpha_pure <- estimate_richness(mock_ps_pure,  measures=c("Observed", "Shannon", "InvSimpson"))
plot_richness(mock_ps, measure = c('Observed','Shannon', 'InvSimpson') )

# calculate the difference between the observed and actual alpha diversity measures 
max_diff <- ((mock_alpha[9,] - mock_alpha_pure) / mock_alpha_pure) 
# Save workspace
save.image("mockDilutionsPrep.RData")