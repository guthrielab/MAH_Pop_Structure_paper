# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(ape)
# Read the Gubbins GFF file (adjust path as needed)
gff_data <- read.gff("/home/idowu/R/x86_64-pc-linux-gnu-library/4.2/RCandy/extdata/core-snp.recombination_predictions.gff")
gff_human <- read.gff("/home/idowu/R/x86_64-pc-linux-gnu-library/4.2/RCandy/extdata/mav_humans.recombination_predictions.gff")
gff_pig <- read.gff("/home/idowu/R/x86_64-pc-linux-gnu-library/4.2/RCandy/extdata/mav_pigs.recombination_predictions.gff")
gff_env <- read.gff("/home/idowu/R/x86_64-pc-linux-gnu-library/4.2/RCandy/extdata/mav_env.recombination_predictions.gff")
gff_map <- read.gff("/home/idowu/Documents/guthrielab/MAP/map.recombination_predictions.gff")

####All Isolates####

# Extract only CDS features
cds_data <- gff_data[gff_data$type == "CDS", ]

# Calculate the length of each recombination event
cds_data$length <- cds_data$end - cds_data$start + 1

# Extract taxa information from attributes
cds_data$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', cds_data$attributes)

# Create a function to count recombination events in windows
calculate_recombination_density <- function(data, window_size = 1000) {
  genome_length <- max(data$end)
  windows <- seq(1, genome_length, by = window_size)
  
  # Count recombination events in each window
  density_df <- data.frame(
    position = windows,
    count = sapply(windows, function(x) {
      sum(data$start <= x + window_size - 1 & data$end >= x)
    }),
    density = sapply(windows, function(x) {
      sum(data$start <= x + window_size - 1 & data$end >= x) / (window_size / 1000)
    })
  )
  
  return(density_df)
}

# Calculate recombination density
recomb_density <- calculate_recombination_density(cds_data)

# Calculate the 99.9th percentile (top 0.1%) threshold
hotspot_threshold <- quantile(recomb_density$density, probs = 0.999, na.rm = TRUE)

# Plot recombination hotspots with threshold line
a <- ggplot(recomb_density, aes(x = position, y = density)) +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = hotspot_threshold, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = (max(recomb_density$position)/2)+5, y = hotspot_threshold * 1.1, 
           label = paste("Top 0.1% threshold =", round(hotspot_threshold, 2)), 
           color = "red") +
  labs(title = "SNP Density in All Samples",
       x = "Genome Position",
       y = "SNP Density per kb") +
  theme_classic() + scale_x_continuous(labels = unit_format(unit = "Mbp", scale = 1e-6)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.margin = margin(t=20, b=20, l=20, r=20, unit = "pt")) + ylim(0,250)

##Human Isolates##

# Extract only CDS features
cds_human <- gff_human[gff_human$type == "CDS", ]

# Calculate the length of each recombination event
cds_human$length <- cds_human$end - cds_human$start + 1

# Extract taxa information from attributes
cds_human$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', cds_human$attributes)

# Calculate recombination density
recomb_density_human <- calculate_recombination_density(cds_human)

# Calculate the 99.9th percentile (top 0.1%) threshold
hotspot_threshold_human <- quantile(recomb_density_human$density, probs = 0.999, na.rm = TRUE)

# Plot recombination hotspots with threshold line
b <- ggplot(recomb_density_human, aes(x = position, y = density)) +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = hotspot_threshold_human, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = (max(recomb_density_human$position)/2)+5, y = hotspot_threshold_human * 1.10, 
           label = paste("Top 0.1% threshold =", round(hotspot_threshold_human, 2)), 
           color = "red") +
  labs(title = "SNP Density in Human Samples",
       x = "Genome Position",
       y = "SNP Density per kb") +
  theme_classic() + scale_x_continuous(labels = unit_format(unit = "Mbp", scale = 1e-6)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.margin = margin(t=20, b=20, l=20, r=20, unit = "pt")) + ylim(0,250)

##Pig Isolates##

# Extract only CDS features
cds_pig <- gff_pig[gff_pig$type == "CDS", ]

# Calculate the length of each recombination event
cds_pig$length <- cds_pig$end - cds_pig$start + 1

# Extract taxa information from attributes
cds_pig$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', cds_pig$attributes)

# Calculate recombination density
recomb_density_pig <- calculate_recombination_density(cds_pig)

# Calculate the 99.9th percentile (top 0.1%) threshold
hotspot_threshold_pig <- quantile(recomb_density_pig$density, probs = 0.999, na.rm = TRUE)

# Plot recombination hotspots with threshold line
c <- ggplot(recomb_density_pig, aes(x = position, y = density)) +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = hotspot_threshold_pig, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = (max(recomb_density_pig$position)/2)+5, y = hotspot_threshold_pig * 1.10, 
           label = paste("Top 0.1% threshold =", round(hotspot_threshold_pig, 2)), 
           color = "red") +
  labs(title = "SNP Density in Pig Samples",
       x = "Genome Position",
       y = "SNP Density per kb") +
  theme_classic() + scale_x_continuous(labels = unit_format(unit = "Mbp", scale = 1e-6)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.margin = margin(t=20, b=20, l=20, r=20, unit = "pt")) + ylim(0, 20)

##Environmental Isolates##

# Extract only CDS features
cds_env <- gff_env[gff_env$type == "CDS", ]

# Calculate the length of each recombination event
cds_env$length <- cds_env$end - cds_env$start + 1

# Extract taxa information from attributes
cds_env$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', cds_env$attributes)

# Calculate recombination density
recomb_density_env <- calculate_recombination_density(cds_env)

# Calculate the 99.9th percentile (top 0.1%) threshold
hotspot_threshold_env <- quantile(recomb_density_env$density, probs = 0.999, na.rm = TRUE)

# Plot recombination hotspots with threshold line
d <- ggplot(recomb_density_env, aes(x = position, y = density)) +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = hotspot_threshold_env, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = (max(recomb_density_env$position)/2)+5, y = hotspot_threshold_env * 1.10, 
           label = paste("Top 0.1% threshold =", round(hotspot_threshold_env, 2)), 
           color = "red") +
  labs(title = "SNP Density in Environmental Samples",
       x = "Genome Position",
       y = "SNP Density per kb") +
  theme_classic() + scale_x_continuous(labels = unit_format(unit = "Mbp", scale = 1e-6)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.margin = margin(t=20, b=20, l=20, r=20, unit = "pt")) + ylim(0,20)

final.plot <- cowplot::plot_grid(a,b,c,d, nrow=2, ncol = 2)
ggsave("/home/idowu/Documents/guthrielab/m-avium/fastq/mav-snp-density-hotspots-kb1.pdf", width = 10, height = 8)




##MAP Isolates##

# Extract only CDS features
gff_map <- gff_map[gff_map$type == "CDS", ]

# Calculate the length of each recombination event
gff_map$length <- gff_map$end - gff_map$start + 1

# Extract taxa information from attributes
gff_map$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', gff_map$attributes)

# Calculate recombination density
gff_map_dens <- calculate_recombination_density(gff_map)

# Calculate the 99.9th percentile (top 0.1%) threshold
#hotspot_threshold_env <- quantile(recomb_density_env$density, probs = 0.999, na.rm = TRUE)

# Plot recombination hotspots with threshold line
ggplot(gff_map_dens, aes(x = position, y = density)) +
  geom_line(color = "steelblue") +
  labs(title = "SNP Density in MAP Isolates",
       x = "Genome Position (bp)",
       y = "SNP Density per kb") +
  theme_classic() + scale_x_continuous(labels = label_comma()) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
write.csv(gff_map_dens, file = "/home/idowu/Documents/guthrielab/MAP/map_snp_dens.csv",row.names = F)
write.csv(recomb_density, file = "/home/idowu/Documents/guthrielab/m-avium/fastq/mav_snp_dens.csv", row.names = F)

mah.map <- read.csv("/home/idowu/Documents/guthrielab/m-avium/fastq/mah-map-snp-density.csv")

## Plot Boxplot to compare MAH and MAP SNP densities

p <- ggplot(mah.map, aes(x=Subsp, y=SNP_density)) +
  geom_boxplot(width=0.5, color='black') +
  theme_classic() + theme(legend.position = "none") + xlab("M. avium subsp.") +
  ylab("SNP Density per kb") + geom_signif(comparisons = list(c("MAH", "MAP")),
                                           map_signif_level = T)
pp <- p + stat_compare_means()
ggsave("SuppFig3.pdf", plot = pp, width = 10, height = 8)

wilcox.test(mah.map$SNP_density[mah.map$Subsp == "MAH"], mah.map$SNP_density[mah.map$Subsp == "MAP"])
summary(mah.map$SNP_density[mah.map$Subsp == "MAH"])
summary(mah.map$SNP_density[mah.map$Subsp == "MAP"])


####Lineage Recombination####

setwd("/Users/idolawoye/Downloads/3_projects/avium/lineage_recomb/")
library(ape)

##lineage1##
lin1 <- read.table("lin_1.recombination_predictions.gff", sep = "\t", header = F,
                                 comment.char = "#", stringsAsFactors = FALSE)
colnames(lin1) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
lin1 <- read.gff("lin_1.recombination_predictions.gff")

# Create a recombination events per base dataset
lin1 <- lin1 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin1$lineage <- "lineage_1"

##lineage2##
lin2 <- read.table("lin_2.recombination_predictions.gff", sep = "\t", header = F,
                   comment.char = "#", stringsAsFactors = FALSE)
colnames(lin2) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Create a recombination events per base dataset
lin2 <- lin2 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin2$lineage <- "lineage_2"

##lineage3##
lin3 <- read.table("lin_3.recombination_predictions.gff", sep = "\t", header = F,
                   comment.char = "#", stringsAsFactors = FALSE)
colnames(lin3) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Create a recombination events per base dataset
lin3 <- lin3 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin3$lineage <- "lineage_3"

##lineage4##
lin4 <- read.table("lin_4.recombination_predictions.gff", sep = "\t", header = F,
                   comment.char = "#", stringsAsFactors = FALSE)
colnames(lin4) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Create a recombination events per base dataset
lin4 <- lin4 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin4$lineage <- "lineage_4"

##lineage5##
lin5 <- read.table("lin_5.recombination_predictions.gff", sep = "\t", header = F,
                   comment.char = "#", stringsAsFactors = FALSE)
colnames(lin5) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Create a recombination events per base dataset
lin5 <- lin5 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin5$lineage <- "lineage_5"

##lineage6##
lin6 <- read.table("lin_6.recombination_predictions.gff", sep = "\t", header = F,
                   comment.char = "#", stringsAsFactors = FALSE)
colnames(lin6) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Create a recombination events per base dataset
lin6 <- lin6 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin6$lineage <- "lineage_6"

##lineage7##
lin7 <- read.table("lin_7.recombination_predictions.gff", sep = "\t", header = F,
                   comment.char = "#", stringsAsFactors = FALSE)
colnames(lin7) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Create a recombination events per base dataset
lin7 <- lin7 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin7$lineage <- "lineage_7"

##lineage8##
lin8 <- read.table("lin_8.recombination_predictions.gff", sep = "\t", header = F,
                   comment.char = "#", stringsAsFactors = FALSE)
colnames(lin8) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Create a recombination events per base dataset
lin8 <- lin8 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin8$lineage <- "lineage_8"

##lineage9##
lin9 <- read.table("lin_9.recombination_predictions.gff", sep = "\t", header = F,
                   comment.char = "#", stringsAsFactors = FALSE)
colnames(lin9) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Create a recombination events per base dataset
lin9 <- lin9 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin9$lineage <- "lineage_9"

##lineage10##
lin10 <- read.table("lin_10.recombination_predictions.gff", sep = "\t", header = F,
                   comment.char = "#", stringsAsFactors = FALSE)
colnames(lin10) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Create a recombination events per base dataset
lin10 <- lin10 %>%
  rowwise() %>%
  mutate(position = list(seq(start, end))) %>%
  unnest(position) %>%
  group_by(position) %>%
  summarise(events = n()) %>%
  ungroup()

lin10$lineage <- "lineage_10"

library(tidyverse)
all_lineages <- bind_rows(lin1,lin2,lin3,lin4,lin5,lin6,lin7,lin8,lin9,lin10)
all_lineages$lineage <- factor(all_lineages$lineage, levels = unique (all_lineages$lineage))

lin1.10 <- ggplot(all_lineages, aes(x=position, y=events, group=lineage, color=lineage)) +
  scale_color_manual(name = "Lineage", values = c("#2196f3","#212bf3","#7f21f3","#e821f3",
                                                  "#610d3b","#f3212b","#f37f21","#f3e831",
                                                  "#94f321","#2ba421","#72be9f","#21f3e8")) +
  geom_line() + facet_wrap(~ lineage) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1), legend.position = "none") + 
  scale_x_continuous(labels = label_comma()) +
  labs(x="Genomic Position", y="Recombination Events")

ggsave("all_lineages_recomb.pdf",plot = lin1.10, height=8, width = 10) 


all_recomb$region_id <- gsub(".*ID=|;.*", "", all_recomb$attributes)
all_recomb$y <- 1
ggplot(all_recomb, aes(xmin=start, xmax=end, ymin=0, ymax=1)) +
  geom_line(fill="black", alpha=0.7) +
  scale_y_continuous(expand = c(0,0)) 


lineage1 <- read.gff("lin_1.recombination_predictions.gff")
lineage1 <- lineage1[lineage1$type == 'CDS', ]
lineage1$length <- lineage1$end - lineage1$start + 1
lineage1$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage1$attributes)

lineage1.dens <- calculate_recombination_density(lineage1)
lineage1.dens$lineage <- "lineage_1"

ggplot(lineage1.dens, aes(x = position, y = density)) +
  geom_line(color = "#2196f3")

lineage2 <- read.gff("lin_2.recombination_predictions.gff")

lineage2 <- lineage2[lineage2$type == 'CDS', ]
lineage2$length <- lineage2$end - lineage2$start + 1
lineage2$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage2$attributes)
lineage2.dens <- calculate_recombination_density(lineage2)
lineage2.dens$lineage <- "lineage_2"

lineage3 <- read.gff("lin_3.recombination_predictions.gff")

lineage3 <- lineage3[lineage3$type == 'CDS', ]
lineage3$length <- lineage3$end - lineage3$start + 1
lineage3$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage3$attributes)
lineage3.dens <- calculate_recombination_density(lineage3)
lineage3.dens$lineage <- "lineage_3"

lineage4 <- read.gff("lin_4.recombination_predictions.gff")

lineage4 <- lineage4[lineage4$type == 'CDS', ]
lineage4$length <- lineage4$end - lineage4$start + 1
lineage4$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage4$attributes)
lineage4.dens <- calculate_recombination_density(lineage4)
lineage4.dens$lineage <- "lineage_4"

lineage5 <- read.gff("lin_5.recombination_predictions.gff")

lineage5 <- lineage5[lineage5$type == 'CDS', ]
lineage5$length <- lineage5$end - lineage5$start + 1
lineage5$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage5$attributes)
lineage5.dens <- calculate_recombination_density(lineage5)
lineage5.dens$lineage <- "lineage_5"

lineage6 <- read.gff("lin_6.recombination_predictions.gff")

lineage6 <- lineage6[lineage6$type == 'CDS', ]
lineage6$length <- lineage6$end - lineage6$start + 1
lineage6$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage6$attributes)
lineage6.dens <- calculate_recombination_density(lineage6)
lineage6.dens$lineage <- "lineage_6"

lineage7 <- read.gff("lin_7.recombination_predictions.gff")

lineage7 <- lineage7[lineage7$type == 'CDS', ]
lineage7$length <- lineage7$end - lineage7$start + 1
lineage7$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage7$attributes)
lineage7.dens <- calculate_recombination_density(lineage7)
lineage7.dens$lineage <- "lineage_7"

lineage8 <- read.gff("lin_8.recombination_predictions.gff")

lineage8 <- lineage8[lineage8$type == 'CDS', ]
lineage8$length <- lineage8$end - lineage8$start + 1
lineage8$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage8$attributes)
lineage8.dens <- calculate_recombination_density(lineage8)
lineage8.dens$lineage <- "lineage_8"

lineage9 <- read.gff("lin_9.recombination_predictions.gff")

lineage9 <- lineage9[lineage9$type == 'CDS', ]
lineage9$length <- lineage9$end - lineage9$start + 1
lineage9$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage9$attributes)
lineage9.dens <- calculate_recombination_density(lineage9)
lineage9.dens$lineage <- "lineage_9"

lineage10 <- read.gff("lin_10.recombination_predictions.gff")

lineage10 <- lineage10[lineage10$type == 'CDS', ]
lineage10$length <- lineage10$end - lineage10$start + 1
lineage10$taxa <- gsub('.*taxa="([^"]+)".*', '\\1', lineage10$attributes)
lineage10.dens <- calculate_recombination_density(lineage10)
lineage10.dens$lineage <- "lineage_10"

merge_lineages <- bind_rows(lineage1.dens,lineage2.dens,lineage3.dens,lineage4.dens,lineage5.dens,
                            lineage6.dens,lineage7.dens,lineage8.dens,lineage9.dens,lineage10.dens)
sc_lineages <- merge_lineages

merge_lineages$lineage <- factor(merge_lineages$lineage, levels = unique (merge_lineages$lineage))

lineages1.10 <- ggplot(merge_lineages, aes(x=position, y=density, group=lineage, color=lineage)) +
  scale_color_manual(name = "Lineage", values = c("#2196f3","#212bf3","#7f21f3","#e821f3",
                                                  "#610d3b","#f3212b","#f37f21","#f3e831",
                                                  "#94f321","#2ba421","#72be9f","#21f3e8")) +
  geom_line() + facet_wrap(~ lineage) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1), legend.position = "none") + 
  scale_x_continuous(labels = label_comma()) +
  labs(x="Genomic Position", y="SNP Density per kb")

ggsave("all_lineages_snp_density.pdf",plot = lineages1.10, height=8, width = 10)


sc_lineages$lineage <- gsub('lineage_1', 'SC5', sc_lineages$lineage)
sc_lineages$lineage <- gsub('lineage_2', 'SC1', sc_lineages$lineage)
sc_lineages$lineage <- gsub('lineage_3', 'MahEA2', sc_lineages$lineage)
sc_lineages$lineage <- gsub('lineage_4', 'SC2', sc_lineages$lineage)
sc_lineages$lineage <- gsub('lineage_5', 'SC3', sc_lineages$lineage)
sc_lineages$lineage <- gsub('lineage_6', 'SC6', sc_lineages$lineage)
sc_lineages$lineage <- gsub('lineage_7', 'SC4', sc_lineages$lineage)
sc_lineages$lineage <- gsub('lineage_8', 'MahEA1', sc_lineages$lineage)
sc_lineages$lineage <- gsub('lineage_9', 'SC7', sc_lineages$lineage)
sc_lineages$lineage <- gsub('lineage_10', 'SC8', sc_lineages$lineage)

#sc_lineages$lineage <- as.factor(sc_lineages$lineage)

fig3 <- ggplot(sc_lineages, aes(x=position, y=density, group=lineage, color=lineage)) +
  scale_color_manual(name = "Lineage", values = c("#2196f3","#212bf3","#7f21f3","#e821f3",
                                                  "#610d3b","#f3212b","#f37f21","#f3e831",
                                                  "#94f321","#2ba421","#72be9f","#21f3e8")) +
  geom_line() + facet_wrap(~ lineage, ncol = 5) + theme_classic() +
  theme(axis.text.x = element_text(size=10, angle = 45,vjust = 1, hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 12), legend.position = "none") + 
  scale_x_continuous(labels = unit_format(unit = "Mbp", scale = 1e-6)) +
  labs(x="Genomic Position", y="SNP Density per kb")

ggsave("Figure3.pdf", plot = fig3, width = 15, height = 8)
