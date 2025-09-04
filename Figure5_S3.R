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

##All Isolates##

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
