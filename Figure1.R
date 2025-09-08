library(treeio)
library(ggtree)
library(ape)
library(tidytree)
library(phangorn)
library(ggplot2)
library(ggtreeExtra)
library(ggnewscale)
library(tidyr)
library(reshape2)
library(dplyr)

# Load tree and metadata files
#### Figure 1 ####
pan.mavium.tree <- read.newick("~/Documents/guthrielab/m-avium/pangenome-20241020/panaroo-tree/core-genome.treefile")
metadata <- read.csv("/home/idowu/Documents/guthrielab/m-avium/pangenome-20240426-161121/iqtree/annotate.csv", header = T, na.strings = "")
metadata$Location <- as.factor(metadata$Location)
metadata$Host <- as.factor(metadata$Host)
avium.baps <- read.csv("/home/idowu/Documents/guthrielab/m-avium/pangenome-20241020/avium_hierbaps_partition1.csv", header = T)
avium.baps$level.1 <- as.factor(avium.baps$level.1)
#avium.baps$level.2 <- as.factor(avium.baps$level.2)

p1 <- ggtree(midpoint(pan.mavium.tree2), layout = "circular") %<+% metadata +
  geom_tippoint(aes(colour=Host), size=1.5) + scale_colour_manual(name = "Host", values = c("#1dcad3", "#d3261d","gold2")) + 
  geom_treescale(width=0.0005)

p2 <- p1 + geom_fruit(
  geom = geom_tile, 
  mapping = aes(fill=Location), 
  pwidth=0.0, 
  offset = 0.1, width=0.0005) +
  scale_fill_manual(name = "Location", values = c("#623cc9","#d0c5ef","#c90c1c","#063970","#abdbe3","#c9840c"))

p4 <- p2 %<+% avium.baps+new_scale_fill() + geom_fruit(
  geom = geom_tile,
  mapping = aes(fill=level.1),
  pwidth = 0.01,
  offset = 0.05, width = 0.0005) +
  scale_fill_manual(name = "Lineage", values=c("#2196f3","#212bf3","#7f21f3","#e821f3","#610d3b","#f3212b",
                                               "#f37f21","#f3e821","#94f321","#2ba421","#72be9f","#21f3e8")) +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size = 14))

ggsave("Figure1A.pdf", plot = p4, device = "pdf", width = 10, height = 8, units = "in")

# Plot pangenome per lineage
lineage.pg <- read.csv('/home/idowu/Documents/guthrielab/m-avium/pangenome-20241020/lineage-pangenome1.csv')
lineage.pg$Lineages <- as.factor(lineage.pg$Lineages)
#lineage.pg$Lineages <- factor(lineage.pg$Lineages, levels = unique(lineage.pg$Lineages))

genes <- ggplot(lineage.pg, aes(x=Lineages, y = Total.number.of.genes)) +
  geom_col(aes(fill = Pan.genome), width = 0.7) + theme_classic() + coord_flip() + 
  scale_fill_manual(name = "Pan-genome", values = c("#045080","#178ED9")) +
  ylab("Total number of genes") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
ggsave("Figure1B.pdf", plot = genes, device = "pdf", width = 10, height = 8, units = "in")

