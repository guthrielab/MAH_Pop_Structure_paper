library(treeio)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(igraph)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(phangorn)
library(ggpubr)

#### Figure 4 ####

subsampled.tree <- read.newick('core_genome_subsmp.aln.treefile')
subsampled.lineage <- read.csv('lineage_information.csv')
subsampled.lineage$Lineage <- as.factor(subsampled.lineage$Lineage)
subsampled.recomb <- read.delim('recent_recombinations.tsv', header = T, sep='\t')
subsampled.recomb$StrainName <- as.factor(subsampled.recomb$StrainName)
subsampled.recomb$DonorLineage <- as.factor(subsampled.recomb$DonorLineage)
subsampled.recomb$length <- subsampled.recomb$End - subsampled.recomb$Start
write.csv(subsampled.recomb, file = "recombination_length.csv",row.names = F)
recomb_summary <- read.csv("recomb_length_summary.csv")
recomb_summary$StrainName <- as.factor(recomb_summary$StrainName)
recomb_summary$DonorLineage <- as.factor(recomb_summary$DonorLineage)

recombtree <- ggtree(subsampled.tree) %<+% subsampled.lineage + geom_tippoint(aes(color=Lineage),size=2.5)+ 
  scale_color_manual(name = "Lineage", values=c("#2196f3","#212bf3","#7f21f3","#e821f3","#610d3b","#f3212b",
                                                "#f37f21","#f3e821","#94f321","#2ba421"))

panel4 <- recombtree + geom_facet(panel = "Recombination Origin Across Core Genomes", data = subsampled.recomb, geom = geom_segment,
                        mapping=aes(x=Start, xend=End, color=DonorLineage),linewidth=4) + 
  scale_x_continuous(labels = label_comma()) +
  geom_facet(panel = "Proportion of Recent Recombination Events", data = recomb_summary, geom = geom_col,
             aes(x = Proportion, fill=DonorLineage), orientation = 'y', width = .6)+
  scale_fill_manual(name = "Lineage", values=c("#2196f3","#212bf3","#7f21f3","#e821f3","#610d3b","#f3212b",
                                               "#f37f21","#f3e821","#94f321","#2ba421")) + theme_tree2(legend.position=c(.05, .75)) +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14)) + labs(x="Genome Position")

ggsave("Figure4.pdf", plot = panel4, device = "pdf", width = 15, height = 10, units = "in")

recomb100 <- subsampled.proportion %>% right_join(subsampled.lineage, by = c('StrainName' = 'Isolate'))
#Total count of lineage 10
sum(recomb100$Count[recomb100$Lineage == 10], na.rm = TRUE)
#Total count of lineage 9
sum(recomb100$Count[recomb100$Lineage == 9], na.rm = TRUE)
#Total count of lineage 8
sum(recomb100$Count[recomb100$Lineage == 8], na.rm = TRUE)
#Total count of lineage 7
sum(recomb100$Count[recomb100$Lineage == 7], na.rm = TRUE)
#Total count of lineage 6
sum(recomb100$Count[recomb100$Lineage == 6], na.rm = TRUE)
#Total count of lineage 5
sum(recomb100$Count[recomb100$Lineage == 5], na.rm = TRUE)
#Total count of lineage 4
sum(recomb100$Count[recomb100$Lineage == 4], na.rm = TRUE)
#Total count of lineage 3
sum(recomb100$Count[recomb100$Lineage == 3], na.rm = TRUE)
#Total count of lineage 2
sum(recomb100$Count[recomb100$Lineage == 2], na.rm = TRUE)
#Total count of lineage 1
sum(recomb100$Count[recomb100$Lineage == 1], na.rm = TRUE)
write.csv(recomb100, file = 'recomb100.csv')

#### Supplementary Figure S4 ####

# Pearson Correlation Coefficient

contigs.recomb <- read.csv("recomb_contigs.csv")

ggplot(contigs.recomb, aes(x=number_of_contigs, y=recombination_events)) +
  geom_point(size=2, color="purple") + theme_classic() +
  geom_smooth(method = "lm") +
  xlab("Number of contigs") +
  ylab("Recombination events")
cor.test(contigs.recomb$recombination_events,contigs.recomb$number_of_contigs,
         method = 'pearson')
s4 <- ggscatter(contigs.recomb,x="number_of_contigs", y="recombination_events",
          add = "reg.line", conf.int = T,
          cor.coef = T, cor.method = "spearman",
          xlab="Number of contigs", ylab="Recombination events")

s4.1 <- ggscatter(contigs.recomb, x = "number_of_contigs", y = "recombination_events", add = "reg.line") +
  stat_cor(method = "spearman", label.x = 3)

p_adjusted <- s4.1 %>% ggadjust_pvalue(p.adjust.method = "bonferroni")

shapiro.test(contigs.recomb$recombination_events)
shapiro.test(contigs.recomb$number_of_contigs)


cor_test_result <- cor.test(contigs.recomb$number_of_contigs, contigs.recomb$recombination_events, method = "spearman", exact = FALSE)
adjusted_p <- p.adjust(cor_test_result$p.value, method = "bonferroni")
SupFig4 <- ggscatter(contigs.recomb, x = "number_of_contigs", y = "recombination_events",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Number of contigs", ylab = "Recombination events") +
  labs(subtitle = paste("Bonferroni-adjusted p-value:", 
                        format(adjusted_p, scientific = TRUE, digits = 3)))

ggsave("SuppFig4.pdf", plot = SupFig4)

#### Supplementary Figure S5 ####

# mmpL4 recombination
mmpl4.tree <- read.newick("output/rnd.fasta.treefile")
big.tree <- read.newick('../../pangenome-20241020/panaroo-tree/core-genome.treefile')
mmpl4.recomb <- read.csv("output/recombinations_recent1.txt")
mmpl4.lin <- read.csv("../../pangenome-20241020/panaroo-tree/figtreeannot1.tsv.csv")
mmpl4.recomb$StrainName <- as.factor(mmpl4.recomb$StrainName)
mmpl4.recomb$DonorLineage <- as.factor(mmpl4.recomb$DonorLineage)
mmpl4.lin$lineage <- as.factor(mmpl4.lin$lineage)

rnd.tree <- ggtree(midpoint(big.tree)) %<+% mmpl4.lin + geom_tippoint(aes(color=lineage),size=2.5)+ 
  scale_color_manual(name = "Lineage", values=c("#2196f3","#212bf3","#7f21f3","#e821f3","#610d3b","#f3212b",
                                                "#f37f21","#f3e821","#94f321","#2ba421"))

S5 <- rnd.tree + geom_facet(panel = "Recombination Origins in the MmpL4 protein", data = mmpl4.recomb, geom = geom_segment,
                        mapping=aes(x=Start, xend=End, color=DonorLineage),linewidth=1) +
  scale_x_continuous(labels = label_comma()) +
  scale_fill_manual(name = "Lineage", values=c("#2196f3","#212bf3","#7f21f3","#e821f3","#610d3b","#f3212b",
                                               "#f37f21","#f3e821","#94f321","#2ba421")) + theme_tree2(legend.position=c(.05, .75)) +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14,),
        plot.margin = margin(t=20, b=20, l=20, r=20, unit = "pt")) + labs(x="Gene Length")

ggsave("SuppFig5.pdf", plot = S5, width = 12, height = 8)
