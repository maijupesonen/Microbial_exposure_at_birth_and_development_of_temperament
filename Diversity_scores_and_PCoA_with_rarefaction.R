
#################################################################
# Microbial exposure at birth and the development of behavioral 
# temperament during the first three years of childhood
# 
# This R-script preprocesses the ASV data, computes alpha and
# beta diversity, and computes principal components and coordinates 
# for alpha and beta diversity
#################################################################


# libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(viridis)
library(haven)
library(phytools)
library(microbiome)
library(phyloseq)
library(picante)
library(GUniFrac)


# phyloseq object
load("~/Data/bact_phyloseq_nocontrols.Rdata")

min(sample_sums(bact.asv.nocont)) # 1002
max(sample_sums(bact.asv.nocont)) # 195565

# questionnaire data
quest.dat <- read_sav("~/Data/Questionnaire_data/OMP3_selection.sav")

# subset of samples with questionnaire data

# sample names from the OTU data
otu.dat <- otu_table(bact.asv.nocont)
# sample names
colnames(otu.dat)

# sample data with KuBiCo id
sample.data <- data.frame(id = colnames(otu.dat), short.id = substr(colnames(otu.dat), 4, nchar(colnames(otu.dat))))

# add the mode of delivery variable and subset only vaginal births (0 = vaginal, 1 = CS)
mode <- quest.dat[,c("id", "Binary_Mode_delivery")]
sample.data.all <- left_join(sample.data, mode, by = "short.id")

# add sample data to the phyloseq object
rownames(sample.data.all) <- sample.data.all$id
sample_data(bact.asv.nocont) <- sample.data.all

# subset with questionnaire data
phy.quest <- subset_samples(bact.asv.nocont, (short.id %in% quest.dat$id) & Binary_Mode_delivery == 0)

# select only bacteria
phy.bacteria <- subset_taxa(phy.quest, Rank1 == "Bacteria")

# phyla.data
phy.phylum <- tax_glom(phy.bacteria, "Rank2")
tax.phylum <- tax_table(phy.phylum)

# genera.data
phy.genera <- tax_glom(phy.bacteria, "Rank6")
tax.genera <- tax_table(phy.genera)
# filter out genera with zero relative abundance
phy.genera.no.zero.rel <- filter_taxa(phy.genera, function(x){mean(x/sample_sums(phy.genera)) > 0}, prune = T)


# relative abundances for phyla and genera data
phy.phylum.rel <- transform_sample_counts(phy.phylum, function(x) x/sum(x))
phy.phylum.rel.df <- psmelt(phy.phylum.rel)
phy.genera.rel <- transform_sample_counts(phy.genera, function(x) x/sum(x))
phy.genera.rel.df <- psmelt(phy.genera.rel)


######### Most prevalent (median relative abundance > 0) genera and phyla

top.genera <- phy.genera.rel.df %>%
  filter(Rank6 != "NA") %>%
  group_by(Rank6) %>%
  mutate(median.abund.genera = median(Abundance)) %>%
  ungroup() %>%
  filter(median.abund.genera > 0) %>%
  arrange(desc(median.abund.genera)) %>%
  distinct(Rank6, median.abund.genera)


top.genera.long <- phy.genera.rel.df %>%
  filter(Rank6 %in% top.genera$Rank6) %>%
  select(short.id, Rank6, Abundance)

names(top.genera.long) <- c("short.id", "genera", "rel.abundance")
top.genera.long$genera <- factor(top.genera.long$genera, levels = c("Lactobacillus", "Staphylococcus", "Corynebacterium", "Prevotella", 
                                                                    "Anaerococcus", "Bacteroides", "Stenotrophomonas", "Streptococcus"))

# wide format, to be merged with questionnaire data
top.genera.wide <- pivot_wider(top.genera.long, names_from = genera, values_from = rel.abundance)

# top phyla
top.phyla <- phy.phylum.rel.df %>%
  filter(Rank2 != "NA") %>%
  group_by(Rank2) %>%
  mutate(median.abund.phyla = median(Abundance)) %>%
  ungroup() %>%
  filter(median.abund.phyla > 0) %>%
  arrange(desc(median.abund.phyla)) %>%
  distinct(Rank2, median.abund.phyla)


# dominance in sample, Supplementary Figure S2
top.8.most.prevalent <- top.genera$Rank6

genera.abundant <- phy.genera.rel.df %>%
                    filter(Rank6 %in% top.8.most.prevalent) %>%
                    arrange(desc(short.id))
# factor levels arranged according to prevalence
genera.abundant$Rank6 <- factor(genera.abundant$Rank6, levels = top.8.most.prevalent)

# cumulative abundances by genera
genera.abundant.sort <- genera.abundant %>%
                        arrange(short.id, Rank6) %>%
                        group_by(short.id) %>%
                        mutate(cum.abundance = cumsum(Abundance)) %>%
                        ungroup() %>%
                        select(OTU:Rank1, Rank2, Rank6, cum.abundance)
cum.plot.dat <- genera.abundant.sort %>%
                group_by(Rank6) %>%
                summarise(as_tibble_row(quantile(cum.abundance, probs = seq(0.1, 0.9, by = 0.1)))) %>%
                pivot_longer(cols = "10%":"90%") %>%
                group_by(Rank6) %>%
                mutate(cum.difference = value-lag(value))

cum.plot.dat$cum.difference[is.na(cum.plot.dat$cum.difference)] <- cum.plot.dat$value[is.na(cum.plot.dat$cum.difference)]
cum.plot.dat$number.of.genera <- NA
cum.plot.dat$number.of.genera[which(cum.plot.dat$Rank6 == "Lactobacillus")] <- 1
cum.plot.dat$number.of.genera[which(cum.plot.dat$Rank6 == "Staphylococcus")] <- 2
cum.plot.dat$number.of.genera[which(cum.plot.dat$Rank6 == "Corynebacterium")] <- 3
cum.plot.dat$number.of.genera[which(cum.plot.dat$Rank6 == "Prevotella")] <- 4
cum.plot.dat$number.of.genera[which(cum.plot.dat$Rank6 == "Anaerococcus")] <- 5
cum.plot.dat$number.of.genera[which(cum.plot.dat$Rank6 == "Bacteroides")] <- 6
cum.plot.dat$number.of.genera[which(cum.plot.dat$Rank6 == "Stenotrophomonas")] <- 7
cum.plot.dat$number.of.genera[which(cum.plot.dat$Rank6 == "Streptococcus")] <- 8

cum.plot.dat.sort <- cum.plot.dat %>%
                      arrange(number.of.genera, name) %>%
                      rename(Percentile = name)

cum.plot.dat.sort$Percentile <- factor(cum.plot.dat.sort$Percentile, levels = rev(unique(cum.plot.dat.sort$Percentile)))

cols <- c("#fffdd0", "#ffe5b4","#e1c5a3", "#ffbf00", "#f24f26", "#9be3f9", "#0097c3",  "#003151","#222021")  

plot <- ggplot(cum.plot.dat.sort,aes(x=number.of.genera, y = cum.difference, fill = Percentile)) +
  geom_area() +
  scale_fill_manual(values=c(rev(cols)))+
  scale_x_continuous(name="Number of most prevalent genera", breaks=seq(1,8,1), expand = c(0, 0))+
  scale_y_continuous(name="Cumulative relative abundance", breaks=seq(0,1,0.1), expand = c(0, 0))+
  theme(axis.text.y   =  element_text(size=12),
        axis.text.x   = element_text(size=12),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title =  element_text(size=12),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill = "white"),
        plot.margin = unit(c(1,1,1,1), "cm"))
plot
ggsave("~/Figures/dominance.area.png", plot, width = 2000, height = 1500, units = "px")




###### Alpha diversity applying rarefaction values of 1002 sequences

# rarefaction
set.seed(150623)
phy.genera.rar <- rarefy_even_depth(phy.genera.no.zero.rel)

# all rarefactioned to 1002 reads

# alpha diversity
alpha.dat <- estimate_richness(phy.genera.rar, measures = c("Chao1", "Shannon", "Observed", "InvSimpson"))
alpha.dat$id <- rownames(alpha.dat)

# phylogenetic diversity
otu <- as.data.frame(phy.genera.rar@otu_table)
tree <- phy.genera.rar@phy_tree
# pd
df.pd <- pd(t(otu), tree, include.root=T)
df.pd$id <- rownames(df.pd)
# combine
alpha.div.df <- left_join(alpha.dat, df.pd, by = "id")
alpha.div.df <- alpha.div.df %>%
                select(id, Shannon, Observed, PD, Chao1)
row.names(alpha.div.df) <- alpha.div.df$id


# PCA on alpha diversity
alpha.pca <- prcomp(alpha.div.df[,c(2:5)], center = TRUE,scale. = TRUE)
summary(alpha.pca)

# scree plot
var_explained <- alpha.pca$sdev^2 / sum(alpha.pca$sdev^2)

scree.pca.alpha <- qplot(c(1:4), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Alpha Diversity Scree plot") +
  ylim(0, 1) +
  scale_x_continuous(breaks=seq(0,4,1)) +
  theme_base() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        #        strip.text.x = element_text(size = 10,face = "bold"), 
        panel.border=element_blank())
scree.pca.alpha


# loadings to PCs
alpha.loadings <- as.data.frame(alpha.pca$rotation[,c(1,2)])
alpha.loadings$score <- rownames(alpha.loadings)
alpha.loadings.long <- pivot_longer(alpha.loadings, cols = c("PC1", "PC2"), 
                                    names_to = "comp", values_to = "loading")

alpha.loadings.long$pos <- ifelse(alpha.loadings.long$loading < 0, FALSE, TRUE)
alpha.loadings.long$pos <- factor(alpha.loadings.long$pos, levels = c(TRUE, FALSE))

alpha.loadings.long$score <- as.factor(alpha.loadings.long$score)
levels(alpha.loadings.long$score) <- c("Chao1", "Observed\nSpecies", "Faith's\nPhylogenetic\nDiversity", "Shannon\nIndex" )

# plot Figure 3
p1 <- ggplot(alpha.loadings.long, aes(score, loading, fill = pos))+ 
  geom_bar(stat="identity", color="black") +
  scale_y_continuous(limits = c(-1,1)) +
  scale_fill_manual(values=c("#ffcc99", "#6699ff")) +
  theme_base() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 0.95), 
        strip.text.x = element_text(size = 12, face = "bold"), 
        panel.border=element_blank(),
        panel.spacing.x = unit(0,"line")) +
  facet_wrap(~comp, labeller = labeller(comp = 
                                          c("PC1" = "Alpha Diversity\nPrincipal Component 1 Loading",
                                            "PC2" = "Alpha Diversity\nPrincipal Component 2 Loading")
  ))  +
  xlab("") + 
  ylab("Correlation coefficient")
p1

ggsave("~/Figures/Alpha_loadings_rarefaction.png", p1, width = 3000, height =1500, units = "px")

# alpha PC data
alpha.pc.df <- data.frame(id = rownames(alpha.pca$x),
                          short.id = substr(rownames(alpha.pca$x), 4, nchar(rownames(alpha.pca$x))),
                          alpha.pc1 = alpha.pca$x[,1], 
                          alpha.pc2 = alpha.pca$x[,2])


######### Beta diversity applying rarefaction values of 1002 sequences

# midpoint rooted tree
tree <- phy.genera.rar@phy_tree
midpoint.tree <- midpoint.root(tree)

# otu transposed
otu.rar <- phy.genera.rar@otu_table
otu.rar.t <- t(otu.rar)

# weighted and unweighted generalised unifrac distances
set.seed(999)
unifracs <- GUniFrac(otu.rar.t, midpoint.tree, alpha = c(0, 1))$unifracs
dw <- unifracs[, , "d_1"] # Weighted UniFrac
du <- unifracs[, , "d_UW"] # Unweighted UniFrac

# principal coordinates of the distance matrix
pcoa.dw <- pcoa(dw, correction="none", rn=NULL)

val <- pcoa.dw$values
pcoa.dw$trace
sum(val$Relative_eig[1], val$Relative_eig[2])

barplot(val$Relative_eig[1:20], names = paste ('PCoA', 1:20), las = 3, ylab = 'eigenvalues')

pcoa.scree <- qplot(c(1:10), val$Relative_eig[1:10]) + 
  geom_line() + 
  xlab("Principal Coordinate") + 
  ylab("Variance Explained") +
  ggtitle("Beta Diversity Scree Plot") +
  ylim(0, 1) +
  scale_x_continuous(breaks=seq(0,10,1)) +
  theme_base() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        #        strip.text.x = element_text(size = 10,face = "bold"), 
        panel.border=element_blank())
pcoa.scree

# Supplementary Figure S1
both.screes <- grid.arrange(scree.pca.alpha, pcoa.scree, ncol = 2)

ggsave("~/Figures/Scree_plots.png", both.screes,width = 3000,height =1500, units = "px")



# beta PCo data
beta.pcoa.df <- data.frame(id = rownames(pcoa.dw$vectors),
                           short.id = substr(rownames(pcoa.dw$vectors), 4, nchar(rownames(pcoa.dw$vectors))),
                           beta.pco1 = pcoa.dw$vectors[,1], 
                           beta.pco2 = pcoa.dw$vectors[,2])


# correlation between rel abundance of genus and principal coordinate
cor.beta <- as.data.frame(cor(otu.rar.t, beta.pcoa.df$beta.pco1))
cor.beta$asv.id <- rownames(cor.beta)

temp <- as.data.frame(phy.genera.rar@tax_table)
temp$asv.id <- rownames(temp)
asv.key <- temp[,c("asv.id", "Rank6")]

cor.asv <- left_join(cor.beta, asv.key, by = "asv.id")
names(cor.asv)[1] <- "pco1"

cor.beta2 <- as.data.frame(cor(t(phy.genera.rar@otu_table), beta.pcoa.df$beta.pco2))
cor.beta2$asv.id <- rownames(cor.beta2)

cor.asv2 <- left_join(cor.beta2, asv.key, by = "asv.id")
names(cor.asv2)[1] <- "pco2"

cor <- left_join(cor.asv, cor.asv2, by = "asv.id")
cor.sel <- cor %>%
          select(Rank6.x, pco1, pco2)
cor.long <- pivot_longer(cor.sel, cols = c("pco1", "pco2"), values_to = "cor", names_to = "pco")
cor.long$pos <- ifelse(cor.long$cor < 0, FALSE, TRUE)
cor.long$pos <- factor(cor.long$pos, levels = c(TRUE, FALSE))
names(cor.long)[1] <- "genera"
cor.long.sort <- cor.long[order(cor.long$cor, decreasing = FALSE),]


# remove na
plot1 <- cor.long[cor.long$pco == "pco1" & cor.long$genera != "NA",]
plot2 <- cor.long[cor.long$pco == "pco2" & cor.long$genera != "NA",]

plot1.sort <- plot1[order(plot1$cor, decreasing = FALSE),]
plot2.sort <- plot2[order(plot2$cor, decreasing = FALSE),]

# full correlation results for supplementary
pco1.genera.all <- plot1.sort %>%
  mutate(cor.val = round(cor, 3)) %>%
  select(genera, cor.val) %>%
  arrange(desc(abs(cor.val)))
pco2.genera.all <- plot2.sort %>%
  mutate(cor.val = round(cor, 3)) %>%
  select(genera, cor.val) %>%
  arrange(desc(abs(cor.val)))

# save full correlation lists
write.csv(pco1.genera.all, "~/Data/cor.genera.beta.pco1.txt", sep="", row.names = F)
write.csv(pco2.genera.all, "~/Data/cor.genera.beta.pco2.txt", sep="", row.names = F)


# plot only loadings > 0.4, Figure 4
plot1.sort.04 <- plot1.sort[which(abs(plot1.sort$cor) > 0.4),]
plot2.sort.04 <- plot2.sort[which(abs(plot2.sort$cor) > 0.4),]


# pcoa1
p1 <- ggplot(plot1.sort.04, aes(x = reorder(genera, cor), y = cor, fill = pos))+ 
  geom_bar(stat="identity", color = "black") +
  scale_y_continuous(limits = c(-1,1)) +
  scale_fill_manual(values=c("#ffcc99", "#6699ff")) +
  theme_base() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 0.95, size = 12), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10,face = "bold"), 
        panel.border=element_blank(),
        panel.spacing.x = unit(0,"line")) +
  labs(title="          Generalized Unifrac PCo1 Genera Loadings >|0.4|", x ="", y= "Correlation coefficient")
p1

# pcoa2
p2 <- ggplot(plot2.sort.04, aes(x = reorder(genera, cor), y = cor, fill = pos))+ 
  geom_bar(stat="identity", color = "black") +
  scale_y_continuous(limits = c(-1,1)) +
  scale_fill_manual(values=c("#ffcc99", "#6699ff")) +
  theme_base() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 0.95, size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10,face = "bold"), 
        panel.border=element_blank(),
        panel.spacing.x = unit(0,"line")) +
  labs(title="          Generalized Unifrac PCo2 Genera Loadings >|0.4|",
       x ="", y= "Correlation coefficient")
p2

p3 <- grid.arrange(p1, p2, nrow = 2)
p3

ggsave("~/Figures/Beta_PCo_cor_04.png", p3,width = 3000,height =3200, units = "px")


# scatter plot of beta PCo1 and PCo2 by AB use, Supplementary Figure S3

ab <- quest.dat %>%
  select(id, birth_AB) %>%
  rename(short.id = id)
beta.pcoa.ab <- left_join(beta.pcoa.df, ab, by = "short.id")
beta.pcoa.ab$beta.pco1 <- as.numeric(beta.pcoa.ab$beta.pco1)
beta.pcoa.ab$beta.pco2 <- as.numeric(beta.pcoa.ab$beta.pco2)
beta.pcoa.ab$birth_AB <- factor(beta.pcoa.ab$birth_AB, levels = c(0,1), labels =c("non-AB", "AB"))


p <- ggplot(beta.pcoa.ab, aes(x = beta.pco2, y = beta.pco1, color = birth_AB)) +
  geom_point(size = 4, alpha=0.6) +
  scale_color_manual(values =  c("#0097c3", "#f24f26")) +
  xlab("Beta PCo2 (18%)") +
  ylab("Beta PCo1 (45%)") +
  theme(axis.text.y   =  element_blank(),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=16),
        axis.title.x  = element_text(size=16),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title = element_blank(),
        legend.position = c(0.12, 0.88),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "white"),
        plot.margin = unit(c(1,3,3,1), "cm"))
p
ggsave("~/Figures/Beta_PCo12.png", p,width = 2100,height =2000, units = "px")

