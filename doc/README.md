
# FeOB Data Analysis Using RStudio

Load all necessary packages:
```{r}
library(dada2)
library(phyloseq)
library(ggplot2)
library(ape)
library(microbiome)
library(tidyverse)
library(kableExtra)
library(phangorn)
library(DECIPHER)
library(reshape2)
library(treeio)
library(ggtree)
library(ggstance)
library(scales)
library(dplyr)
library(ggpattern)
library(vegan)
library(MASS)
library(ecodist)
library(scatterplot3d)
library(fso)
library(dplyr)
library(svglite)
library(ggstatsplot)
```

# dada2 fastq file analysis
1) Identify path of files
```{r}
path <- "/Users/maggieshostak/Desktop/Shipwreck/data/fastq_files"
list.files(path)
```

2) Forward & Reverse Strands
```{r}
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
```

3) Extract sample names
```{r}
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)
```

4) Inspect read quality scores
```{r}
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

5) Filter & trim
```{r}
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)
```

6) Learning & Plotting Error rates
```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
```

7) Dereplication
This combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence. It substantially reduces computation time by eliminating redundant comparisons.
```{r}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

8) Sample inference
```{r}
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]
```

9) Merging paired reads
The mergers object is a list of data.frames from each sample. Each data.frame contains the merged sequence, its abundance, and the indices of the forward and reverse sequence variants that were merged. Paired reads that did not exactly overlap were removed by mergePairs, further reducing spurious output.
```{r}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])
```

10) Sequence table construction
We can now construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods. The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. 
```{r}
seqtab <- makeSequenceTable(mergers)

dim(seqtab)

table(nchar(getSequences(seqtab)))
```

11) Removing chimeras
```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
```

12) Tracking reads through the pipeline
```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

13) Assigning taxonomy
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/Shipwreck/data/fastq_files/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

# Formatting files for further analysis

```{r}
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))

asv_otu <- t(seqtab.nochim)
row.names(asv_otu) <- sub(">", "", asv_headers)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)

otu_tax_table <- merge(asv_otu, asv_tax, by=0)
```

Write output files:
```{r}
write(asv_fasta, "asv_fasta.fa")
write.table(asv_otu, "asv_otu.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "asv_tax.csv", sep=",", quote=F, col.names=NA)
write.table(otu_tax_table, "otu_tax_table.csv", sep=",", quote=F, col.names=NA)
```

# Formating files to join data.frames

```{r}
metadata
```

```{r}
metadata %>% count(location)
```

```{r}
metadata %>% count(depth)
```

```{r}
metadata %>% count(water_level)
```

# Formating files to join data.frames
Split *_asv_otu file_* into 2 separate dataframe to manually transpose in excel, then merged: 
```{r}
df1 <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/df_1.csv")
df1

df2 <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/df_2.csv")
df2

df_list <- list(df1, df2)
otu_count <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list) %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count
```

Check column headers for joining data.frames:
```{r}
taxonomy
metadata
otu_count
```

# Generate an OTU Relative Abundance Table:
This will be used for NMDS plots & statistical testing

```{r}
otu_rel_abund <- inner_join(metadata, otu_count, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Genus", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund
write.table(otu_rel_abund, "otu_rel_abund.csv", sep=",", quote=F, col.names=NA)
```

# Stacked Barcharts: All Samples
Using stacked barcharts will allow us to visualize differences in microbial communities for each of the sample locations. A great first step to compare samples!
```{r}
otu_rel_abund %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("phylum_stacked_bar.tiff", width=15, height=7)
```

## Starboard vs Port Side Samples

### Starboard:
```{r}
otu_rel_abund %>%
  filter(level=="Phylum", location == "Starboard") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("Starboard_stacked_bar.tiff", width=5, height=7)
```

### Port:
```{r}
otu_rel_abund %>%
  filter(level=="Phylum", location == "Port") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("Port_stacked_bar.tiff", width=5, height=7)
```

### Sediment Sample Comparison
```{r}
otu_rel_abund_sed <- otu_rel_abund %>% 
  filter(location == "Sediment")

write.csv(otu_rel_abund_sed, "otu_rel_abund_sed.csv")
```

```{r}
otu_rel_abund_sed %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("Sediment_stacked_bar.tiff", width=10, height=7)
```

# NMDS Plots: Location
A way to condense information from multidimensional data (multiple variables/species/ASVs), into a 2D representation or ordination. The closer 2 points are, the more similar the corresponding samples are with respect to the variables that went into making the NMDS plot.

```{r}
pc = read.csv("/Users/maggieshostak/Desktop/Shipwreck/data/metadata_asv.csv")
pc
```

```{r}
com = pc[,5:ncol(pc)]
com
m_com <- as.matrix(com)
m_com
```

```{r}
set.seed(100)
nmds = metaMDS(m_com, distance = "bray")
```

```{r}
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$location = pc$location
data.scores$sample_id = pc$sample_id
head(data.scores)
```

Make NMDS Plot:
```{r}
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 5, aes(shape = location, colour = location))+
 theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2") +
 scale_colour_manual(values = c("#7CB342", "#039BE5", "#EC407A", "Black", "#E53935", "#212121", "#FF9900", "#FFCC99", "Orange"))
xx
```

```{r}
ggsave("Shipwreck_NMDS.tiff")
```

# NMDS Depth
```{r}
data.scores2 <- as.data.frame(scores(nmds)$sites)
data.scores2$depth = pc$depth
data.scores2$sample_id = pc$sample_id
head(data.scores2)
```

```{r}
xx = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 5, aes(depth = depth, shape=factor(depth), colour=factor(depth))) +
  scale_fill_discrete()
 theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Depth", y = "NMDS2")
xx
```

```{r}
ggsave("Shipwreck_NMDS_depth.tiff")
```

# NMDS Waterlevel
```{r}
data.scores3 <- as.data.frame(scores(nmds)$sites)
data.scores3$water_level = pc$water_level
data.scores3$sample_id = pc$sample_id
head(data.scores3)
```

```{r}
xx = ggplot(data.scores3, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 5, aes(shape = water_level, colour = water_level))+
 theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "water_level", y = "NMDS2")
xx
```

```{r}
ggsave("Shipwreck_NMDS_water_level.tiff")
```

# Diversity Index Value Generating
```{r}
otu_table <- read.csv("/Users/maggieshostak/Desktop/Shipwreck/data/asv_otu_no_chloroplast.csv", header=T, row.names=1, check.names=FALSE)
```

Transpose the data to have sample names on rows
```{r}
t_otu_table <-t(otu_table)
data(t_otu_table)
H<-diversity(t_otu_table)
simp<-diversity(t_otu_table, "simpson")
invsimp<-diversity(t_otu_table, "inv")
```

Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
```{r}
unbias.simp <- rarefy(t_otu_table, 2) - 1
```

Fisher alpha
```{r}
alpha <- fisher.alpha(t_otu_table)
```

Species richness (S) and Pielou's evenness (J):
```{r}
S <- specnumber(t_otu_table)
J <- H/log(S)
```

Plot all
```{r}
pairs(cbind(H, simp, invsimp, unbias.simp, alpha), pch="+", col="red")
alpha
```

```{r}
write.table(J, "/Users/maggieshostak/desktop/pielou_evenness.txt", sep="\t")
```

# ANOSIM Test

```{r}
pc = read.csv("/Users/maggieshostak/Desktop/Shipwreck/data/metadata_asv.csv", header= TRUE)
```

```{r}
#make community matrix - extract columns with abundance information, turn data frame into matrix
com = pc[,5:ncol(pc)]
m_com = as.matrix(com)
```

```{r}
ano = anosim(m_com, pc$location, distance = "bray", permutations = 9999)
ano
```

Results: 
