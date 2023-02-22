# FeOB Data Anlysis Using RStudio

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
library(grid)
library(labdsv)
library(indicspecies)
library(viridis)
library(ggplot2)
library('cowplot')
```

# P-value Formatter
```{r}
formatPvalues <- function(pvalue) {
  ra<- ""
  if(pvalue <= 0.1) ra<- "."
  if(pvalue <= 0.05) ra<- "*"
  if(pvalue <= 0.01) ra<- "**"
  if(pvalue <= 0.001) ra<- "***"
  return(ra)
}
```

# Dada2 Pipeline
Fastq file analysis

1) Identify path of files
```{r}
path <- "/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/fastq_files"
list.files(path)
```

2) Forward & Reverse Strands
We read in the names of the fastq files, and perform some string manipulation to get lists of the forward and reverse fastq files in matched order
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
Now we visualize the quality profile of the reverse reads. 

In gray-scale is a heat map of the frequency of each quality score at each base position. The median quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. The red line shows the scaled proportion of reads that extend to at least that position (this is more useful for other sequencing technologies, as Illumina reads are typically all the same lenghth, hence the flat red line).

The forward reads are good quality. We generally advise trimming the last few nucleotides to avoid less well-controlled errors that can arise there. These quality profiles do not suggest that any additional trimming is needed. We will truncate the forward reads at position 240 (trimming the last 10 nucleotides).
```{r}
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

5) Filter & trim
The standard filtering parameters are starting points, not set in stone. 

If you want to speed up downstream computation, consider tightening maxEE. If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads and reducing the truncLen to remove low quality tails. 

Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later.

The maxEE=c(8,8) setting is saying that there can be a max of 8 ambiguous nucleotides in a row for each forward and reverse read before the read is tossed out.

```{r}
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT"
trimLeft = c(FWD, REV)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = c(19,20))
head(out)
```

6) Learning & Plotting Error rates
Points are the observed error rates for each consensus quality score. The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score
```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```
104655519 total bases in 400979 reads from 24 samples will be used for learning the error rates.
101224800 total bases in 562360 reads from 35 samples will be used for learning the error rates.

Visualize Errors
```{r}
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
dadaFs[[1]]
```
dada-class: object describing DADA2 denoising results
319 sequence variants were inferred from 4141 
input unique sequences.
Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

```{r}
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaRs[[1]]
```
dada-class: object describing DADA2 denoising results
245 sequence variants were inferred from 4096 
input unique sequences.
Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

9) Merging paired reads
The mergers object is a list of data.frames from each sample. Each data.frame contains the merged sequence, its abundance, and the indices of the forward and reverse sequence variants that were merged. Paired reads that did not exactly overlap were removed by mergePairs, further reducing spurious output.
```{r}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE, verbose=TRUE)
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
The core dada method corrects substitution and indel errors, but chimeras remain. Fortunately, the accuracy of the sequence variants after denoising makes identifying chimeras simpler than it is when dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.
```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```
Identified 31390 bimeras out of 46629 input sequences.
[1]    89 15239
[1] 0.9192429

12) Tracking reads through the pipeline
```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

13) Assigning taxonomy
It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to assign taxonomy to the sequence variants. 
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```
# Formatting files for further analysis
These files will be easier to work with for different statistical testing, plotting and other uses further down in analysis
d
## Seq Table
```{r}
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
```

## Fasta Header as ASV_
```{r}
for (i in 1:dim(seqtab.nochim)[2]) {
asv_headers[i] <- paste(">ASV", i, sep="_")
}
```

## ASV Sequences
```{r}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
```

## ASV Abundance
```{r}
asv_otu <- t(seqtab.nochim)
row.names(asv_otu) <- sub(">", "", asv_headers)
```

## ASV Taxonomy
```{r}
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
```

## Merging Abundance and Tax Table
```{r}
otu_tax_table <- merge(asv_otu, asv_tax, by=0)
```

## Write output files:
```{r}
write(asv_fasta, "asv_fasta_adjusted.fa")
write.table(asv_otu, "asv_otu_adjusted.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "asv_tax_adjusted.csv", sep=",", quote=F, col.names=NA)
write.table(otu_tax_table, "otu_tax_table_adjusted.csv", sep=",", quote=F, col.names=NA)
```

# Formating files to join data.frames
```{r}
metadata_bio_water <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/metadata_bio_water.csv")
metadata_bio_water
```

# Generate a OTU Count File:
```{r}
otu_counts <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/asv_otu.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts
```

```{r}
taxonomy <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/asv_tax.csv")
taxonomy
```

# Generate an OTU Relative Abundance Table:
This will be used for NMDS plots & statistical testing
```{r}
otu_rel_abund_bio_water <- inner_join(metadata_bio_water, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Genus", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_bio_water
```

```{r}
write.table(otu_rel_abund_bio_water, "otu_rel_abund_adjusted_biofilm_water.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu_rel_abund_bio_water <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/otu_rel_abund_adjusted_biofilm_water.csv")
```

# Determining Taxonomic Abundances & Community Composition
## Stacked Barcharts
Using stacked barcharts will allow us to visualize differences in microbial communities for each of the sample locations. A great first step to compare samples!

### Stacked Barcharts: All Samples
```{r}
otu_rel_abund_bio_water %>%
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
ggsave("phylum_stacked_bar_adjusted_water_biofilm.tiff", width=20, height=7)
```

### Stacked Barcharts: Biofilms Samples Only
```{r}
metadata_bio <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/metadata_biofilm.csv")
metadata_bio
```

```{r}
otu_rel_abund_bio <- inner_join(metadata_bio, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Genus", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_bio
```

```{r}
otu_rel_abund_bio %>%
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
```

```{r}
ggsave("phylum_stacked_bar_adjusted_biofilm.tiff", width=15, height=7)
```

### Stacked Barcharts: Water Samples Only
```{r}
metadata_water <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/metadata_water.csv")
metadata_water
```

```{r}
otu_rel_abund_water <- inner_join(metadata_water, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Genus", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_water
```

```{r}
otu_rel_abund_water %>%
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
```

```{r}
ggsave("phylum_stacked_bar_adjusted_water.tiff", width=12, height=7)
```

### Stacked Barcharts: Starboard vs Port Side
```{r}
metadata_star_port <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/metadata_star_port_only.csv")
metadata_star_port
```

```{r}
otu_rel_abund_SP <- inner_join(metadata_star_port, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Genus", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_SP
```

```{r}
write.table(otu_rel_abund_SP, "otu_rel_abund_adjusted_SP.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu_rel_abund_SP %>%
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
```

```{r}
ggsave("phylum_stacked_bar_adjusted_SP.tiff", width=10, height=8)
```

## NMDS Plots
A way to condense information from multidimensional data (multiple variables/species/ASVs), into a 2D representation or ordination. The closer 2 points are, the more similar the corresponding samples are with respect to the variables that went into making the NMDS plot.

### NMDS: Location
```{r}
pc = read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/nmds_asv_otu_biofilm_all_water_together.csv")
pc
```

```{r}
com = pc[,6:ncol(pc)]
com
m_com <- as.matrix(com)
m_com
```

```{r}
set.seed(1000)
nmds = metaMDS(m_com, distance = "bray")
```

```{r}
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$location = pc$location
data.scores$sample_id = pc$sample_id
head(data.scores)
```

Make NMDS Plot
```{r}
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Sample Location")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx
ggsave("Shipwreck_NMDS_location_adjusted_bio_water.tiff")
```

### Water All Combined vs Biofilm
```{r}
pcBW <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/nmds_asv_otu_biofilm_all_water_together.csv")

comBW = pcBW[,6:ncol(pcBW)]
m_comBW <- as.matrix(comBW)

set.seed(1000)
nmdsBW = metaMDS(m_comBW, distance = "bray")

data.scores.BW <- as.data.frame(scores(nmdsBW)$sites)
data.scores.BW$location = pcBW$location
data.scores.BW$sample_id = pcBW$sample_id
head(data.scores.BW)
```

```{r}
xxBW = ggplot(data.scores.BW, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Sample Location")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xxBW

ggsave("Shipwreck_NMDS_location_adjusted_biofilm_water.tiff")
```

### NMDS Depth
```{r}
data.scores.BW2 <- as.data.frame(scores(nmdsBW)$sites)
data.scores.BW2$depth = pcBW$depth
data.scores.BW2$sample_id = pcBW$sample_id
head(data.scores.BW2)
```

```{r}
data.scores.BW2 <- as.data.frame(scores(nmdsBW)$sites)
data.scores.BW2$depth = pcBW$depth
data.scores.BW2$sample_id = pcBW$sample_id
head(data.scores.BW2)
write.table(data.scores.BW2, "nmds_data.scores_all_depth_bio_water.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xxBW2 = ggplot(data.scores.BW2, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(depth = depth, colour=factor(depth))) +
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Sample Depth") +
 theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Depth (in ft)", y = "NMDS2")
xxBW2
ggsave("Shipwreck_NMDS_depth_adjusted.tiff")
```

### NMDS Starboard vs Port Side
```{r}
pc2 <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/nmds_asv_otu_star_port_biofilm_water.csv")
pc2
```

```{r}
com2 = pc2[,7:ncol(pc2)]
com2

m_com2 <- as.matrix(com2)
m_com2
```

```{r}
set.seed(1000)
nmds2 = metaMDS(m_com2, distance = "bray")
```

```{r}
data.scores2 <- as.data.frame(scores(nmds2)$sites)
data.scores2$location = pc2$location
data.scores2$sample_id = pc2$sample_id
head(data.scores2)
```

```{r}
xx4 = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(shape = location, colour = location))+
  ggtitle("NMDS Ordination - Sample Location") +
  scale_fill_discrete()+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx4
```

```{r}
ggsave("Shipwreck_NMDS_Star_Port_adjusted_bio_water.tiff")
```

### NMDS Waterlevel
```{r}
data.scores.BW3 <- as.data.frame(scores(nmdsBW)$sites)
data.scores.BW3$water_level = pcBW$water_level
data.scores.BW3$sample_id = pcBW$sample_id
head(data.scores.BW3)
```

```{r}
xxSPBW = ggplot(data.scores.BW3, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(shape = water_level, colour = water_level))+
  ggtitle("NMDS Ordination - Sample Water-level") +
  scale_fill_discrete()+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "water_level", y = "NMDS2")
xxSPBW
```

```{r}
ggsave("Shipwreck_NMDS_water_level_adjusted.tiff")
```

# NMDS: Biofilm Samples Only

```{r}
pcB = read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/nmds_asv_otu_biofilm.csv")
pcB
```

```{r}
comB = pcB[,6:ncol(pcB)]
comB
m_comB <- as.matrix(comB)
m_comB

set.seed(1000)
nmdsB = metaMDS(m_comB, distance = "bray")
```

```{r}
data.scoresB <- as.data.frame(scores(nmdsB)$sites)
data.scoresB$location = pcB$location
data.scoresB$sample_id = pcB$sample_id
head(data.scoresB)
```

Make NMDS Plot
```{r}
xxB = ggplot(data.scoresB, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Sample Location")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xxB
```

```{r}
ggsave("Shipwreck_NMDS_location_adjusted_biofilm_only.tiff")
```

### Biofilm Only: Star vs Port
```{r}
pcSP <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/nmds_asv_otu_star_port_biofilm.csv")
pcSP
```

```{r}
comSP = pcSP[,7:ncol(pcSP)]
comSP

m_comSP <- as.matrix(comSP)
m_comSP
```

```{r}
set.seed(1000)
nmdsSP = metaMDS(m_comSP, distance = "bray")
```

```{r}
data.scoresSP <- as.data.frame(scores(nmdsSP)$sites)
data.scoresSP$location = pcSP$location
data.scoresSP$sample_id = pcSP$sample_id
head(data.scoresSP)
```

```{r}
xx4 = ggplot(data.scoresSP, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(shape = location, colour = location))+
  ggtitle("NMDS Ordination - Sample Location") +
  scale_fill_discrete()+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx4
```

```{r}
ggsave("Shipwreck_NMDS_Star_Port_adjusted_star_port_biofilm.tiff")
```

# Phyloseq
```{r}
theme_set(theme_bw())
```

```{r}
samples.out <- rownames(seqtab.nochim)
samples <- sapply(strsplit(samples.out, "Shos-"), `[`, 2)
samples
```

```{r}
map <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/metadata_biofilm.csv")
map
```

```{r}
location <-map[,2,ncol(map)]
location
```

```{r}
depth <- map[,3,ncol(map)]
depth
```

```{r}
waterlevel <- map[,4,ncol(map)]
waterlevel
```

```{r}
samdf <- data.frame(Location=location, Depth=depth, Waterlevel=waterlevel)
samdf$Waterlevel <- "Waterline"
samdf$Waterlevel[samdf$Depth>0] <- "Below Water"
rownames(samdf) <- samples.out
```

```{r}
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))
ps <- prune_samples(sample_names(ps), ps)
```

```{r}
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 15239 taxa and 89 samples ]
sample_data() Sample Data:       [ 89 samples by 3 sample variables ]
tax_table()   Taxonomy Table:    [ 15239 taxa by 6 taxonomic ranks ]
refseq()      DNAStringSet:      [ 15239 reference sequences ]

```{r}
plot_richness(ps, x="Depth", measures=c("Shannon", "Simpson"), color="Waterlevel")
```

```{r}
plot_richness(ps, x="Depth", measures=c("Shannon", "Simpson"), color="Location")
```

```{r}
plot_richness(ps, x="Location", measures=c("Shannon", "Simpson"), color="Depth")
```

```{r}
plot_richness(ps, x="Location", measures=c("Shannon", "Simpson"), color="Waterlevel")
```

```{r}
plot_richness(ps, x="Waterlevel", measures=c("Shannon", "Simpson"), color="Location")
```

```{r}
plot_richness(ps, x="Waterlevel", measures=c("Shannon", "Simpson"), color="Depth")
```

```{r}
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Waterlevel", title="Bray NMDS")
```

```{r}
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Location", title="Bray NMDS")
```

```{r}
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Depth", title="Bray NMDS")
```

# Statistical Analysis: ANOSIM, ANOVA, SIMPER
# ANOSIM Statistical Testing
## Location Specific Biofilm & Water
Make community matrix - extract columns with abundance information, turn data frame into matrix
```{r}
pcL <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/Anosim_Sed_S_P_RDP_BK_AFST_W.csv")
pcL

comL = pcL[,3:ncol(pcL)]
comL

m_comL = as.matrix(comL)
m_comL
```

```{r}
anoL = anosim(m_comL, pcL$location, distance = "bray", permutations = 9999)
anoL
```

## Depth Biofilm
```{r}
pcD <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/ANOSIM_Depth_Biofilm.csv")
pcD

comD = pcD[,3:ncol(pcD)]
comD

m_comD = as.matrix(comD)
m_comD
```

```{r}
anoD = anosim(m_comD, pcD$depth, distance = "bray", permutations = 9999)
anoD
```

## Water-level Biofilm
```{r}
pcWL <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/Anosim_Waterlevel_Biofilm.csv")
pcWL

comWL = pcWL[,3:ncol(pcWL)]
comWL

m_comWL = as.matrix(comWL)
m_comWL
```

```{r}
anoWL = anosim(m_comWL, pcWL$water_level, distance = "bray", permutations = 9999)
anoWL
```

## Testing significance between microenvironments
### Shipwreck_Biofilm, Sediment_Biofilm and Water_Samples
```{r}
pcSSW <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/Anosim_Ship_Sed_Water.csv")
pcSSW

comSSW = pcSSW[,3:ncol(pcSSW)]
comSSW

m_comSSW = as.matrix(comSSW)
m_comSSW
```

```{r}
anoSSW <- anosim(m_comSSW, pcSSW$location, distance = "bray", permutations = 9999)
anoSSW
```

### Shipwreck_biofilm (Starboard & Port), Sediment_biofilm, Bulkhead_biofilm, Rudder Post_biofilm, Aft_Starboard_biofilm, Water_sample
```{r}
pcSSBRAW <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/Anosim_Ship_Sed_Water_Bulk_Aft_Rud.csv")
pcSSBRAW

comSSBRAW = pcSSBRAW[,3:ncol(pcSSBRAW)]
comSSBRAW

m_comSSBRAW = as.matrix(comSSBRAW)
m_comSSBRAW
```

```{r}
anoSSBRAW <- anosim(m_comSSBRAW, pcSSBRAW$location, distance = "bray", permutations = 9999)
anoSSBRAW
```

## Testing significance between Starboard side & Port side
### SP Location Biofilm
```{r}
pcSPL <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/ANOSIM_Location_SP_Biofilm.csv")
pcSPL

comSPL = pcSPL[,3:ncol(pcSPL)]
comSPL

m_comSPL = as.matrix(comSPL)
m_comSPL
```

```{r}
anoSPL = anosim(m_comSPL, pcSPL$location, distance = "bray", permutations = 9999)
anoSPL
```

### SP Depth Biofilm
```{r}
pcSPD <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/ANOSIM_Depth_SP_Biofilm.csv")
pcSPD

comSPD = pcSPD[,3:ncol(pcSPD)]
comSPD

m_comSPD = as.matrix(comSPD)
m_comSPD
```

```{r}
anoSPD = anosim(m_comSPD, pcSPD$depth, distance = "bray", permutations = 9999)
anoSPD
```

### SP Water-level Biofilm
```{r}
pcSPWL <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/ANOSIM_Waterlevel_SP_Biofilm.csv")
pcSPWL

comSPWL = pcSPWL[,3:ncol(pcSPWL)]
comSPWL

m_comSPWL = as.matrix(comSPWL)
m_comSPWL
```

```{r}
anoSPWL = anosim(m_comSPWL, pcSPWL$water_level, distance = "bray", permutations = 9999)
anoSPWL
```

When interpreting these results you want to look at the ANOSIM statistic R and the Significance values. 

“The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups to the mean of ranked dissimilarities within groups. 

**An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high and low ranks within and between groups”** (GUSTAME). In other words, the **higher the R value, the more dissimilar your groups are in terms of microbial community composition!**

**A Significance value less than 0.05 is generally considered to be statistically significant**, and means the null hypothesis can be rejected. Therefore, there is a **statistically significant difference in the microbial communities between your groups.** **Greater than 0.05, means that there is no statistical difference between the microbial communities in your groups.**

# Diversity Index Value Generating
## Biofilm & Water Samples
```{r}
otu_table <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/asv_otu.csv", header=T, row.names=1, check.names=FALSE)
```

Transpose the data to have sample names on rows
```{r}
otu.table.diver <- t(otu_table)
otu.table.diver <- as.data.frame(otu.table.diver)
head(otu.table.diver)
```

### Shannon
```{r}
data(otu.table.diver)
H<-diversity(otu.table.diver)
H
```

### Richness
```{r}
richness <- specnumber(otu.table.diver)
```

### Pielou Evenness
```{r}
eveness <- H/log(richness)
```

```{r}
metadata_bio_water <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/metadata_bio_water.csv")
metadata_bio_water
```

```{r}
alphaBW <- cbind(shannon = H, richness = richness, pielou = eveness, metadata_bio_water)
write.csv(alphaBW, "diversity_indices_bio_water.csv")
head(alphaBW)
```

```{r}
plot.shan <- ggplot(alphaBW, aes(x = location, y = shannon, colour = location)) +
  geom_point(size = 3) +
  ylab("Shannon's H'") + 
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan
ggsave("Shannon_Location_bio_water.tiff")
```

```{r}
plot.rich <-ggplot(alphaBW, aes(x = location, y = richness, colour = location)) +
  geom_point(size = 3) +
  ylab("Species Richness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich
ggsave("Richness_Location_bio_water.tiff")
```

```{r}
plot.even <- ggplot(alphaBW, aes(x = location, y = pielou, colour = location)) +
  geom_point(size = 3) +
  ylab("Pielou's Evenness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even
ggsave("Pielou's_Evenness_Location_bio_water.tiff")
```

```{r}
legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("Shannon_Richness_Eveness_bio_water.tiff")
```

## Biofilm Only
```{r}
otu_table_bio <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/asv_otu_biofilm.csv", header=T, row.names=1, check.names=FALSE)
head(otu_table_bio)
```

Transpose the data to have sample names on rows
```{r}
otu.table.diver.bio <- t(otu_table_bio)
otu.table.diver.bio <- as.data.frame(otu.table.diver.bio)
head(otu.table.diver.bio)
#write.csv(otu.table.diver.bio, "otu.table.diversity.biofilm.csv")
```

```{r}
data(otu.table.diver.bio)
H<-diversity(otu.table.diver.bio)
H
```

```{r}
richness <- specnumber(otu.table.diver.bio)

eveness <- H/log(richness)
```

```{r}
metadata_bio <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/metadata_biofilm.csv")
metadata_bio
```

```{r}
alphaBio <- cbind(shannon = H, richness = richness, pielou = eveness, metadata_bio)
write.csv(alphaBio, "diversity_indices_biofilm.csv")
head(alphaBio)
```

```{r}
plot.shan <- ggplot(alphaBio, aes(x = location, y = shannon, colour = location)) +
  geom_point(size = 3) +
  ylab("Shannon's H'") + 
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan
ggsave("Shannon_Location_biofilm.tiff")
```

```{r}
plot.rich <-ggplot(alphaBio, aes(x = location, y = richness, colour = location)) +
  geom_point(size = 3) +
  ylab("Species Richness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich
ggsave("Richness_Location_biofilm.tiff")
```

```{r}
plot.even <- ggplot(alphaBio, aes(x = location, y = pielou, colour = location)) +
  geom_point(size = 3) +
  ylab("Pielou's Evenness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even
ggsave("Pielou's_Evenness_Location_biofilm.tiff")
```

```{r}
legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("Shannon_Richness_Eveness_biofilm.tiff")
```

# ANOVA
## Location
```{r}
anova <- aov(shannon ~ location, alpha)
summary(anova)
```

## Depth
```{r}
anova <- aov(shannon ~ depth, alpha)
summary(anova)
```

## Water_Level
```{r}
anova <- aov(shannon ~ water_level, alpha)
summary(anova)
```  

# SIMPER: Pairwise Dissimilarity
```{r}
otu.table.diver.mdf.bio <- as.matrix.data.frame(otu.table.diver.bio)
rownames(otu.table.diver.mdf.bio) <- metadata_bio$location

otu.table.diver.bray.bio <- vegdist(otu.table.diver.mdf.bio, method="bray")
otu.table.diver.bray.bio
```

## Location
```{r}
simper2 <- simper(otu.table.diver.bio, metadata_bio$location, permutations=999)
options(max.print=999999)
summary(simper2)
dput(simper2, file = "simp_location.txt")
sim2 <- dget("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/doc/simp_location.txt")
summary(sim2)
```

## Depth
```{r}
simper <- simper(otu.table.diver.bio, metadata_bio$depth, permutations=999)
summary(simper)
options(max.print=999999)
```

## Water Level
```{r}
simper3 <- simper(otu.table.diver.bio, metadata_bio$water_level, permutations=999)
summary(simper3)
options(max.print=999999)

# PCoA
```{r}
#calculate principle coordinate analysis (Bray-Curtis)
pcoa.otu.table.diver.bray.bio <- cmdscale(otu.table.diver.bray.bio, k=2, eig=T)
```

```{r}
#extract axis positions for each location from cmdscale object & create dataframe for plotting
pcoa.otu.table.diver.bray.plotting.bio <- as.data.frame(pcoa.otu.table.diver.bray.bio$points)
colnames(pcoa.otu.table.diver.bray.plotting.bio) <- c("axis_1", "axis_2")
pcoa.otu.table.diver.bray.plotting.bio$location <- rownames(pcoa.otu.table.diver.bray.plotting.bio)
```

```{r}
#calculate proportion of variance in the data which is explained by the first 2 PCOA axes
pcoa.otu.table.diver.bray.bio$eig[1]/(sum(pcoa.otu.table.diver.bray.bio$eig))
```

```{r}
pcoa.otu.table.diver.bray.bio$eig[2]/(sum(pcoa.otu.table.diver.bray.bio$eig))
```

```{r}
#create PCoA plot
pcoa.otu.bray.plot.bio <- ggplot(pcoa.otu.table.diver.bray.plotting.bio, aes(x = axis_1, y = axis_2, colour = location)) +
  geom_point(size = 3) +
  theme_bw() + 
  xlab("PCoA 1 (11.1%)") +
  ylab("PCoA 2 (6.14%)") +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1)
pcoa.otu.bray.plot.bio
ggsave("PCoA_Location_Biofilm.tiff", width=15, height=7) 
```

# PERMANOVA
Homogeneity of dispersion test
```{r}
permutest(betadisper(otu.table.diver.bray.bio, metadata_bio$location))
```

```{r}
permutest(betadisper(otu.table.diver.bray.bio, metadata_bio$depth))
```

```{r}
permutest(betadisper(otu.table.diver.bray.bio, metadata_bio$water_level))
```

# Adonis Analysis
```{r}
adonis(otu.table.diver.bray.bio ~ location, data = metadata_bio, permutations = 999)
```

```{r}
adonis(otu.table.diver.bray.bio ~ depth, data = metadata_bio, permutations = 999)
```

```{r}
adonis(otu.table.diver.bray.bio ~ water_level, data = metadata_bio, permutations = 999)
```
