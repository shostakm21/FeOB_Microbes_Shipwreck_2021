# FeOB Data Anlysis Using RStudio

## Load all necessary packages:
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
library(skimr)
library(kableExtra)
library(data.table)
library(taxa)
```

## dada2 fastq file analysis
1) Identify path of files
```{r}
path <- "/Users/maggieshostak/Desktop/Shipwreck/data/fastq_files"
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

The forward reads are good quality. We generally advise trimming the last few nucleotides to avoid less well-controlled errors that can arise there. These quality profiles do not suggest that any additional trimming is needed. 

```{r}
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

5) Filter & trim
The standard filtering parameters are starting points, not set in stone. 

If you want to speed up downstream computation, consider tightening maxEE. If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads and reducing the truncLen to remove low quality tails. 

Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later.

The maxEE=c(8,8) setting is saying that there can be a max of 8 ambiguous nucleotides in a row for each forward and reverse read before the read is tossed out.

FWD & REV are the sequences for the universal primers (V4-V5) listed on IMRs website.

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
```{r}
dada-class: object describing DADA2 denoising results
319 sequence variants were inferred from 4141 
input unique sequences.
Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 
16
```

```{r}
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaRs[[1]]
```

```{r}
dada-class: object describing DADA2 denoising results
245 sequence variants were inferred from 4096 input unique sequences.
Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
```

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
```{r}
[1]    89 46629

  451 
46629 
```

11) Removing chimeras
The core dada method corrects substitution and indel errors, but chimeras remain. Fortunately, the accuracy of the sequence variants after denoising makes identifying chimeras simpler than it is when dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.

```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```
```{r}
Identified 31390 bimeras out of 46629 input sequences.
  89 15239
  0.9192429
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
It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to assign taxonomy to the sequence variants. 
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/Shipwreck/data/fastq_files/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

# Formatting Files for Further Analysis
```{r}
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
```

```{r}
for (i in 1:dim(seqtab.nochim)[2]) {
asv_headers[i] <- paste(">ASV", i, sep="_")
}
```

```{r}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
```

```{r}
asv_otu <- t(seqtab.nochim)
row.names(asv_otu) <- sub(">","", asv_headers)
```

```{r}
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
```

```{r}
otu_tax_table <- merge(asv_otu, asv_tax, by=0)
```

## Output Files:
```{r}
write(asv_fasta, "asv_fasta_adjusted.fa")
write.table(asv_otu, "asv_otu_adjusted.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "asv_tax_adjusted.csv", sep=",", quote=F, col.names=NA)
write.table(otu_tax_table, "OTU_TAX_table_adjusted.csv", sep=",", quote=F, col.names=NA)
```

## Checking counts for different variables:
```{r}
metadata %>% count(location)
```

```{r}
metadata %>% count(depth)
```

```{r}
metadata %>% count(water_level)
```

## Formating files to join data.frames

```{r}
df1 <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/df_1.csv")
df1
df2 <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/df_2.csv")
df2
df_list <- list(df1, df2)
df_list

otu_count <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list) %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count

write.table(otu_count, "otu_count.csv", sep=",", quote=F, col.names=NA)
```

```{r}
taxonomy
metadata
otu_count
```


## Generate an OTU Relative Abundance Table:
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

## Stacked Barcharts: All Samples
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

### Sediment Sample Comparison
```{r}
otu_rel_abund_sed <- otu_rel_abund %>% 
  filter(location == "Sediment")

write.csv(otu_rel_abund_sed, "otu_rel_abund_sed.csv")

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

## NMDS Plots
A way to condense information from multidimensional data (multiple variables/species/ASVs), into a 2D representation or ordination. The closer 2 points are, the more similar the corresponding samples are with respect to the variables that went into making the NMDS plot.

```{r}
pc <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/nmds_asv_otu.csv")
pc
```

You want to make a data frame that only includes your species/OTU abundance columns
```{r}
com = pc[,7:ncol(pc)]
com
```

Often in R you will get errors because your data is not in the right format. Right now our data is in the “data frame” format, but it needs to be in the “matrix” format for use in the vegan package. The following code is how to convert it
```{r}
m_com <- as.matrix(com)
m_com
```

Now we can run the metaMDS command from the vegan package to generate an NMDS plot. 

Calling your nmds object in R, will give you some information about your analysis. An important number to note is the stress, which is roughly the “goodness of fit” of your NMDS ordination. For a good representation of your data, the stress value should ideally be less than 0.2. If the stress value is 0, it might mean you have an outlier sample that is very different than all your other samples.
```{r}
set.seed(100)
nmds = metaMDS(m_com, distance = "bray")
```

You can easily plot a simple NMDS in base R:
```{r}
plot(nmds)
```

To make a nicer plot, use the ggplot2 package.
First you need to obtain the coordinates for your NMDS1 and NMDS2 axes and put them in a new data frame: I’ve called this new data frame, “data.scores”

Next, you can add columns from your original data (pc) to your new NMDS coordinates data frame. This will come in handy when you plot your data and want to differentiate groups or treatments:
```{r}
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$location = pc$location
data.scores$sample_id = pc$sample_id
head(data.scores)
write.table(data.scores, "nmds_data.scores_all_samples.csv", sep=",", quote=F, col.names=NA)
```

Make NMDS Plot using ggplot2:
```{r}
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 2, aes(colour = location))+
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
ggsave("Shipwreck_NMDS_all.tiff")
```

### NMDS Depth
```{r}
data.scores2 <- as.data.frame(scores(nmds)$sites)
data.scores2$depth = pc$depth
data.scores2$sample_id = pc$sample_id
head(data.scores2)
write.table(data.scores2, "nmds_data.scores_all_depth.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx2 = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 2, aes(depth = depth, colour=factor(depth))) +
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
xx2
ggsave("Shipwreck_NMDS_depth.tiff")
```

### NMDS Waterlevel
```{r}
data.scores3 <- as.data.frame(scores(nmds)$sites)
data.scores3$water_level = pc$water_level
data.scores3$sample_id = pc$sample_id
head(data.scores3)
write.table(data.scores3, "nmds_data.scores_water_level.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx3 = ggplot(data.scores3, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 2, aes(colour = water_level))+
 theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "water_level", y = "NMDS2")
xx3
ggsave("Shipwreck_NMDS_waterlevel.tiff")
```

## Stacked Barchart & NMDS: Starboard vs Port Side

### Prepping files into data.frame
```{r}
df1_ship <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/df1_star_port_samples_only.csv")
df1_ship
df2_ship <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/df2_star_port_samples_only.csv")
df2_ship
df_list2 <- list(df1_ship, df2_ship)
df_list2

otu_count_ship <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list2) %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")
otu_count_ship

write.table(otu_count_ship, "otu_count_ship.csv", sep=",", quote=F, col.names=NA)
write.table(nmds_asv_otu, "nmds_asv_otu.csv", sep=",", quote=F, col.names=NA)
```

### Gaining OTU Relative Abundance Values
```{r}
otu_rel_abund_ship <- inner_join(metadata_star_port_only, otu_count_ship, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Genus", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_ship
write.table(otu_rel_abund_ship, "otu_rel_abund_ship.csv", sep=",", quote=F, col.names=NA)
```

### Stacked Barchart of Starboard vs Port Side Samples Only
```{r}
otu_rel_abund_ship %>%
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
ggsave("star_port_stacked_bar.tiff", width=12, height=10)
```

### Prepping NMDS files into workable data.frame
```{r}
pc2 <- read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/nmds_asv_otu_ship.csv")
pc2
```

```{r}
com2 = pc2[,7:ncol(pc2)]
com2
```

```{r}
m_com2 <- as.matrix(com2)
m_com2
```

```{r}
set.seed(100)
nmds2 = metaMDS(m_com2, distance = "bray")
plot(nmds2)
```

```{r}
data.scores2 <- as.data.frame(scores(nmds2)$sites)
data.scores2$location = pc2$location
data.scores2$sample_id = pc2$sample_id
head(data.scores2)
write.table(data.scores2, "nmds_data.scores_star_port.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx4 = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 2, aes(shape = location, colour = location))+
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
ggsave("Shipwreck_star_port_only.tiff")
```


## Diversity Index Value Generating
```{r}
otu_table <- read.csv("/Users/maggieshostak/Desktop/Shipwreck/data/nmds_asv_otu.csv", header=T, row.names=1, check.names=FALSE)
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

## ANOSIM Test

```{r}
pc_ano = read.csv("/Users/maggieshostak/Desktop/FeOB_Shipwreck_Analysis/data/nmds_asv_otu.csv", header= TRUE)
```

Make community matrix - extract columns with abundance information, turn data frame into matrix
```{r}
com3 = pc_ano[,7:ncol(pc_ano)]
m_com3 = as.matrix(com3)
```

```{r}
ano = anosim(m_com3, pc$location, distance = "bray", permutations = 9999)
ano

ano2 = anosim(m_com3, pc$water_level, distance = "bray", permutations = 9999)
ano2

ano3 = anosim(m_com3, pc$depth, distance = "bray", permutations = 9999)
ano3
```

### Results
*Location*: (R = 0.3495, Significance = 0.0001)

*Water_Level*: (R = 0.4067, Significance = 0.0001)

*Depth*: (R = 0.4701, Significance = 0.0001)

When interpreting these results you want to look at the ANOSIM statistic "R" & the Significance values. 
  - A Significance value <0.05 is generally considered to be statistically significant, & means the null hypothesis can be rejected. Therefore, there is a statistically significant difference in the microbial communities between your groups. 
  - A Significance value >0.05, means that there is no statistical difference between the microbial communities in your groups.

The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups to the mean of ranked dissimilarities within groups. An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high & low ranks within & between groups. In other words, the higher the R value, the more dissimilar your groups are in terms of microbial community composition!

# Distance Matrix & PCA
```{r}

```
