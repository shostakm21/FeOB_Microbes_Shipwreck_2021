# Thesis Shipwreck Biofilm, Sediment and Water Samples Analysis
## Using packages vegan, ggplot2, and dada2 pipeline in R-Studio


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
library(dbplyr)
library(svglite)
library(ggstatsplot)
library(grid)
library(labdsv)
library(indicspecies)
library(viridis)
library(ggplot2)
library('cowplot')
library(permute)
library(lattice)
library(breakaway)
library(dtplyr)
```

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
```{r}
path <- "/Users/maggieshostak/Desktop/RStudio Work/fastq_files"
list.files(path)
```

```{r}
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
```

```{r}
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)
```

```{r}
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])
```

```{r}
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT"
trimLeft = c(FWD, REV)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = c(19,20))
head(out)
```

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```

```{r}
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)
```

```{r}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

```{r}
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#dadaFs[[1]]

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#dadaRs[[1]]
```

```{r}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE, verbose=TRUE)
head(mergers[[1]])
```

```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
```

```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/RStudio Work/fastq_files/silva_nr99_v138.1_train_set.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

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
row.names(asv_otu) <- sub(">", "", asv_headers)
```

```{r}
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
```

```{r}
otu_tax_table <- merge(asv_otu, asv_tax, by=0)
```

```{r}
write(asv_fasta, "asv_fasta_final.fa")
write.table(asv_otu, "asv_otu_final.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "asv_tax_final.csv", sep=",", quote=F, col.names=NA)
write.table(otu_tax_table, "otu_tax_table_final.csv", sep=",", quote=F, col.names=NA)
```

# Formating Files for Further Analysis
```{r}
metadata_all <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_biofilm_sediment_water.csv")
metadata_all

otu_counts <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_original.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts

taxonomy <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_tax.csv")
taxonomy
```

```{r}
asv_tax <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_tax_final.csv")
asv_tax
```

# Generate an OTU Relative Abundance:
```{r}
otu_rel_abund <- inner_join(metadata_all, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund

write.table(otu_rel_abund, "otu_rel_abund_final.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu_rel_abund <- read.csv("/Users/maggieshostak/Desktop/RStudio\ Work/data/otu_rel_abund_final.csv")
otu_rel_abund
```

```{r}
metadata_SSW <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_ship_sed_water.csv")
metadata_SSW
```

```{r}
metadata_SP <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_SP_biofilm.csv")
metadata_SP
```

```{r}
metadata_BS <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_BS_biofilm.csv")
metadata_BS
```

```{r}
metadata_BS_SP <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_BS_SP.csv")
metadata_BS_SP
```

```{r}
metadata_S <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_sediment.csv")
metadata_S
```

# Stacked Barcharts
## All Samples: Biofilms Separated
```{r}
## Phylum
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
#ggsave("phylum_stacked_barchart_all.tiff", width=20, height=7)

## Class
otu_rel_abund %>%
  filter(level=="Class") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("class_stacked_barchart_all.tiff", width=25, height=10)
```

## Only Biofilms Separated (SIMPER)
```{r}
## Relative Abundance
metadata_ship_simper <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_biofilm_separate.csv")
metadata_ship_simper

otu_counts_simper <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_ship_simper_sig.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts_simper

taxonomy_simper <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_tax_ship_simper_sig.csv")
taxonomy_simper

otu_rel_abund_simper <- inner_join(metadata_ship_simper, otu_counts_simper, by="sample_id") %>%
  inner_join(., taxonomy_simper, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_simper

write.table(otu_rel_abund_simper, "otu_rel_abund_simper.csv", sep=",", quote=F, col.names=NA)

## Phylum
otu_rel_abund_simper %>%
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
ggsave("phylum_stacked_barchart_simper.tiff", width=13, height=15)

## Class
otu_rel_abund_simper %>%
  filter(level=="Class") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("class_stacked_barchart_simper.tiff", width=13, height=15)

## Order
otu_rel_abund_simper %>%
  filter(level=="Order") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("order_stacked_barchart_simper.tiff", width=13, height=15)

## Family
otu_rel_abund_simper %>%
  filter(level=="Family") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("family_stacked_barchart_simper.tiff", width=20, height=15)
```

## Starboard vs Port
```{r}
metadata_SP <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_SP_biofilm.csv")
metadata_SP

## Phylum
otu_rel_abund_SP <- inner_join(metadata_SP, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_SP

write.table(otu_rel_abund_SP, "otu_rel_abund_SP.csv", sep=",", quote=F, col.names=NA)

## Class
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

ggsave("phylum_stacked_bar_SP.tiff", width=10, height=8)

otu_rel_abund_SP %>%
  filter(level=="Class") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()

ggsave("class_stacked_bar_SP.tiff", width=17, height=7)
```

## Bow vs Stern
```{r}
## Relative Abundances
metadata_BS <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_BS_biofilm.csv")
metadata_BS

otu_rel_abund_BS <- inner_join(metadata_BS, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_BS

write.table(otu_rel_abund_BS, "otu_rel_abund_BS.csv", sep=",", quote=F, col.names=NA)

## Phylum
otu_rel_abund_BS %>%
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

ggsave("phylum_stacked_bar_BS.tiff",  width=10, height=7)

## Class
otu_rel_abund_BS %>%
  filter(level=="Class") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()

ggsave("class_stacked_bar_BS.tiff",  width=15, height=7)

```

## Sediment Samples
```{r}
## Relative Abundances
metadata_S <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_sediment.csv")
metadata_S

otu_counts_S <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_sediment.csv") %>%
 pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts_S

taxonomy_S <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_tax_sediment.csv")
taxonomy_S
```

```{r}
otu_rel_abund_S <- inner_join(metadata_S, otu_counts_S, by="sample_id") %>%
  inner_join(., taxonomy_S, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_S

write.table(otu_rel_abund_S, "otu_rel_abund_S.csv", sep=",", quote=F, col.names=NA)

## Phylum
otu_rel_abund_S %>%
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

ggsave("phylum_stacked_bar_sediment.tiff", width=17, height=8)

## Class
otu_rel_abund_S %>%
  filter(level=="Class") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()

ggsave("class_stacked_bar_sediment.tiff", width=20, height=7)
```

## Depth
```{r}
metadata_depth <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_depth.csv")
metadata_depth

otu_counts_depth <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_depth.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts_depth

taxonomy_depth <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_tax_depth.csv")
taxonomy_depth

otu_rel_abund_depth <- inner_join(metadata_depth, otu_counts_depth, by="sample_id") %>%
  inner_join(., taxonomy_depth, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund_depth

write.table(otu_rel_abund_depth, "otu_rel_abund_depth.csv", sep=",", quote=F, col.names=NA)

## Phylum
otu_rel_abund_depth %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, depth, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(depth, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=depth, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=depth, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()

ggsave("phylum_stacked_bar_depth.tiff", width=10, height=7)

## Class
otu_rel_abund_depth %>%
  filter(level=="Class") %>%
  group_by(sample_id, depth, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(depth, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=depth, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=depth, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()

ggsave("class_stacked_bar_depth.tiff", width=15, height=7)
```

# NMDS
```{r}
# All Samples: Biofilm Separated
pc1 = read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/nmds_asv_otu_biofilm_sep_sediment_water.csv")
pc1

com1 = pc1[,4:ncol(pc1)]
com1

m_com1 <- as.matrix(com1)
m_com1
```

```{r}
set.seed(1000)
nmds1 = metaMDS(m_com1, distance = "bray") #stress 0.08247397 

data.scores1 <- as.data.frame(scores(nmds1)$sites)
data.scores1$location = pc1$location
data.scores1$sample_id = pc1$sample_id
head(data.scores1)

#write.table(data.scores1, "nmds_scores_biofilm_sed_water.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx1 = ggplot(data.scores1, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx1
ggsave("NMDS_biofilm_sediment_water.tiff", width = 10, height = 10)
```

```{r}
# All Samples: Shipwreck, Sediment & Water
pc1b = read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/nmds_asv_otu_shipwreck_sediment_water.csv")
pc1b

com1b = pc1b[,3:ncol(pc1b)]
com1b

m_com1b <- as.matrix(com1b)
m_com1b
```

```{r}
set.seed(1000)
nmds1b = metaMDS(m_com1b, distance = "bray") #stress 0.06829258 

data.scores1b <- as.data.frame(scores(nmds1b)$sites)
data.scores1b$location = pc1b$location
data.scores1b$sample_id = pc1b$sample_id
head(data.scores1b)

#write.table(data.scores1b, "nmds_scores_SSW.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx1b = ggplot(data.scores1b, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site By Type")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx1b
ggsave("NMDS_SSW.tiff", width = 10, height = 10)
```

```{r}
# Only Biofilm Samples
pc1c = read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/nmds_asv_otu_biofilm.csv")
pc1c

com1c = pc1c[,4:ncol(pc1c)]
com1c

m_com1c <- as.matrix(com1c)
m_com1c
```

```{r}
set.seed(1000)
nmds1c = metaMDS(m_com1c, distance = "bray") #stress 0.1418903 

data.scores1c <- as.data.frame(scores(nmds1c)$sites)
data.scores1c$location = pc1c$location
data.scores1c$sample_id = pc1c$sample_id
head(data.scores1c)

#write.table(data.scores1c, "nmds_scores_biofilm_only.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx1c = ggplot(data.scores1c, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx1c
ggsave("NMDS_biofilm_only.tiff", width = 10, height = 10)
```

## Star vs Port
```{r}
# Starboard vs Port
pc2 = read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/nmds_asv_otu_SP_biofilm.csv")
pc2

com2 = pc2[,4:ncol(pc2)]
com2

m_com2 <- as.matrix(com2)
m_com2
```

```{r}
set.seed(1000)
nmds2 = metaMDS(m_com2, distance = "bray") #stress 0.1809246 

data.scores2 <- as.data.frame(scores(nmds2)$sites)
data.scores2$location = pc2$location
data.scores2$sample_id = pc2$sample_id
head(data.scores2)

#write.table(data.scores2, "nmds_scores_biofilm_SP.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx2 = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx2
ggsave("NMDS_biofilm_SP.tiff", width = 10, height = 10)
```

## Bow vs Stern
```{r}
# Bow vs Stern
pc3 = read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/nmds_asv_otu_BS_biofilm.csv")
pc3

com3 = pc3[,4:ncol(pc3)]
com3

m_com3 <- as.matrix(com3)
m_com3
```

```{r}
set.seed(1000)
nmds3 = metaMDS(m_com3, distance = "bray") #stress 0.1702069 

data.scores3 <- as.data.frame(scores(nmds3)$sites)
data.scores3$location = pc3$location
data.scores3$sample_id = pc3$sample_id
head(data.scores3)

#write.table(data.scores3, "nmds_scores_biofilm_BS.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx3 = ggplot(data.scores3, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx3
ggsave("NMDS_biofilm_BS.tiff", width = 10, height = 10)
```

## Starboard vs Port vs Bow vs Stern
```{r}
# Starboard vs Port vs Bow vs Stern
pc4 = read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/nmds_asv_otu_BS_SP_biofilm.csv")
pc4

com4 = pc4[,4:ncol(pc4)]
com4

m_com4 <- as.matrix(com4)
m_com4
```

```{r}
set.seed(1000)
nmds4 = metaMDS(m_com4, distance = "bray") #stress 0.1702069 

data.scores4 <- as.data.frame(scores(nmds4)$sites)
data.scores4$location = pc4$location
data.scores4$sample_id = pc4$sample_id
head(data.scores4)

#write.table(data.scores4, "nmds_scores_biofilm_SP_BS.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx4 = ggplot(data.scores4, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
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
ggsave("NMDS_biofilm_BS_SP.tiff", width = 10, height = 10)
```

## Sediment
```{r}
# Sediment Samples
pc5 = read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/nmds_asv_otu_sediment.csv")
pc5

com5 = pc5[,3:ncol(pc5)]
com5

m_com5 <- as.matrix(com5)
m_com5
```

```{r}
set.seed(1000)
nmds5 = metaMDS(m_com5, distance = "bray") #stress 0.08533898

data.scores5 <- as.data.frame(scores(nmds5)$sites)
data.scores5$location = pc5$location
data.scores5$sample_id = pc5$sample_id
head(data.scores5)

#write.table(data.scores5, "nmds_scores_biofilm_SEDIMENT.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx5 = ggplot(data.scores5, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx5
ggsave("NMDS_biofilm_sediment.tiff", width = 10, height = 10)
```

## Depth
```{r}
# Depth
pc6 = read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/nmds_asv_otu_depth.csv")
pc6

com6 = pc6[,3:ncol(pc6)]
com6

m_com6 <- as.matrix(com6)
m_com6
```

```{r}
set.seed(1000)
nmds6 = metaMDS(m_com6, distance = "bray") #stress 0.2003507 

data.scores6 <- as.data.frame(scores(nmds6)$sites)
data.scores6$depth = pc6$depth
data.scores6$sample_id = pc6$sample_id
head(data.scores6)

#write.table(data.scores6, "nmds_scores_biofilm_depth.csv", sep=",", quote=F, col.names=NA)
```

```{r}
xx6 = ggplot(data.scores6, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = depth))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples At Varying Depths")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "depth", y = "NMDS2")
xx6
ggsave("NMDS_biofilm_depth.tiff", width = 10, height = 10)
```

# Statistical Analyses
## Anosim
```{r}
#Location Specific Biofilm: Since data must be numeric, a number was assigned to each different location- Port(1), Starboard(2), Aft_Star(3), Bulkhead(4), and Rudder Post(5)
pcL <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/ANOSIM_Location_Biofilm.csv")
pcL

comL = pcL[,3:ncol(pcL)]
comL

m_comL = as.matrix(comL)
m_comL

anoL = anosim(m_comL, pcL$location, distance = "bray", permutations = 9999)
anoL
```
Call:
anosim(x = m_comL, grouping = pcL$location, permutations = 9999, distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.4937 
      Significance: 1e-04 

Permutation: free
Number of permutations: 9999

```{r}

```{r}
#Depth Specific Biofilm
pcD <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/ANOSIM_Depth_Biofilm.csv")
pcD

comD = pcD[,3:ncol(pcD)]
comD

m_comD = as.matrix(comD)
m_comD

anoD = anosim(m_comD, pcD$depth, distance = "bray", permutations = 9999)
anoD
```
Call:
anosim(x = m_comD, grouping = pcD$depth, permutations = 9999, distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.5569 
      Significance: 1e-04 

Permutation: free
Number of permutations: 9999

```{r}
#Depth: Waterline vs Below Waterline (Star & Port)
pcDW <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/ANOSIM_Depth_SP_Biofilm.csv")
pcDW

comDW = pcDW[,3:ncol(pcDW)]
comDW

m_comDW = as.matrix(comDW)
m_comDW

anoDW = anosim(m_comDW, pcDW$depth, distance = "bray", permutations = 9999)
anoDW
```
Call:
anosim(x = m_comDW, grouping = pcDW$depth, permutations = 9999,      distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.4665 
      Significance: 1e-04 

Permutation: free
Number of permutations: 9999

```{r}
#Sediment & Depth vs Below waterline vs Waterline: (5ft vs 3ft vs 0.3ft vs 0ft)
pcS <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/Anosim_sediment.csv")
pcS

comS = pcS[,3:ncol(pcS)]
comS

m_comS = as.matrix(comS)
m_comS

anoS = anosim(m_comS, pcS$depth, distance = "bray", permutations = 9999)
anoS
```
Call:
anosim(x = m_comS, grouping = pcS$depth, permutations = 9999,      distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.6182 
      Significance: 1e-04 

Permutation: free
Number of permutations: 9999

```{r}
# Location: Starboard vs Port Biofilm
pcSP <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/ANOSIM_Location_SP_Biofilm.csv")
pcSP

comSP = pcSP[,3:ncol(pcSP)]
comSP

m_comSP = as.matrix(comSP)
m_comS

anoSP = anosim(m_comSP, pcSP$location, distance = "bray", permutations = 9999)
anoSP
```
Call:
anosim(x = m_comSP, grouping = pcSP$location, permutations = 9999,      distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.2887 
      Significance: 3e-04 

Permutation: free
Number of permutations: 9999

```{r}
#Bow vs Stern (Port vs Star) Specific Biofilm
pcBSSP <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/ANOSIM_bow_stern.csv")
pcBSSP

comBSSP = pcBSSP[,3:ncol(pcBSSP)]
comBSSP

m_comBSSP = as.matrix(comBSSP)
m_comBSSP

anoBSSP = anosim(m_comBSSP, pcBSSP$location, distance = "bray", permutations = 9999)
anoBSSP
```
Call:
anosim(x = m_comBS, grouping = pcBS$location, permutations = 9999,      distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.1684 
      Significance: 0.0283 

Permutation: free
Number of permutations: 9999

```{r}
#Bow vs Stern Specific Biofilm
pcBS <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/ANOSIM_BS_only.csv")
pcBS

comBS = pcBS[,3:ncol(pcBS)]
comBS

m_comBS = as.matrix(comBS)
m_comBS

anoBS = anosim(m_comBS, pcBS$location, distance = "bray", permutations = 9999)
anoBS
```
Call:
anosim(x = m_comBS, grouping = pcBS$location, permutations = 9999,      distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.02141 
      Significance: 0.297 

Permutation: free
Number of permutations: 9999

```{r}
#Type of Samples: Biofilm (1),Sediment(2) & Water(3)
pcBSW <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/Anosim_Ship_Sed_Water.csv")
pcBSW

comBSW = pcBSW[,3:ncol(pcBSW)]
comBSW

m_comBSW = as.matrix(comBSW)
m_comBSW

anoBSW = anosim(m_comBSW, pcBSW$location, distance = "bray", permutations = 9999)
anoBSW
```
Call:
anosim(x = m_comBSW, grouping = pcBSW$location, permutations = 9999,      distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.7852 
      Significance: 1e-04 

Permutation: free
Number of permutations: 9999

```{r}
# Sediment
pcSed <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/ANOSIM_sediment_distance.csv")
pcSed

comSed = pcSed[,3:ncol(pcSed)]
comSed

m_comSed = as.matrix(comSed)
m_comSed

anoSed = anosim(m_comSed, pcSed$location, distance = "bray", permutations = 9999)
anoSed
```
Call:
anosim(x = m_comSed, grouping = pcSed$location, permutations = 9999, distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.0934 
      Significance: 0.2213 

Permutation: free
Number of permutations: 9999

# Diversity Index Value Generating
```{r}
## Biofilm & Sediment & Water Samples
otu_table <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_original.csv", header=T, row.names=1, check.names=FALSE)

## Transpose the data to have sample names on rows
otu.table.diver <- t(otu_table)
otu.table.diver <- as.data.frame(otu.table.diver)
head(otu.table.diver)
```

```{r}
## Biofilm Only
otu_table_bio <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_biofilm.csv", header=T, row.names=1, check.names=FALSE)
otu_table_bio

## Transpose the data to have sample names on rows
otu.table.diver.bio <- t(otu_table_bio)
otu.table.diver.bio <- as.data.frame(otu.table.diver.bio)
head(otu.table.diver.bio)
```

# Ecological Distance Matrices & Rarefaction
```{r}
library(tidyverse)

df1 <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/df_1.csv")
df1
df2 <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/df_2.csv")
df2
df_list <- list(df1, df2)
df_list

otu_count_all <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list) %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")

#write.table(otu_count_all, "otu_count_all.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu_count_all %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 5000) %>%
  group_by(ASV) %>%
  mutate(total=sum(count)) %>% 
  filter(total != 0) %>%
  as.data.frame()
#Going to set threshold at 5000
```

lowest total read count make sure its even

min(rowSums(otu_final))
otus.r <- rrarefy(otu_final, 546)
fisher ,- fisher.alpha(otus.r)



## Shannon
```{r}
data(otu.table.diver)
H<-diversity(otu.table.diver)
H
```

## Richness
```{r}
richness <- specnumber(otu.table.diver)
richness
```

## Pielou Evenness
```{r}
evenness <- H/log(richness)
evenness
```

```{r}
metadata_all <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_biofilm_sediment_water.csv")
metadata_all
```

```{r}
alpha <- cbind(shannon = H, richness = richness, pielou = eveness, metadata_all)
write.csv(alpha, "diversity_indices_bio_sed_water.csv")
head(alpha)
```


```{r}
plot.shan <- ggplot(alpha, aes(x = location, y = shannon, colour = location)) +
geom_point(size = 3) +
ylab("Shannon's H'") + 
xlab("") +
ggtitle("Shannon's Diversity - Samples Across Site")+
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan
ggsave("Shannon_Location_bio_sed_water.tiff")
```


```{r}
plot.rich <-ggplot(alpha, aes(x = location, y = richness, colour = location)) +
geom_point(size = 3) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich
ggsave("Richness_Location_bio_sed_water.tiff")
```


```{r}
plot.even <- ggplot(alpha, aes(x = location, y = pielou, colour = location)) +
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

ggsave("Shannon_Richness_Eveness_bio_sed_water.tiff")
```

## Biofilm Only
```{r}
otu_table_bio <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/asv_otu_biofilm.csv", header=T, row.names=1, check.names=FALSE)
head(otu_table_bio)
```

Transpose the data to have sample names on rows
```{r}
otu.table.diver.bio <- t(otu_table_bio)
otu.table.diver.bio <- as.data.frame(otu.table.diver.bio)
head(otu.table.diver.bio)
write.csv(otu.table.diver.bio, "otu.table.diversity.biofilm.csv")
```

```{r}
data(otu.table.diver.bio)
H<-diversity(otu.table.diver.bio)
H

richness <- specnumber(otu.table.diver.bio)

eveness <- H/log(richness)
```

```{r}
metadata_bio <- read.csv("/Users/maggieshostak/Desktop/RStudio Work/data/metadata_biofilm_separate.csv")
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

## Simper
```{r}
otu.table.diver.mdf.bio <- as.matrix.data.frame(otu.table.diver.bio)
rownames(otu.table.diver.mdf.bio) <- metadata_bio$location

otu.table.diver.bray.bio <- vegdist(otu.table.diver.mdf.bio, method="bray")
otu.table.diver.bray.bio
```

```{r}
#Location
simper <- simper(otu.table.diver.bio, metadata_bio$location, permutations=999)
options(max.print=999999)
summary(simper)
dput(simper, file = "simp_location.txt")
sim <- dget("/Users/maggieshostak/Desktop/RStudio Work/data/simp_location.txt")
summary(sim)
```

# Rarefaction
```{r}
library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
library(phyloseq) ; packageVersion("phyloseq") # 1.22.3
library(vegan) ; packageVersion("vegan") # 2.5.4
library(DESeq2) ; packageVersion("DESeq2") # 1.18.1
library(viridis) ; packageVersion("viridis") # 0.5.1
library(phyloseq)
library(ggplot2)
library(devtools)
library(ggpubr)
```

# TESTING ALL SAMPLES
```{r}
count_tab <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_otu_original_copy.csv", header=T, row.names=1, check.names=F, sep=",")
count_tab

tax_tab <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_tax_original_copy.csv", header=T, row.names=1, check.names=F, sep=",")
tax_tab

metadata <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/metadata_copy.csv")
metadata
```

```{r}
deseq_counts <- DESeqDataSetFromMatrix(countData = count_tab, colData = metadata, design = ~sample_id) 
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab))

euc_clust <- hclust(euc_dist, method="ward.D2")
plot(euc_clust) 
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(metadata$color[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.", width = 35, height = 35 )
```

```{r}
rarecurve(t(count_tab), step=500, col=metadata$color, lwd=2, ylab="ASVs", label=F, cex=0.6)
  abline(v=(min(rowSums(t(count_tab)))), col = "red")
```

# TESTING WITHOUT SAMPLES 16 & 48 & 51
```{r}
count_tab_a <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_otu_no_16_48.csv", header=T, row.names=1, check.names=F, sep=",")
count_tab_a

tax_tab_a <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_tax_original_copy.csv", header=T, row.names=1, check.names=F, sep=",")
tax_tab_a

metadata_a <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/metadata_no_16_48.csv")
metadata_a
```

```{r}
deseq_counts_a <- DESeqDataSetFromMatrix(countData = count_tab_a, colData = metadata_a, design = ~sample_id) 
deseq_counts_a <- estimateSizeFactors(deseq_counts_a, type = "poscounts")
deseq_counts_vst_a <- varianceStabilizingTransformation(deseq_counts_a)

vst_trans_count_tab_a <- assay(deseq_counts_vst_a)

euc_dist_a <- dist(t(vst_trans_count_tab_a))

euc_clust_a <- hclust(euc_dist_a, method="ward.D2")
plot(euc_clust_a) 

euc_dend_a <- as.dendrogram(euc_clust_a, hang=0.1)
dend_cols_a <- as.character(metadata_a$color[order.dendrogram(euc_dend_a)])
labels_colors(euc_dend_a) <- dend_cols_a

plot(euc_dend_a, ylab="VST Euc. dist.", width = 35, height = 35 )
```

```{r}
rarecurve(t(count_tab_a), step=500, lwd=2, col=metadata_a$color, ylab="ASVs", label=T, cex=0.6)
  abline(v=(min(rowSums(t(count_tab_a)))), col = "red")
```

# TESTING WITHOUT SAMPLES 16 & 48 & 51
```{r}
count_tab_b <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_otu_no_16_48_51.csv", header=T, row.names=1, check.names=F, sep=",")
count_tab_b

tax_tab_b <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_tax_original_copy.csv", header=T, row.names=1, check.names=F, sep=",")
tax_tab_b

metadata_b <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/metadata_no_16_48_51.csv")
metadata_b
```

```{r}
deseq_counts_b <- DESeqDataSetFromMatrix(countData = count_tab_b, colData = metadata_b, design = ~sample_id) 
deseq_counts_b <- estimateSizeFactors(deseq_counts_b, type = "poscounts")
deseq_counts_vst_b <- varianceStabilizingTransformation(deseq_counts_b)

vst_trans_count_tab_b <- assay(deseq_counts_vst_b)

euc_dist_b <- dist(t(vst_trans_count_tab_b))

euc_clust_b <- hclust(euc_dist_b, method="ward.D2")
plot(euc_clust_b) 

euc_dend_b <- as.dendrogram(euc_clust_b, hang=0.1)
dend_cols_b <- as.character(metadata_b$color[order.dendrogram(euc_dend_b)])
labels_colors(euc_dend_b) <- dend_cols_b

plot(euc_dend_b, ylab="VST Euc. dist.", width = 35, height = 35 )
```

```{r}
rarecurve(t(count_tab_b), step=500, lwd=2, col=metadata_b$color, ylab="ASVs", label=T, cex=0.6)
  abline(v=(min(rowSums(t(count_tab_b)))), col = "red")
```

# TESTING BASED ON SAMPLE TYPE: Biofilm
```{r}
count_tab_B <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_otu_biofilm.csv", header=T, row.names=1, check.names=F, sep=",")
count_tab_B

tax_tab_B <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_tax_original_copy.csv", header=T, row.names=1, check.names=F, sep=",")
tax_tab_B

metadata_B <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/metadata_biofilm.csv")
metadata_B
```

```{r}
deseq_counts_B <- DESeqDataSetFromMatrix(countData = count_tab_B, colData = metadata_B, design = ~sample_id) 
deseq_counts_B <- estimateSizeFactors(deseq_counts_B, type = "poscounts")
deseq_counts_vst_B <- varianceStabilizingTransformation(deseq_counts_B)

vst_trans_count_tab_B <- assay(deseq_counts_vst_B)

euc_dist_B <- dist(t(vst_trans_count_tab_B))

euc_clust_B <- hclust(euc_dist_B, method="ward.D2")
plot(euc_clust_B) 

euc_dend_B <- as.dendrogram(euc_clust_B, hang=0.1)
dend_cols_B <- as.character(metadata_B$color[order.dendrogram(euc_dend_B)])
labels_colors(euc_dend_B) <- dend_cols_B

plot(euc_dend_B, ylab="VST Euc. dist.", width = 35, height = 35 )
```

```{r}
rarecurve(t(count_tab_B), step=500, lwd=2, col=metadata_B$color, ylab="ASVs", label=T, cex=0.6)
  abline(v=(min(rowSums(t(count_tab_B)))), col = "red")
```

# TESTING BASED ON SAMPLE TYPE: Sediment
```{r}
count_tab_S <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_otu_sed.csv", header=T, row.names=1, check.names=F, sep=",")
count_tab_S

tax_tab_S <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_tax_original_copy.csv", header=T, row.names=1, check.names=F, sep=",")
tax_tab_S

metadata_S <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/metadata_sed.csv")
metadata_S
```

```{r}
deseq_counts_S <- DESeqDataSetFromMatrix(countData = count_tab_S, colData = metadata_S, design = ~sample_id) 
deseq_counts_S <- estimateSizeFactors(deseq_counts_S, type = "poscounts")
deseq_counts_vst_S <- varianceStabilizingTransformation(deseq_counts_S)

vst_trans_count_tab_S <- assay(deseq_counts_vst_S)

euc_dist_S <- dist(t(vst_trans_count_tab_S))

euc_clust_S <- hclust(euc_dist_S, method="ward.D2")
plot(euc_clust_S) 

euc_dend_S <- as.dendrogram(euc_clust_S, hang=0.1)
dend_cols_S <- as.character(metadata_S$color[order.dendrogram(euc_dend_S)])
labels_colors(euc_dend_S) <- dend_cols_S

plot(euc_dend_S, ylab="VST Euc. dist.", width = 35, height = 35 )
```

```{r}
rarecurve(t(count_tab_S), step=500, lwd=2, col=metadata_S$color, ylab="ASVs", label=T, cex=0.6)
  abline(v=(min(rowSums(t(count_tab_S)))), col = "red")
```

# TESTING BASED ON SAMPLE TYPE: Water
```{r}
count_tab_W <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_otu_water.csv", header=T, row.names=1, check.names=F, sep=",")
count_tab_W

tax_tab_W <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/asv_tax_original_copy.csv", header=T, row.names=1, check.names=F, sep=",")
tax_tab_W

metadata_W <- read.csv("/Users/maggieshostak/Desktop/rarefaction_testing/metadata_water.csv")
metadata_W
```

```{r}
deseq_counts_W <- DESeqDataSetFromMatrix(countData = count_tab_W, colData = metadata_W, design = ~sample_id) 
deseq_counts_W <- estimateSizeFactors(deseq_counts_W, type = "poscounts")
deseq_counts_vst_W <- varianceStabilizingTransformation(deseq_counts_W)

vst_trans_count_tab_W <- assay(deseq_counts_vst_W)

euc_dist_W <- dist(t(vst_trans_count_tab_W))

euc_clust_W <- hclust(euc_dist_W, method="ward.D2")
plot(euc_clust_W) 

euc_dend_W <- as.dendrogram(euc_clust_W, hang=0.1)
dend_cols_W <- as.character(metadata_W$color[order.dendrogram(euc_dend_W)])
labels_colors(euc_dend_W) <- dend_cols_W

plot(euc_dend_W, ylab="VST Euc. dist.", width = 35, height = 35 )
```

```{r}
rarecurve(t(count_tab_W), step=500, lwd=2, col=metadata_W$color, ylab="ASVs", label=T, cex=0.6)
  abline(v=(min(rowSums(t(count_tab_W)))), col = "red")
```

## NMDS Ordination Plot Check
```{r}
df1 <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/df_1.csv")
df1
df2 <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/df_2.csv")
df2
df_otu_5000 <- list(df1, df2)
df_otu_5000

#write.table(df_otu_5000, "/Users/maggieshostak/Desktop/Mallows_5000_cutoff/df_otu_5000.csv", sep=",", quote=F, col.names=NA)

otu_count_5000 <- Reduce(function(x, y) merge(x, y, all=TRUE), df_otu_5000) %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")

#write.table(otu_count_5000, "otu_count_5000.csv", sep=",", quote=F, col.names=NA)

otu_count_5000 %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 5000) %>%
  group_by(ASV) %>%
  mutate(total=sum(count)) %>% 
  filter(total != 0) %>%
  as.data.frame()
#Going to set threshold at 5000
```

```{r}
df_meta_5000 <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/metadata_5000_cutoff.csv")
df_meta_5000

df_otu_5000 <-read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/df_otu_5000.csv", row.names=1 )
df_otu_5000

nmds_asv_otu_all_5000 <- inner_join(df_meta_5000, df_otu_5000, by="sample_id")

write.table(nmds_asv_otu_all_5000, "/Users/maggieshostak/Desktop/Mallows_5000_cutoff/nmds_asv_otu_5000.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## All Samples: By Location
pcX <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/nmds_asv_otu_biofilm_sep_sed_water_5000.csv")
pcX
```

```{r}
comX = pcX[,4:ncol(pcX)]
comX

m_comX <- as.matrix(comX)
#m_comX
```

```{r}
set.seed(1000)
nmdsX = metaMDS(m_comX, distance = "bray")

data.scoresX <- as.data.frame(scores(nmdsX)$sites)
data.scoresX$location = pcX$location
data.scoresX$sample_id = pcX$sample_id
head(data.scoresX)

xxX = ggplot(data.scoresX, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xxX
ggsave("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/NMDS_all_samples_location_5000.tiff", width = 10, height = 10)
```

## NMDS Plot Check: Star vs Port
```{r}
df1_SP <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/df_1_SP.csv")
df1_SP

df2_SP <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/df_2_SP.csv")
df2_SP

df_otu_5000_SP <- list(df1_SP, df2_SP)
df_otu_5000_SP

write.table(df_otu_5000_SP, "/Users/maggieshostak/Desktop/Mallows_5000_cutoff/df_otu_5000_SP.csv", sep=",", quote=F, col.names=NA)
```

```{r}
df_meta_5000_SP <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/metadata_5000_SP.csv")
df_meta_5000_SP

df_otu_5000_SP <-read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/df_otu_5000_SP.csv", row.names=1 )
df_otu_5000_SP

nmds_asv_otu_SP_5000 <- inner_join(df_meta_5000_SP, df_otu_5000_SP, by="sample_id")

write.table(nmds_asv_otu_SP_5000, "/Users/maggieshostak/Desktop/Mallows_5000_cutoff/nmds_asv_otu_SP_5000.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## All Samples: By Location
pcX2 <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/nmds_asv_otu_SP_5000.csv")
pcX2
```

```{r}
comX2 = pcX2[,4:ncol(pcX2)]
comX2

m_comX2 <- as.matrix(comX2)
#m_comX2
```

```{r}
set.seed(1000)
nmdsX2 = metaMDS(m_comX2, distance = "bray")

data.scoresX2 <- as.data.frame(scores(nmdsX2)$sites)
data.scoresX2$location = pcX2$location
data.scoresX2$sample_id = pcX2$sample_id
head(data.scoresX2)

xxX2 = ggplot(data.scoresX2, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xxX2
ggsave("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/NMDS_SP_location_5000.tiff", width = 10, height = 10)
```

## NMDS Plot Check: PB vs PS vs SB vs SS
```{r}
df_meta_5000_SP_BS <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/metadata_5000_SP_BS.csv")
df_meta_5000_SP_BS

df_otu_5000_SP_BS <-read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/df_otu_5000.csv", row.names=1 )
df_otu_5000_SP_BS

nmds_asv_otu_SP_BS_5000 <- inner_join(df_meta_5000_SP_BS, df_otu_5000_SP_BS, by="sample_id")

write.table(nmds_asv_otu_SP_BS_5000, "/Users/maggieshostak/Desktop/Mallows_5000_cutoff/nmds_asv_otu_SP_BS_5000.csv", sep=",", quote=F, col.names=NA)
```

```{r}
## All Samples: By Location
pcX3 <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/nmds_asv_otu_SP_BS_5000.csv")
pcX3
```

```{r}
comX3 = pcX3[,4:ncol(pcX3)]
comX3

m_comX3 <- as.matrix(comX3)
#m_comX3
```

```{r}
set.seed(1000)
nmdsX3 = metaMDS(m_comX3, distance = "bray")

data.scoresX3 <- as.data.frame(scores(nmdsX3)$sites)
data.scoresX3$location = pcX3$location
data.scoresX3$sample_id = pcX3$sample_id
head(data.scoresX3)

xxX3 = ggplot(data.scoresX3, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xxX3
ggsave("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/NMDS_SP_BS_location_5000.tiff", width = 10, height = 10)
```

# ANOSIM: Statistical Test
```{r}
#All Samples
pcQ <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/nmds_asv_otu_biofilm_sep_sed_water_5000.csv")
pcQ

comQ = pcQ[,4:ncol(pcQ)]
comQ

m_comQ = as.matrix(comQ)
m_comQ

anoQ = anosim(m_comQ, pcQ$location, distance = "bray", permutations = 9999)
anoQ
```
ANOSIM statistic R: 0.8448 
      Significance: 1e-04 -> 0.0001

```{r}
# Star vs Port Biofilm
pcR <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/nmds_asv_otu_SP_5000.csv")
pcR

comR = pcR[,4:ncol(pcR)]
comR

m_comR = as.matrix(comR)
m_comR

anoR = anosim(m_comR, pcR$location, distance = "bray", permutations = 9999)
anoR
```
ANOSIM statistic R: 0.2479 
      Significance: 9e-04 -> 0.0009

```{r}
# Star vs Port vs Bow vs Stern Biofilm
pcT <- read.csv("/Users/maggieshostak/Desktop/Mallows_5000_cutoff/nmds_asv_otu_SP_BS_5000.csv")
pcT

comT = pcT[,4:ncol(pcT)]
comT

m_comT = as.matrix(comT)
m_comT

anoT = anosim(m_comT, pcT$location, distance = "bray", permutations = 9999)
anoT
```
ANOSIM statistic R: 0.04828 
      Significance: 0.2456 









