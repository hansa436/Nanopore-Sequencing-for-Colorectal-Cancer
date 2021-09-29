# Methylation

## Megalodon


script all run on sample crc669
Sample CRC669 was used as an example, 
Each script is altered depending on the sample but base code for the command/process is the same
Run in bash unless stated otherwise

Megalodon needs to be run in a python envrionment

```{bash, eval=F}
#!/bin/bash

# Add ont stuff to path
PATH=$PATH:/Volumes/archive/blacklab/SarahHannah/ont-guppy-gpu/bin/


#Path to fast5
FAST5=/Volumes/archive/blacklab/SarahHannah/HCS-transfer/20210622_mblack/crc669/20210622_1456_X5_FAP36895_73fbf969/fast5_pass
#Output directory
OUTPUT=/Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc669/

megalodon \
	$FAST5 \
	--guppy-server-path /Volumes/archive/blacklab/SarahHannah/ont-guppy-gpu/bin/guppy_basecall_server \
	--outputs basecalls mappings mod_mappings mods \
	--overwrite \
	--output-directory $OUTPUT \
	--reference genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --mod-motif m CG 0 \
	--devices 0 --processes 40
```

this creates a bed file with all the methylation (modified_bases.5mC.bed)

Need to remove the NAs from the array (HM450.hg38_manifest.bed) bed file

```{bash, eval=F}
grep -v "NA" /Volumes/archive/blacklab/SarahHannah/crc-data/crc-meth-data/HM450.hg38.manifest.bed > /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/HM450_cut_NA.hg38_manifest.bed

```


## Intersect the megalodon bed with the 450k array bed

```{bash, eval=F}
#gives read coverage
awk '{print "chr"$0}' /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc636/modified_bases.5mC.bed > /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc669/modified_bases_labled.5mC.bed

#produces coverage
bedtools intersect -a /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc669/modified_bases_labled.5mC.bed -b /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/HM450_cut_NA.hg38_manifest.bed> /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc669/intercept_669

#extract wanted reads
cut -f1-3,5-6,11  /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc636/intersect_636 > /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc636/intercept_cov

#produces cpg site
bedtools intersect -a bed -b /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/HM450_cut_NA.hg38_manifest.bed -b/Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc669/modified_bases_labled.5mC.bed> /Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc669/intercept_cpg

#combine coverage and cpg sites
cut -f5 intercept_cpg.bed | paste intersect_cov.bed -> intercept_636.bed

```

produces a file containing both the methylation of the CpG site by nanopore and the annotation of the CpG site




## Extract the clinical data

Load in the clinical data
```{r}
clin_meth<-readRDS("/Volumes/archive/blacklab/SarahHannah/crc-data/crc-meth-data/crc-meth-data.RDS")
clin_meth <- as.data.frame(clin_meth)
#class(clin_meth)
```

Extract the Beta values
```{r}
crc669_clin_B_values<- clin_meth[,"CRC669T", drop=F]
#write.table(crc669_clin_B_values, file = "/Volumes/archive/blacklab/SarahHannah/crc-data/crc-meth-data/crc669_clin_B-values")

crc636_clin_B_values<- clin_meth[,"CRC636T", drop=F]
#write.table(crc636_clin_B_values, file = "/Volumes/archive/blacklab/SarahHannah/crc-data/crc-meth-data/crc636_clin_B-values")

crc557_clin_B_values<- clin_meth[,"CRC557T", drop=F]
#write.table(crc557_clin_B_values, file = "/Volumes/archive/blacklab/SarahHannah/crc-data/crc-meth-data/crc557_clin_B-values")
```

Remove the cpg sites from being rownames so they can be used as the merging variable
```{r}
crc669_clin_B_values <- cbind(rownames(crc669_clin_B_values), data.frame(crc669_clin_B_values, row.names=NULL))
colnames(crc669_clin_B_values) <- c("cpg", "beta")
crc636_clin_B_values <- cbind(rownames(crc636_clin_B_values), data.frame(crc636_clin_B_values, row.names=NULL))
colnames(crc636_clin_B_values) <- c("cpg", "beta")
crc557_clin_B_values <- cbind(rownames(crc557_clin_B_values), data.frame(crc557_clin_B_values, row.names=NULL))
colnames(crc557_clin_B_values) <- c("cpg", "beta")
```


## R-studio Methylation analysis

load packages
```
library(data.table)
library(ggplot2)
library(DEGreport)
```
### CRC669

Methylation for CRC669 at all readdepths

```{r}
# orignial intercept
int_meth_669<- fread("/Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc669/intercept_669")
colnames(int_meth_669) <- c("chr", "start", "end", "read_depth", "strand", "coverage", "cpg")

#convert coverage on a 0-1 scale for nanopore to be comparable with array
int_meth_669$coverage <- int_meth_669$coverage/100
head(int_meth_669)

# merge the B-values to the rest of the file based on CpG site
int_meth_669<- merge(int_meth_669, crc669_clin_B_values,  by="cpg")
head(int_meth_669)

# plot the methylatyion comparison for NAnopore vs CpG
#ggplot(int_meth_669, mapping = aes(x=coverage, y=beta)) + geom_point() + ggtitle("full methylation: array vs nanopore for crc669") + labs(x="nanopore", y= "array")
full_meth_669 <- ggplot(int_meth_669, ,mapping = aes(x=coverage, y=beta, colour=read_depth)) + geom_point() + ggtitle("full methylation: array vs nanopore for crc669") + labs(x="Methylation at each CpG site: Nanopore", y= "Methylation at each CpG site: Array")
print(full_meth_669 + scale_colour_gradient(low="light blue", high = "dark blue")) + geom_cor(method = "spearman")
```

Subset methylation based on read depth >=8 

```{r}
# subset via read depth
cov_8_669 <- subset(int_meth_669, read_depth >= 8)
head(cov_8_669)

#plot nanopore vs array

ggplot(cov_8_669, mapping=aes(x=coverage, y=beta)) + geom_point() + labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle( "methylated regions for crc669 with read depth >=8: nanopore vs array") + geom_cor()

ggplot(cov_8_669, mapping = aes(x= coverage, y= beta)) + geom_violin() + geom_point() +  labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc669 with read depth >=8: nanopore vs array")

ggplot(cov_8_669, mapping = aes(x= coverage, y= beta)) + geom_jitter(width = 0.01, alpha=0.5) +  labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc669 with read depth >=8: nanopore vs array")

#ggplot(cov_8_669, mapping = aes(x= factor(coverage), y= beta)) +  geom_violin() + geom_point(size= 0.1) +
  #labs(x= "methylated regions for crc669 with read depth >=8: nanopore vs array") + ggtitle("coverage of methylated regions for sample crc636 comparing nanopore vs array")
```

### CRC636

compare methylation for CRC636 with all read depth

```{r}
intercept_636<- fread("/Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc636/intersect_636")
colnames(intercept_636) <- c("chr", "start", "end", "read_depth", "strand", "coverage", "cpg")
intercept_636$coverage <- intercept_636$coverage/100

intercept_636<-merge(intercept_636, crc636_clin_B_values, by= "cpg")
head(intercept_636)

#ggplot(intercept_636, ,mapping = aes(x=coverage, y=beta, colour=read_depth)) + geom_point() + ggtitle("full methylation: array vs nanopore for crc636") + labs(x="nanopore", y= "array")

full_meth_636 <- ggplot(intercept_636, ,mapping = aes(x=coverage, y=beta, colour=read_depth)) + geom_point() + ggtitle("Methylation: array vs nanopore for crc636") + labs(x="Methylation at each CpG site:Nanopore", y= "Methylation at each CpG site: Array")
print(full_meth_636 + scale_colour_gradient(low="light blue", high = "dark blue")) + geom_cor()
```



Compare methylation for CRC636 with readdepth >8

```{r}
#subset read depth
cpg_636 <- subset(intercept_636, read_depth >= 8)
head(cpg_636)

Plot nanopore vs array

#ggplot(cpg_636, mapping = aes(x= coverage, y= beta )) + geom_point() + labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc636 with read depth >8: nanopore vs array") 

point636<-ggplot(cpg_636, mapping = aes(x= coverage, y= beta, colour = read_depth )) + geom_point() + labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc636 with read depth >=8: nanopore vs array") + geom_cor()

print(point636 + scale_colour_gradient(low= "light pink", high= "dark red"))

ggplot(cpg_636, mapping = aes(x= coverage, y= beta)) + geom_violin() + geom_point(alpha=0.1) +  labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc636 with read depth >=8: nanopore vs array")

ggplot(cpg_636, mapping = aes(x= coverage, y= beta)) + geom_jitter(width = 0.01, alpha=0.5) +  labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc636 with read depth >=8: nanopore vs array")

#ggplot(cpg_636, mapping = aes(x= factor(coverage), y= beta)) +  geom_violin() + geom_point(size= 0.1) + labs(x= "methylated regions for crc636 with read depth >=8: nanopore vs array") + ggtitle("coverage of methylated regions for sample crc636 comparing nanopore vs array")

```

Compare methylation for CRC636 with readdepth >=12
```{r}
cov_12_636 <- subset(cpg_636, read_depth >=12)
head(cov_12_636)
ggplot(cov_12_636, mapping = aes(x= coverage, y= beta)) + geom_point() + labs(x="Methylation at each CpG site: Nanopore", y= "Methylation at each CpG site: Array") + ggtitle("Methylation for crc636 with read depth >=12: Nanopore vs Array") + geom_cor()

```

## CRC557

Plot methylation crc557 for all read depth
```{r}

intercept_557<- fread("/Volumes/archive/blacklab/SarahHannah/megalodon/outputs/crc557/intercept_557.bed")
colnames(intercept_557) <- c("chr", "start", "end", "read_depth", "strand", "coverage", "cpg")
intercept_557$coverage <- intercept_557$coverage/100
intercept_557<-merge(intercept_557, crc557_clin_B_values, by= "cpg")
head(intercept_557)


full_meth_557 <- ggplot(intercept_557, mapping= aes(x=coverage, y=beta, colour= read_depth)) + geom_point() + ggtitle("Methylation: array vs nanopore for crc557") + labs(x="Methylation at each CpG site: Nanopore", y= "Methylation at each CpG site: Array")
print(full_meth_557 + scale_color_gradient(low= "light blue", high= "dark blue")) + geom_cor()
```

Compare samples for CRC557 with read depth >=8
```{r}
# subset
cpg_557 <- subset(intercept_557, read_depth>=8)
head(cpg_557)

#plot Nanopore vs arrays

ggplot(cpg_557, mapping = aes(x= coverage, y= beta)) + geom_point() + labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc557 with read depth >=8: nanopore vs array") + geom_cor()

RD_8_plot <-ggplot(cpg_557, mapping = aes(x= coverage, y= beta, colour=read_depth)) + geom_point() + labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc557 with read depth >=8: nanopore vs array")
print(RD_8_plot + scale_colour_gradient(low= "light pink", high= "dark red")) + geom_cor()

ggplot(cpg_557, mapping = aes(x= coverage, y= beta)) + geom_violin() + geom_point(alpha =0.3) +  labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc557 with read depth >=8: nanopore vs array")

#ggplot(cpg_557, mapping = aes(x= factor(coverage), y= beta)) +  geom_violin() + geom_point(size= 0.3) +
  #labs(x= "methylated regions for crc557 with read depth >=8: nanopore vs array") + ggtitle("coverage of methylated regions for sample crc557 comparing nanopore vs array")

ggplot(cpg_557, mapping = aes(x= coverage, y= beta)) + geom_jitter(width = 0.01, alpha=0.1) +  labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc557 with read depth >=8: nanopore vs array")

```

Subset for readdepth >=12

```{r}
cov_12_557 <- subset(cpg_557, read_depth >=12)
head(cov_12_557)
ggplot(cov_12_557, mapping = aes(x= coverage, y= beta)) + geom_point() + labs(x= "Nanopore methylation coverage", y= "array based methylation coverage") + ggtitle("methylated regions for crc636 with read depth >=12: nanopore vs array") + geom_cor ()

```

