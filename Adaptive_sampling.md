# Adpative Sampling
 
Sample CRC669 was used as an example, 
Each script is altered depending on the sample but base code for the command/process is the same
Run in bash unless stated otherwise

## Extracting Reads

Produce the text files and combined FASTQ dile used in the formation of the FASTQ containing the selected and unselected reads

```{bash, eval=FALSE}

#!/bin/bash

Input=/Volumes/archive/blacklab/SarahHannah/HCS-transfer/20210622_mblack/crc669/20210622_1456_X5_FAP36895_73fbf969/other_reports/adaptive_sampling_FAP36895_e2fefa30.csv
Output=/Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669
FASTQ= /Volumes/archive/blacklab/SarahHannah/fastq/crc669-fastq/pass/*.fastq.gz

#-d','  = use a comma as the delimiter to split the file into columns (fields)
#-f7    = select the 7th field (which contains the read IDs)

# Get a list of how many reads there were for each decision:
cat $Input | cut -d',' -f7 | sort | uniq -c

# Get the IDs of the selected reads:
cat $Input | grep stop_receiving | cut -d',' -f5 | sort | uniq > $Output/selected-ids.txt

# Get the IDs of the unselected (unblocked) reads:
cat $Input | grep unblock | cut -d',' -f5 | sort | uniq > $Output/unselected-ids.txt

#combine fastqs
cat $FASTQ | unpigz -p 48 > /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/big.fastq
```

Extract selected reads

```{Terminal, eval=FALSE}

export PATH=$PATH:/Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/seqtk

#choose only the selected reads
seqtk subseq /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/big.fastq selected-ids.txt > /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/selected-reads.fastq


#map the reads
./minimap2/minimap2 -t 24 -I8G -ax map-ont genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/selected-reads.fastq | samtools sort -@24 -o /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/selected_reads_669.bam

```
Extract Unselected reads

```{bash, eval=FALSE}

#choose only the unselected reads
seqtk subseq /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/big.fastq unselected-ids.txt > /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/unselected-reads.fastq


#map the reads
./minimap2/minimap2 -t 24 -I8G -ax map-ont genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/unselected-reads.fastq | samtools sort -@24 -o /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/unselected_reads_669.bam


```

## Nanoplot

selected reads
```{bash, eval=FALSE}
# needs to be run in a python environment

NanoPlot --bam selected_reads_669.bam -o /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/Nanoplot/selected/
```

unselected reads
```{bash, eval=FALSE}

# needs to be run in a python environment
#(run in the /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669 folder)

NanoPlot --bam unselected_reads_669.bam -o /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/Nanoplot/unselected/
```

## Generating counts
 
 selected reads
```{bash, eval=FALSE}
samtools index /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/selected_reads_669.bam

samtools bedcov /Volumes/archive/blacklab/SarahHannah/HCS-transfer/achen_black_ASfasta_bedfile/HS_GRCh38-coding_exons-promoters-4kflank.bed  /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/selected_reads_669.bam | gzip > /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/crc669-selected -read-depth-fast.txt.gz
```

unselected reads
```{bash, eval=FALSE}
samtools index /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/unselected_reads_669.bam

#complement bedfile to get the unselected reads since the bed file has the selected reads 
bedtools complement -i /Volumes/archive/blacklab/SarahHannah/HCS-transfer/achen_black_ASfasta_bedfile/HS_GRCh38-coding_exons-promoters-4kflank.bed -g genome/chrom.sizes > /Volumes/archive/blacklab/SarahHannah/Complement_bed/crc669_unselected.bed
 

samtools bedcov /Volumes/archive/blacklab/SarahHannah/Complement_bed/crc669_unselected.bed /Volumes/archive/blacklab/SarahHannah/output_minimap/crc669/crc669_sorted.bam | gzip > /Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/crc669_unselcted_all_read.txt.gz
```

## R-studio analysis

load packages

```{r}
library(dplyr)
library(ggplot2)
```
### CRC669

Crc669 Selected region
```{r}
Sel_crc669 <- read.table(gzfile("/Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/crc669-selected-read-depth-fast.txt.gz"))
Sel_crc669 <- subset(Sel_crc669, select = c("V1", "V2", "V3", "V7"))
colnames(Sel_crc669) <- c("Chr", "Start_position", "End_position", "Counts")
Sel_crc669$length <- Sel_crc669$End_position - Sel_crc669$Start_position 
Sel_crc669$avg_cov <- Sel_crc669$Counts / Sel_crc669$length
head (Sel_crc669) 
```
average coverage
```(r)
median(Sel_crc669$avg_cov)
mean(Sel_crc669$avg_cov)
```

selected regions coverage via length of region
```{r}
# separate into bins for coverage
breaks = c(0,250,500,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000, Inf) 
tags = c("0-250","250-500", "500-1000", "1000-1500", "1500-2000", "2000-3000", "3000-4000","4000-5000", "5000-6000", "6000-7000", "7000-8000", "8000-9000", "9000-10000", "10000+")
sel_len_bin_669<- cut(Sel_crc669$length, breaks = breaks, include.lowest = F, right = F, labels = tags)
summary(sel_len_bin_669)
Sel_crc669$bins <- factor(sel_len_bin_669, levels = tags, ordered=FALSE)
ggplot(Sel_crc669, mapping = aes(x=bins, y=avg_cov)) +geom_boxplot() +labs(title = "average coverage for length of regions for sample crc669 selected regions", x= "length of region (bp)", y="average coverage" )
```

Remove the extreme outliers as these are in repetitive regions and confound the data
set the cut off threshold at 40x coverage as avgerage coverage is only around 2

```{r}
Sel_crc669 <- subset(Sel_crc669, avg_cov <40)
ggplot(Sel_crc669, mapping = aes(x=bins, y=avg_cov)) +geom_boxplot() +labs(title = "average coverage for length of regions for sample crc669 selected regions", x= "length of region (bp)", y="average coverage" )
```


CRC669 unsl;ected regions
```{r}
Unsel_all_crc669 <- read.table(gzfile("/Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc669/crc669_unselcted_all_read.txt.gz"))
head(Unsel_all_crc669)
#Unsel_all_crc669 <- subset(Unsel_crc669, select= c("V1", "V2", "V3", "V4"))
colnames(Unsel_all_crc669) <- c("Chr", "Start_position", "End_position", "Counts")
#add length and average coverage
Unsel_all_crc669$length <- Unsel_all_crc669$End_position - Unsel_all_crc669$Start_position
Unsel_all_crc669$avg_cov <- Unsel_all_crc669$Counts / Unsel_all_crc669$length
head(Unsel_all_crc669)
median(Unsel_all_crc669$avg_cov)
mean(Unsel_all_crc669$avg_cov)
```

Unselected regions coverage via length of region

```{r}
breaks = c(0,250,500,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000, Inf) 
tags = c("0-250","250-500", "500-1000", "1000-1500", "1500-2000", "2000-3000", "3000-4000","4000-5000", "5000-6000", "6000-7000", "7000-8000", "8000-9000", "9000-10000", "10000+")
unsel_all_len_bin_669<- cut(Unsel_all_crc669$length, breaks = breaks, include.lowest = F, right = F, labels = tags)
summary(unsel_all_len_bin_669)
Unsel_all_crc669$bins <- factor(unsel_all_len_bin_669, levels = tags, ordered=FALSE)
ggplot(Unsel_all_crc669, mapping = aes(x=bins, y=avg_cov)) +geom_boxplot() +labs(title = "average coverage for length of regions for sample crc669 unselected regions", x= "length of region (bp)", y="average coverage" )

```

Remove outliers at coverage of 30x

```{r}
Unsel_all_crc669<- subset(Unsel_all_crc669, avg_cov <30)
ggplot(Unsel_all_crc669, mapping = aes(x=bins, y=avg_cov)) +geom_boxplot() +labs(title = "average coverage for length of regions for sample crc669 unselected regions", x= "length of region (bp)", y="average coverage" )
```

### CRC636


CRC636 Selected regions
```{r}
Sel_crc636 <- read.table(gzfile("/Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc636/crc636-selected-read-depth-fast.txt.gz"))
Sel_crc636 <- subset(Sel_crc636, select = c("V1", "V2", "V3", "V7"))
colnames(Sel_crc636) <- c("Chr", "Start_position", "End_position", "Counts")

# add length and average coverage
Sel_crc636$length <- Sel_crc636$End_position - Sel_crc636$Start_position 
Sel_crc636$avg_cov <- Sel_crc636$Counts / Sel_crc636$length
head (Sel_crc636) 
```

average coverage
```{r}
median(Sel_crc636$avg_cov)
mean(Sel_crc636$avg_cov)
```

sort coverage based on read length
```{r}
breaks = c(0,250,500,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000, Inf) 
tags = c("0-250","250-500", "500-1000", "1000-1500", "1500-2000", "2000-3000", "3000-4000","4000-5000", "5000-6000", "6000-7000", "7000-8000", "8000-9000", "9000-10000", "10000+")
sel_len_bin_636 <- cut(Sel_crc636$length, breaks = breaks, include.lowest = T, right = F, labels = tags)
summary(sel_len_bin_636)
Sel_crc636$bins <- factor(sel_len_bin_636, levels = tags, ordered=FALSE)
ggplot(Sel_crc636, mapping = aes(x=bins, y=avg_cov)) + geom_boxplot() + labs(title = "average coverage for length of selected regions for sample crc636", x= "length of region (bp)", y="average coverage" )
```

remove outliers, set outlier threshold at 50x coverage
```{r}
Sel_crc636 <- subset(Sel_crc636, avg_cov <50)

ggplot(Sel_crc636, mapping = aes(x=bins, y=avg_cov)) + geom_boxplot() + labs(title = "average coverage for length of selected regions for sample crc636", x= "length of region (bp)", y="average coverage" )
```

Unselected regions (all read bed)

```{r}
Unsel_all_crc636 <- read.table(gzfile("/Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc636/crc636_unselcted_all_read.txt.gz"))
head(Unsel_all_crc636)
#Unsel_all_crc669 <- subset(Unsel_crc669, select= c("V1", "V2", "V3", "V4"))
colnames(Unsel_all_crc636) <- c("Chr", "Start_position", "End_position", "Counts")
#add length and average coverage
Unsel_all_crc636$length <- Unsel_all_crc636$End_position - Unsel_all_crc636$Start_position
Unsel_all_crc636$avg_cov <- Unsel_all_crc636$Counts / Unsel_all_crc636$length
head(Unsel_all_crc636)
median(Unsel_all_crc636$avg_cov)
mean(Unsel_all_crc636$avg_cov)
```

unselected regions coverage via length of region

```{r}
breaks = c(0,250,500,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000, Inf) 
tags = c("0-250","250-500", "500-1000", "1000-1500", "1500-2000", "2000-3000", "3000-4000","4000-5000", "5000-6000", "6000-7000", "7000-8000", "8000-9000", "9000-10000", "10000+")
unsel_all_len_bin_636<- cut(Unsel_all_crc636$length, breaks = breaks, include.lowest = F, right = F, labels = tags)
summary(unsel_all_len_bin_636)
Unsel_all_crc636$bins <- factor(unsel_all_len_bin_636, levels = tags, ordered=FALSE)
ggplot(Unsel_all_crc636, mapping = aes(x=bins, y=avg_cov)) +geom_boxplot() +labs(title = "average coverage for length of regions for sample crc669 unselected regions", x= "length of region (bp)", y="average coverage" )

```
remove outliers, threshold is 50x coverage
```{r}
Unsel_all_crc636 <- subset(Unsel_all_crc636, avg_cov <50)
ggplot(Unsel_all_crc636, mapping = aes(x=bins, y=avg_cov)) + geom_boxplot() + labs(title = "average coverage for length of unselected regions for sample crc636", x= "length of region (bp)", y="average coverage" )

```

### CRC557

CRC557 Selected regions
```{r}
Sel_crc557 <- read.table(gzfile("/Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc557/crc557-selected-read-depth-fast.txt.gz"))
Sel_crc557 <- subset(Sel_crc557, select = c("V1", "V2", "V3", "V7"))
colnames(Sel_crc557) <- c("Chr", "Start_position", "End_position", "Counts")

# add length and average coverage
Sel_crc557$length <- Sel_crc557$End_position - Sel_crc557$Start_position 
Sel_crc557$avg_cov <- Sel_crc557$Counts / Sel_crc557$length
head (Sel_crc557) 
```

average coverage
```{r}
median(Sel_crc557$avg_cov)
mean(Sel_crc557$avg_cov)
```

sort avgerage coverage based on read length
```{r}
breaks = c(0,250,500,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000, Inf) 
tags = c("0-250","250-500", "500-1000", "1000-1500", "1500-2000", "2000-3000", "3000-4000","4000-5000", "5000-6000", "6000-7000", "7000-8000", "8000-9000", "9000-10000", "10000+")
sel_len_bin_557 <- cut(Sel_crc557$length, breaks = breaks, include.lowest = F, right = F, labels = tags)
summary(sel_len_bin_557)
Sel_crc557$bins <- factor(sel_len_bin_557, levels = tags, ordered=FALSE)

#plot the average coverage against length groups  
ggplot(Sel_crc557, mapping = aes(x=bins, y=avg_cov)) + geom_boxplot() + labs(title = "average coverage for length of selected regions for sample crc557", x= "length of region (bp)", y="average coverage" )
```

remove outliers, threshold set at 50x coverage
```{r}
Sel_crc557 <- subset(Sel_crc557, avg_cov <50)
ggplot(Sel_crc557, mapping = aes(x=bins, y=avg_cov)) + geom_boxplot() + labs(title = "average coverage for length of selected regions for sample crc557", x= "length of region (bp)", y="average coverage" )
```


CRC557 Unselected regions 

```{r}
Unsel_all_crc557 <- read.table(gzfile("/Volumes/archive/blacklab/SarahHannah/adaptive_sampling/crc557/crc557_unselcted_all_read.txt.gz"))
head(Unsel_all_crc557)
#Unsel_all_crc669 <- subset(Unsel_crc669, select= c("V1", "V2", "V3", "V4"))
colnames(Unsel_all_crc557) <- c("Chr", "Start_position", "End_position", "Counts")

#add length and average coverage
Unsel_all_crc557$length <- Unsel_all_crc557$End_position - Unsel_all_crc557$Start_position
Unsel_all_crc557$avg_cov <- Unsel_all_crc557$Counts / Unsel_all_crc557$length
head(Unsel_all_crc557)
median(Unsel_all_crc557$avg_cov)
mean(Unsel_all_crc557$avg_cov)
```

unselected regions coverage via length of region

```{r}
breaks = c(0,250,500,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000, Inf) 
tags = c("0-250","250-500", "500-1000", "1000-1500", "1500-2000", "2000-3000", "3000-4000","4000-5000", "5000-6000", "6000-7000", "7000-8000", "8000-9000", "9000-10000", "10000+")
unsel_all_len_bin_557<- cut(Unsel_all_crc557$length, breaks = breaks, include.lowest = F, right = F, labels = tags)
summary(unsel_all_len_bin_557)
Unsel_all_crc557$bins <- factor(unsel_all_len_bin_557, levels = tags, ordered=FALSE)
ggplot(Unsel_all_crc557, mapping = aes(x=bins, y=avg_cov)) +geom_boxplot() +labs(title = "average coverage for length of regions for sample crc669 unselected regions", x= "length of region (bp)", y="average coverage" )
```

remove outliers, threshold set at 50x coverage
```{r}
Unsel_all_crc557 <- subset(Unsel_all_crc557, avg_cov <50)
head(Unsel_all_crc557)
ggplot(Unsel_all_crc557, mapping = aes(x=bins, y=avg_cov)) + geom_boxplot() + labs(title = "average coverage for length of unselected regions for sample crc557", x= "length of region (bp)", y="average coverage" )

```