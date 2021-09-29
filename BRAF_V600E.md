# Variant calling

## Extract BRAF V660E muation

### CRC669

```{bash, eval=F}
samtools view -b /Volumes/archive/blacklab/SarahHannah/output_minimap/crc669/crc669-sort.bam  "7:140753335-140753337" > /Volumes/archive/blacklab/SarahHannah/output_minimap/crc669/crc669-v600e.bam
samtools index /Volumes/archive/blacklab/SarahHannah/output_minimap/crc669/crc669-v600e.bam
```

### CRC636

```{bash, eval=F}
samtools view -b /Volumes/archive/blacklab/SarahHannah/output_minimap/crc636/crc636-sort.bam  "7:140753335-140753337" > /Volumes/archive/blacklab/SarahHannah/output_minimap/crc636/crc636-v600e.bam
samtools index /Volumes/archive/blacklab/SarahHannah/output_minimap/crc636/crc636-v600e.bam
```

### CRC557

```{bash, eval=F}
samtools view -b /Volumes/archive/blacklab/SarahHannah/output_minimap/crc557/crc557-sort.bam  "7:140753335-140753337" > /Volumes/archive/blacklab/SarahHannah/output_minimap/crc557/crc557-v600e.bam
samtools index /Volumes/archive/blacklab/SarahHannah/output_minimap/crc557/crc557-v600e.bam
```