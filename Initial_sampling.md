Intial generation

# Bedfiles

NB - the bed file was derived from information in this file:

HM450.hg38.manifest.tsv

from:

https://zwdzwd.github.io/InfiniumAnnotation

To generate the bed file the following command was used:

```{bash, eval=F}
cat HM450.hg38.manifest.tsv | cut -d$'\t' -f1-5 > HM450.hg38.manifest.bed
```

# Basecalling


Sample CRC669 was used as an example, 
Each script is altered depending on the sample but base code for the command/process is the same
Run in bash unless stated otherwise

Basecalled using Guppy produced by Oxford Nanopore Techonologies

```{bash, eval=FALSE} 

# PATH TO FAST5
FAST5=/Volumes/archive/blacklab/SarahHannah/HCS-transfer/20210622_mblack/crc669/20210622_1456_X5_FAP36895_73fbf969/fast5_pass

# DIRECTORY FOR OUTPUT
OUTPUT=/Volumes/archive/blacklab/SarahHannah/fastq/crc669-fastq

VERS=r9.4.1_450bps_$5
CFG=dna_${VERS}.cfg
MOD=template_${VERS}.jsn

# num_callers
CAL=$1
# gpu_runners_per_device
RUN=$2
# chunks_per_runner
CHK=$3
# chunk_size
CSZ=$4
DIR=./ont-guppy-gpu/data/

./ont-guppy-gpu/bin/guppy_basecaller --disable_pings --compress_fastq -i $FAST5 -s $OUTPUT -c ${DIR}${CFG} -m ${DIR}${MOD} -x 'cuda:0' --recursive --num_callers $CAL --gpu_runners_per_device $RUN --chunks_per_runner $CHK --chunk_size $CSZ
``` 

```{bash, eval=FALSE}
cd /Volumes/archive/blacklab/SarahHannah
#  - num_callers: 4
#  - gpu_runners_per_device: 8
#  - chunks_per_runner: 1024 
#  - chunk_size: 1000
#  - fast or hac (high accuracy): hac

./run-guppy-gpu.sh 4 8 1024 1000 hac
```

# Read Mapping

Reads were mapped using Minimap2 produced by Heng Li

```{bash, eval=FALSE}
#!/bin/bash

#ln -s /Volumes/archive/blacklab/MikBlack/Research/Nanopore/minimap2 
#ln -s /Volumes/archive/blacklab/MikBlack/Research/Nanopore/OGF/genome 

./minimap2/minimap2 -t 24 -I8G -ax map-ont genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa /Volumes/archive/blacklab/SarahHannah/fastq/crc669-fastq/pass/*.fastq.gz | samtools sort -@24 -o /Volumes/archive/blacklab/SarahHannah/output_minimap/crc669/crc669_sorted.bam

```
