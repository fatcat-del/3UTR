#/bin/bash

conda activate utr

gtf=/data/fengw/software/refdata-gex-mm10-2020-A/genes/genes.gtf
hfatbam=/data/raw/singlecell/liversc/force/HFAT/outs/possorted_genome_bam.bam
controlbam=/data/raw/singlecell/liversc/force/control/outs/possorted_genome_bam.bam

#Extract the UTR regions from the gtf file
cd /data/fengw/utr
python /data/fengw/software/SCAPE/scripts/main.py \
--gtf $gtf
--prefix mm10

#Process 10X bam files to be compatible with the input of the SCAPE
cd /data/fengw/utr/control
samtools view -H $controlbam > header.sam
sed '/^@SQ/{s/SN:chr\([0-9XY]\)/SN:\1/g;s/SN:chrM/SN:MT/g}' header.sam > newhead
samtools reheader newhead /data/raw/singlecell/liversc/force/control/outs/possorted_genome_bam.bam > output.bam
samtools sort -@ 8 output.bam -o out.sorted.bam
samtools index -@ 8 out.sorted.bam

#ident control single cell apa events
python /data/fengw/software/SCAPE/scripts/main.py apamix \
--bed mm10_utr.bed \
--bam out.sorted.bam\
--out control/ \
--cores 12 \
--cb barcodes.tsv.gz

#Process 10X bam files to be compatible with the input of the SCAPE
cd /data/fengw/utr/HFAT
samtools view -H $hfatbam > header.sam
sed '/^@SQ/{s/SN:chr\([0-9XY]\)/SN:\1/g;s/SN:chrM/SN:MT/g}' header.sam > newhead
samtools reheader newhead /data/raw/singlecell/liversc/force/control/outs/possorted_genome_bam.bam > output.bam
samtools sort -@ 8 output.bam -o out.sorted.bam
samtools index -@ 8 out.sorted.bam

#ident HFAT single cell apa events
python /data/fengw/software/SCAPE/scripts/main.py apamix \
--bed mm10_utr.bed \
--bam out.sorted.bam \
--out control/ \
--cores 12 \
--cb barcodes.tsv.gz

#merge control and HFAT apa events
cd /data/fengw/utr
python /data/fengw/software/SCAPE/scripts/group_pa.py \
--files /data/fengw/utr/HFAT/pasite.csv.gz,/data/fengw/utr/control/control/pasite.csv.gz \
--labels HFAT,control \
--outfile collapse_pa.tsv.gz



