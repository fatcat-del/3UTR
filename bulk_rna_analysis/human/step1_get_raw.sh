echo get rawdata
echo $1
cd /data/fengw/utr/bulk_raw/human
/data/fengw/software/sratoolkit.3.0.10-ubuntu64/bin/prefetch $1
cd $1
/data/fengw/software/sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-3 ${1}.sra

r1=${1}_1.fastq
r2=${1}_2.fastq
c1=${1}_R1.clean.fq.gz
c2=${1}_R2.clean.fq.gz
sam=${1}.sam

echo QC processing 
fastp -i $r1 -o $c1 -I $r2 -O $c2

echo Mapping processing
hisat2 -p 12 -t -x /data/fengw/utr/bulk_raw/human/humanDB/hg38_index -1 $c1 -2 $c2 -S $sam

echo sort bam
samtools view -@ 24 -b $sam > ${1}.bam
samtools sort -@ 24 ${1}.bam -o ${1}.sorted.bam
samtools index -@ 24 ${1}.sorted.bam

rm -rf ${1}.sra $r1 $r2 $c1 $c2 $sam
