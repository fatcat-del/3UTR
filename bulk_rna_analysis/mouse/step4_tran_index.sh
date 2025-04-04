for i in `less test`;do
	samtools view -@8 -b /data/fengw/utr/bulk_raw/bam_data/${i}.sam > /data/fengw/utr/bulk_raw/bam_data/${i}.bam
	echo $i
	samtools sort -@ 8 /data/fengw/utr/bulk_raw/bam_data/${i}.bam -o /data/fengw/utr/bulk_raw/bam_data/${i}.sorted.bam
	cd /data/fengw/utr/bulk_raw/bam_data/
	samtools index -@ 8 ${i}.sorted.bam
	cd /data/fengw/utr/bulk_raw/
done
