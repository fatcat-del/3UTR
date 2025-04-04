for i in `less file.list`;do
	echo $i
	fastp -i /data/fengw/utr/bulk_raw/${i}/${i}_1.fastq -o /data/fengw/utr/bulk_raw/clean_data/${i}_R1.clean.fq.gz -I /data/fengw/utr/bulk_raw/${i}/${i}_2.fastq -O /data/fengw/utr/bulk_raw/clean_data/${i}_R2.clean.fq.gz
done
