path=/data/fengw/utr/bulk_raw/
for i in `less file.list`;do
	echo $i
	cd $path/$i
	/data/fengw/software/sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-3 $i
	cd $path
done
