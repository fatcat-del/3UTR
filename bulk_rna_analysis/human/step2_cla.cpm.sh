#!/bin/bash

cd /data/fengw/utr/bulk_raw/human
for i in `ls */*.sorted.bam`;do
	echo $i
	bamCoverage -b $i -o temp_output.bw --binSize 100 --normalizeUsing CPM --region chr1:211743457:211749898
	/data/fengw/software/bigWigToBedGraph temp_output.bw temp_output_bedGraph
	mkdir -p /data/fengw/utr/code/data/bed/human_bed/$i
	cp temp_output_bedGraph /data/fengw/utr/code/data/bed/human_bed/$i.bed
	rm temp_output.bw temp_output_bedGraph
done

