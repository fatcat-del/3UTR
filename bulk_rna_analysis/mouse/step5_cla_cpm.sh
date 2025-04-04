#!/bin/bash

cd /data/fengw/utr/bulk_raw/20240223/bam_data/sortedbam
for i in `ls *.sorted.bam`;do
	echo $i
	bamCoverage -b $i -o temp_output.bw --binSize 100 --normalizeUsing CPM --region chr1:191778535:191784255
	/data/fengw/software/bigWigToBedGraph temp_output.bw temp_output_bedGraph
	cp temp_output_bedGraph /data/fengw/utr/code/data/bed/mouse_bed/$i.bed
	rm temp_output.bw temp_output_bedGraph
done



