#!/bin/sh
# Name : Sara Nicholson
# Date : February 2024
# Description: Run RSeQC after aligning reads to Genome in order to remove ribosomal genes. Reads not aligned to Ribosomal genes will be filtered into *.ex.bam. Must have samtools & RSeQc installed. 

for i in /leyao/Active/RNA_SEQ_DATA_2024/ANALYSIS/STAR/*sortedByCoord.out.bam; do name=$(basename ${i} .sortedByCoord.out.bam);
       	echo ${name};
	# Index Bam file via samtools (installed in path)
	samtools index /leyao/Active/RNA_SEQ_DATA_2024/ANALYSIS/STAR/${name}.sortedByCoord.out.bam
	# Remove Ribosomal Gene Reads (desired output will be *.ex.bam)
	python3 /home/sara/Downloads/RSeQC-5.0.1/scripts/split_bam.py -i /leyao/Active/RNA_SEQ_DATA_2024/ANALYSIS/STAR/${name}.sortedByCoord.out.bam -r /home/sara/Downloads/RSeQC-5.0.1/scripts/hg38_rRNA.bed -o /leyao/Active/RNA_SEQ_DATA_2024/ANALYSIS/STAR/RSeQC_Out/${name}_ ; done

