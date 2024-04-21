for i in /storage1/fs1/leyao.wang/Active/RNA_SEQ_DATA_2024/RAW_FASTQ_FILES/combined/*R1_001.fastq; do name=$(basename ${i} R1_001.fastq);
       	echo ${name}; 
	STAR --genomeDir /storage1/fs1/leyao.wang/Active/Users/run/rna.attempts/genomeDir --runThreadN 12 --readFilesIn /storage1/fs1/leyao.wang/Active/RNA_SEQ_DATA_2024/RAW_FASTQ_FILES/combined/${name}R1_001.fastq  /storage1/fs1/leyao.wang/Active/RNA_SEQ_DATA_2024/RAW_FASTQ_FILES/combined/${name}R2_001.fastq --outFileNamePrefix /storage1/fs1/leyao.wang/Active/RNA_SEQ_DATA_2024/ANALYSIS/STAR/${name} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterMatchNmin  40; done


