#!/bin/bash

datadir='/Users/lenatrnovec/scTotalRNA/'

lib_id='Dicty_D0_P3_TTTCGACA' # each sub-library is first demultiplexed based on i5 index using standard Illumina demultiplex pipeline

cd ${datadir}

read1='Dicty_D0_P3_TTTCGACA_1'
read2='Dicty_D0_P3_TTTCGACA_2'

gunzip -d ${read2}.fq.gz

python /Users/lenatrnovec/scTotalRNA/snapTotal-seq/scripts/Extract_cell_index.py ${read2}.fq ${lib_id}_read_index.txt 

gzip ${read2}.fq

cat /Users/lenatrnovec/scTotalRNA/snapTotal-seq/scripts/cell_index_list.txt | while read LINE
do
	grep ${LINE} ${lib_id}_read_index.txt > ${lib_id}_${LINE}.txt
	
	seqtk subseq ${read1}.fq.gz ${lib_id}_${LINE}.txt > ${lib_id}_${LINE}_1.fq
	gzip ${lib_id}_${LINE}_1.fq
	
	seqtk subseq ${read2}.fq.gz ${lib_id}_${LINE}.txt > ${lib_id}_${LINE}_2.fq
	gzip ${lib_id}_${LINE}_2.fq
done


