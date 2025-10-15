sample=Dicty_D0_P3_TTGTGTAA

wd=${sample}/split_exon_intron

mkdir $wd

samtools view -b -q 250 ${sample}/Aligned.r1.bam | samtools sort - > ${wd}/r1.bam
samtools index ${wd}/r1.bam

cutadapt -g AGGAGAGTGTGAGTGATGGTTGAGGATGTGTGGAG -o ${sample}/${sample}_trimmed_2.fq.gz ${sample}/${sample}_2.fq.gz > ${sample}/read2_trimming_report.txt

htseq-count -s no -t exon -m intersection-strict -f bam ${wd}/r1.bam genes.gtf -q -o ${wd}/htseq_mappedr1.exon.sam > ${wd}/htseqr1.exon.out

python snapTotal-seq/scripts/readtaggene_Htseqcount_py38.py ${wd}/htseq_mappedr1.exon.sam > ${wd}/r1tagmapped.exon.dat

seqtk subseq ${sample}/${sample}_trimmed_2.fq.gz ${wd}/r1tagmapped.exon.dat > ${wd}/${sample}-mapped-r2.exon.fastq

python snapTotal-seq/scripts/barcodeseq_extract_standard_run.py ${wd}/${sample}-mapped-r2.exon.fastq ${wd}/r1tagmapped.exon.dat 108 > ${wd}/${sample}_mapped_barcode.exon.dat

htseq-count -s no -t mRNA -m intersection-strict -f bam ${wd}/r1.bam genes.gtf -q -o ${wd}/htseq_mappedr1.transcript.sam > ${wd}/htseqr1.transcript.out

python snapTotal-seq/scripts/readtaggene_Htseqcount_py38.py ${wd}/htseq_mappedr1.transcript.sam > ${wd}/r1tagmapped.transcript.dat

seqtk subseq ${sample}/${sample}_trimmed_2.fq.gz ${wd}/r1tagmapped.transcript.dat > ${wd}/${sample}-mapped-r2.transcript.fastq 

python snapTotal-seq/scripts/barcodeseq_extract_standard_run.py ${wd}/${sample}-mapped-r2.transcript.fastq ${wd}/r1tagmapped.transcript.dat 108 > ${wd}/${sample}_mapped_barcode.transcript.dat

python snapTotal-seq/scripts/count_exon_intron_UMI.py ${wd}/${sample}

gzip ${wd}/${sample}-mapped-r2.exon.fastq
gzip ${wd}/${sample}-mapped-r2.transcript.fastq


