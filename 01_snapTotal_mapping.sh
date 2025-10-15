#!/bin/bash
sample=Dicty_D0_P3_TTGTGTAA

# Optional: remove adapter sequence
# Recommended if the read length of read 1 exceeds the insert size of your library. 
cutadapt -a CTCCACACATCCTCAACCATCACTCAC -o ${sample}/${sample}_adapter_trimmed_1.fq.gz ${sample}/${sample}_1.fq.gz > ${sample}/adapter_trimming_report.txt

# Optional: remove potential UMI sequence
# Recommended if the read length of read 1 exceeds the insert size of your library. 
cutadapt -u -10 -m 30 -o ${sample}/${sample}_UMItrimmed_1.fq.gz ${sample}/${sample}_adapter_trimmed_1.fq.gz > ${sample}/UMI_trimming_report.txt
rm ${sample}/${sample}_adapter_trimmed_1.fq.gz



# generate genome

# docker run --rm -it --platform linux/amd64 -v /Users/lenatrnovec/scTotalRNA:/data -w /data quay.io/biocontainers/star:2.7.11b--h43eeafb_1 STAR   --runMode genomeGenerate      --runThreadN 1   --genomeDir Dicty_STAR_index   --genomeFastaFiles genome.fa      --genomeSAindexNbases 12   --sjdbGTFfile genes.gtf   --sjdbGTFfeatureExon exon   --sjdbOverhang 100

# mapping read1

docker run --rm -it --platform linux/amd64 -v /Users/lenatrnovec/scTotalRNA:/data -w /data quay.io/biocontainers/star:2.7.11b--h43eeafb_1 STAR --runThreadN 8 --runMode alignReads --genomeDir Dicty_STAR_index --outFilterMismatchNmax 5 --readFilesCommand zcat --readFilesIn ${sample}/${sample}_UMItrimmed_1.fq.gz --outSAMtype BAM Unsorted --outFileNamePrefix ${sample}/

####
# Use ${sample}_1.fq.gz for mapping if no trimming was performed
####

samtools sort ${sample}/Aligned.out.bam -o  ${sample}/Aligned.r1.bam
rm -rf ${sample}/Aligned.out.bam


