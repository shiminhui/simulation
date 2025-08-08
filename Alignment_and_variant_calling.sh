#!/bin/bash

#step1: Sequence alignment
bwa mem -t 48 -k 30 -M -R "@RG\tID:AT01\tPL:ILLUMINA\tLB:AT01\tSM:AT01" refgenome.fa AT01.fq1.gz AT01.fq2.gz| samtools view -bS -o AT01.bam

#step2: Sort
samtools sort -@ 36 AT01.bam -o AT01.sort.bam

#step3: MarkDuplicates
java -jar picard.jar MarkDuplicates I=AT01.sort.bam O=AT01.dedup.bam M=AT01.dup_metrics.txt VALIDATION_STRINGENCY=SILENT

#step4: Build index
samtools index AT01.dedup.bam

#step5: Calculate the mapping rate
samtools flagstat AT01.dedup.bam > AT12.flagstat

#step6: Calculate sequencing depth
samtools depth -aa AT01.dedup.bam | awk '{chr[$1]+=$3; count[$1]++} END{print "Chromosome\tMeanDepth"; for(c in chr) print c "\t" chr[c]/count[c]}' > AT01.chr_mean_depth.tsv

#step7: HaplotypeCaller
java -Xmx8G -jar gatk-package-4.6.1.0-local.jar HaplotypeCaller -R refgenome.fa -I AT01.dedup.bam -O AT01.g.vcf.gz -ERC GVCF --native-pair-hmm-threads 8


#stepX: The above steps can be replaced by MegaBOLT
echo "AT01 AT01.fq1.gz AT01.fq2.gz" >fq.list
MegaBOLT --type full --runtype WGS -X-ref refgenome.fa --list fq.list --bwa 1 --no-bqsr 1 --hc4 1 --run-genotypegvcfs 1 --run-vqsr 0 --ERC GVCF --no-fastq-output 1 --no-bam-output-for-sort 1
