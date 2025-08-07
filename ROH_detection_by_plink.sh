#!/bin/bash

#step1: alignment and variant calling
MegaBOLT --type full --runtype WGS --ref refgenome.fa --list fq.list --bwa 1 --no-bqsr 1 --hc4 1 --run-genotypegvcfs 1 --run-vqsr 0 --ERC GVCF --no-fastq-output 1 --no-bam-output-for-sort 1

#step2: CombineGVCFs
gatk --java-options "-Xmx30g" CombineGVCFs -R refgenome.fa \
 -variant AT01.bwa.sortdup.hc4.g.vcf.gz \
 -variant AT02.bwa.sortdup.hc4.g.vcf.gz \
 -variant AT06.bwa.sortdup.hc4.g.vcf.gz \
 -variant AT07.bwa.sortdup.hc4.g.vcf.gz \
 -variant AT08.bwa.sortdup.hc4.g.vcf.gz \
 -variant AT10.bwa.sortdup.hc4.g.vcf.gz \
 -variant AT11.bwa.sortdup.hc4.g.vcf.gz \
 -variant AT12.bwa.sortdup.hc4.g.vcf.gz \
 -O ATpopulation.g.vcf.gz 

#step3: GenotypeGVCFs
gatk --java-options "-Xmx100g  -XX:ParallelGCThreads=8 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true "  GenotypeGVCFs -R refgenome.fa -V ATpopulation.g.vcf.gz -O  ATpopulation.raw.vcf.gz

#step4: selectSNPs
gatk --java-options "-Xmx10g" SelectVariants --select-type-to-include SNP -R refgenome.fa -V ATpopulation.raw.vcf.gz -O all4.snp.vcf.gz

#step5: bgzip
gzip -cd all4.snp.vcf.gz > Tiger.SNP_hardfil_true.vcf
bgzip -c Tiger.SNP_hardfil_true.vcf > Tiger.SNP_hardfil_true.vcf.gz
tabix -p vcf Tiger.SNP_hardfil_true.vcf.gz






#step: VCF to PLINK格式
vcftools --gzvcf Tiger.SNP_true.DP0.5.PL20.maf01.recode.vcf.gz --plink --out Tiger_Leopard
plink --noweb --file Tiger_Leopard --make-bed --out Tiger_Leopard

#step: LD pruning
plink --bfile Tiger_Leopard --indep-pairwise 50 1 0.8
plink --bfile Tiger_Leopard --extract plink.prune.in --make-bed --out Tiger_Leopard_pruned
plink --bfile Tiger_Leopard_pruned --pca 10 --allow-extra-chr --threads 6

#step: detection
plink --bfile Tiger_Leopard_pruned --homozyg --homozyg-window-snp 80 --homozyg-density 50 --homozyg-kb 10
