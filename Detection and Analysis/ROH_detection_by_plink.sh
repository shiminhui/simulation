#!/bin/bash

#step1: VCF to PLINK format
vcftools --gzvcf Tiger.SNP_true.DP0.5.PL20.maf01.recode.vcf.gz --plink --out Tiger_Leopard
plink --noweb --file Tiger_Leopard --make-bed --out Tiger_Leopard

#step2: LD pruning
plink --bfile Tiger_Leopard --indep-pairwise 50 1 0.8
plink --bfile Tiger_Leopard --extract plink.prune.in --make-bed --out Tiger_Leopard_pruned
plink --bfile Tiger_Leopard_pruned --pca 10 --allow-extra-chr --threads 6

#step3: ROH detection
plink --bfile Tiger_Leopard_pruned --homozyg --homozyg-window-snp 80 --homozyg-density 50 --homozyg-kb 10
