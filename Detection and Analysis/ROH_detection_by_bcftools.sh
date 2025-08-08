#!/bin/bash

#step1: Use the same vcf.gz file
bcftools view ../roh_plink/Tiger.SNP_true.DP0.5.PL20.maf01.recode.vcf.gz -Ob -o Tiger.bcf
bcftools index Tiger.bcf

#step2: ROH detection, the parameters are as consistent as possible with the plink method
bcftools roh --AF-dflt 0.01 --GTs-only 20 --rec-rate 1e-8 --hw-to-az 0.05 --output-type r --threads 6 -o roh_output.txt Tiger.bcf

#step3: Screening ROH fragment length
awk ‘BEGIN{OFS=“\t”} $1==“RG” && $3-$2 >= 10000’ roh_output.txt > final_roh.bed
