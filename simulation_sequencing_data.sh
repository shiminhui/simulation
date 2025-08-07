#!/bin/bash

#step1: Haploid genomes with different FROH as input templates for simulation of WGS data
perl simulation_haploidgenomes_with_SNP.2.pl -v Tiger.SNP_true.DP0.5.PL20.maf01.recode.vcf.gz -r refgenome.fa -s AT01 -o AT01 #vcf.gz from the variants calling step
sed -i 's/^>chr/>hap1chr/g' AT01_hap1.fa
sed -i 's/^>chr/>hap2chr/g' AT01_hap2.fa
cat AT01_hap1.fa AT01_hap2.fa >AT01.fa

#step2: Simulating sequencing reads, 30bp, 50bp, 75bp, 100bp, 150bp 
wgsim ../AT01.fa AT01.30bp.fq1 AT01.30bp.fq2 -N 1282340133 -1 30 -2 30 
wgsim ../AT01.fa AT01.50bp.fq1 AT01.50bp.fq2 -N 2137233556 -1 50 -2 50 
wgsim ../AT01.fa AT01.75bp.fq1 AT01.75bp.fq2 -N 1424822370 -1 75 -2 75 
wgsim ../AT01.fa AT01.100bp.fq1 AT01.100bp.fq2 -N 1068616778 -1 100 -2 100 
wgsim ../AT01.fa AT01.150bp.fq1 AT01.150bp.fq2 -N 712411185 -1 150 -2 150 

#step3: Get sequencing data at different depths, taking 100bp 20X as an example
seqtk sample -s 123 AT01.100bp.fq1 0.01 > AT01.100bp.20x.fq1
seqtk sample -s 123 AT01.100bp.fq2 0.01 > AT01.100bp.20x.fq2
gzip AT01.100bp.20x.fq1.gz
gzip AT01.100bp.20x.fq2.gz
