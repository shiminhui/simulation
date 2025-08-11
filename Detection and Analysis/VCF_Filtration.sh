#!/bin/bash

#step1: CombineGVCFs
gatk --java-options "-Xmx30g" CombineGVCFs -R refgenome.fa \
  --variant /home/admin/genome/07.MP/02.bam/07.gvcf/AT01.g.vcf.gz \
  --variant /home/admin/genome/07.MP/02.bam/07.gvcf/AT02.g.vcf.gz \
  --variant /home/admin/genome/07.MP/02.bam/07.gvcf/AT06.g.vcf.gz \
  --variant /home/admin/genome/07.MP/02.bam/07.gvcf/AT07.g.vcf.gz \
  --variant /home/admin/genome/07.MP/02.bam/07.gvcf/AT08.g.vcf.gz \
  --variant /home/admin/genome/07.MP/02.bam/07.gvcf/AT10.g.vcf.gz \
  --variant /home/admin/genome/07.MP/02.bam/07.gvcf/AT11.g.vcf.gz \
  --variant /home/admin/genome/07.MP/02.bam/07.gvcf/AT12.g.vcf.gz \
  -O Tiger.gvcf.gz

#step2: GenotypeGVCFs
gatk --java-options "-Xmx30g" GenotypeGVCFs  -R refgenome.fa  --variant Tiger.gvcf.gz  -O Tiger.raw.vcf.gz

#step3: SelectVariants (SNP)
gatk --java-options "-Xmx15g" SelectVariants --select-type-to-include SNP -R refgenome.fa  -V Tiger.raw.vcf.gz -O Tiger.raw.snp.vcf.gz

#step4: Hard Filtering
gatk --java-options "-Xmx50g"  VariantFiltration -V  Tiger.raw.snp.vcf.gz --filter-expression "QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0  || SB > -1.0" --filter-name "snp_filter" -O Tiger.SNP_hardfil.vcf.gz
gatk  SelectVariants --exclude-filtered true  -R refgenome.fa  -V Tiger.SNP_hardfil.vcf.gz -O Tiger.SNP_hardfil_true.vcf.gz

#step5: Biallelic filtering
zcat Tiger.SNP_hardfil_true.vcf.gz |awk '$5!~","' | gzip > Tiger_SNP_hardfil_ture_2allelic.vcf.gz

#step6: Depth Filtering
bcftools query -f '%CHROM\t%POS\t%DP\n' Tiger_SNP_hardfil_ture_2allelic.vcf.gz > Tiger_2allelic.DP.txt
cat Tiger_2allelic.DP.txt |wc -l > Tiger_2allelic.DP.txt.line
gzip -c Tiger_2allelic.DP.txt > all.chr_snp.DP.txt.gz
perl /path/to/find_dp.pl all.chr_snp.DP.txt.gz all.DP0.25.txt <Number from Tiger_2allelic.DP.txt.line>
perl /path/to/filter_dp.pl Tiger_SNP_hardfil_ture_2allelic.vcf.gz  Tiger_SNP_hardfil_true.2allelic.dp2.5.vcf.gz <Sample size> <Number1 from all.DP0.25.txt>  <Number2 from all.DP0.25.txt>

#step7: PL value filtering
perl /path/to/filter_PLvalue.pl Tiger_SNP_hardfil_true.2allelic.dp2.5.vcf.gz Tiger.SNP_true.DP0.5.PL20.vcf.gz 0.2 20 8

#step8: MAF filtering
vcftools --gzvcf  Tiger.SNP_true.DP0.5.PL20.vcf.gz --maf 0.01 --recode --out Tiger.SNP_true.DP0.5.PL20.maf01 --recode-INFO-all

#step9: Rename
echo "NC_056660.1 chr1
NC_056666.1 chr2
NC_056672.1 chr3
NC_056674.1 chr4" >list
bcftools annotate --rename-chrs list Tiger.SNP_true.DP0.5.PL20.maf01.recode.vcf -Oz -o Tiger.SNP_true.DP0.5.PL20.maf01.recode.vcf.gz

#The final variant set can be used for ROH detection using two methods.
