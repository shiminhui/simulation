==> step1.sh <==
cnvnator -root AT01.root -tree AT01.marked.bam -chrom NC_056660.1 NC_056666.1 NC_056672.1 NC_056674.1

==> step2.sh <==
 cnvnator -root AT01.root -his 1000 -d dir/genome/

==> step3.sh <==
cnvnator -root AT01.root -stat 1000

==> step4.sh <==
cnvnator -root AT01.root -partition 1000

==> step5.sh <==
cnvnator -root AT01.root -call 1000 > AT01.root.call.txt
