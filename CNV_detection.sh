#!/bin/bash

#step1
cnvnator -root AT01.root -tree AT01.marked.bam -chrom chr1 chr2 chr3 chr4

#step2
 cnvnator -root AT01.root -his 1000 -d dir/genome/

#step3
cnvnator -root AT01.root -stat 1000

#step4
cnvnator -root AT01.root -partition 1000

#step5
cnvnator -root AT01.root -call 1000 > AT01.root.call.txt
