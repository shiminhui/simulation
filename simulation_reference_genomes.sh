#!/bin/bash

#The first refgenome.fa consists of four chromosomes from GCF_018350195.1 (A1: NC_056660.1, B4: NC_056666.1, D4: NC_056672.1, E2: NC_056674).
#Gaps are then inserted to form reference genomes with different contig N50.

perl simulation_refgenome_with_continuity.pl -i refgenome.fa -o refgenome.contig1.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig1.fa -o refgenome.contig2.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig2.fa -o refgenome.contig3.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig3.fa -o refgenome.contig4.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig4.fa -o refgenome.contig5.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig5.fa -o refgenome.contig6.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig6.fa -o refgenome.contig7.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig7.fa -o refgenome.contig8.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig8.fa -o refgenome.contig9.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig9.fa -o refgenome.contig10.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig10.fa -o refgenome.contig11.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig11.fa -o refgenome.contig12.fa
perl simulation_refgenome_with_continuity.pl -i refgenome.contig12.fa -o refgenome.contig13.fa

