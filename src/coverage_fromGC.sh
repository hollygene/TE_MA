#!/bin/bash

#$ -q rcc-30d

# S.c. Gene Conversion
date

AR=(90 91 92)

for i in "${AR[@]}"
do
    /usr/local/bedtools/latest/bin/bamToBed -i Sample${i}.sorted.fixed.marked.realigned.bam > BedFiles/Sample${i}.bed
    cat CoverageHead.txt > Coverage/Sample${i}.Coverage10K.txt
    /usr/local/bedtools/latest/bin/coverageBed -a BedFiles/Sample${i}.bed -b ../../GeneConv/ReferenceGenome/Sc10K.bed|cut -f 1-2,4 | sort -k 1,1 -k2,2n >> Coverage/Sample${i}.Coverage10K.txt


done
