#!/bin/bash
#PBS -N depth
#PBS -q highmem_q
#PBS -l nodes=1:ppn=1
#PBS -l walltime=480:00:00
#PBS -l mem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

H0_bams="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Muver/H0/bams"
D0_bams="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Muver/D0/bams"
D1_bams="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Muver/D1/bams"
D20_bams="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Muver/D20/bams"

module load SAMtools/1.9-foss-2016b

#H0
for file in ${H0_bams}/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}


samtools sort ${H0_bams}/${BASE}.bam \
-o ${H0_bams}/${BASE}.sorted.bam

samtools depth \
${H0_bams}/${BASE}.sorted.bam \
|  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${H0_bams}/${BASE}.txt

done


### D0
for file in ${D0_bams}/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}


samtools sort ${D0_bams}/${BASE}.bam \
-o ${D0_bams}/${BASE}.sorted.bam

samtools depth \
${D0_bams}/${BASE}.sorted.bam \
|  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${D0_bams}/${BASE}.txt

done

#D1
for file in ${D1_bams}/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}


samtools sort ${D1_bams}/${BASE}.bam \
-o ${D1_bams}/${BASE}.sorted.bam

samtools depth \
${D1_bams}/${BASE}.sorted.bam \
|  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${D1_bams}/${BASE}.txt

done

#D20
for file in ${D20_bams}/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}


samtools sort ${D20_bams}/${BASE}.bam \
-o ${D20_bams}/${BASE}.sorted.bam

samtools depth \
${D20_bams}/${BASE}.sorted.bam \
|  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${D20_bams}/${BASE}.txt

done
