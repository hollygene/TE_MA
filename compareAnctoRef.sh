#!/bin/bash
#PBS -N CompareAncestorsToReference
#PBS -q highmem_q
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae


H0_bams="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Muver/H0/bams/"
ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa"
GATK_module="GATK/4.0.3.0-Java-1.8.0_144"

module load ${GATK_module}

time gatk HaplotypeCaller \
     -R ${ref_genome} \
     -I ${H0_bams}/H0-A.bam \
     -ploidy 1 \
     -O ${H0_bams}/H0-A_variants.vcf
