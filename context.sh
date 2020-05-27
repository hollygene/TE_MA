#PBS -S /bin/bash
#PBS -q batch
#PBS -N getfasta
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=12:00:00
#PBS -l mem=20gb
#PBS -M hcm14449@uga.edu
#PBS -m abe


# script for finding nucleotide context of SNPs

# using bedops to go from vcf to bed file
# then using getfasta to get nucleotide context

module load BEDOPS/2.4.30
# time vcf2bed < /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0_variants_final.vcf > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0_variants_final.bed

## use awk to subtract 1 from first position and add 1 to second position
# this will give you the trinucleotide context

awk '{$2-=1;$3+=1}1' OFS='\t' H0_variants_final.bed > H0_variants_final_tri.bed

## use bedtools getfasta to get trinucleotide context

module load BEDTools/2.29.2-GCC-8.2.0-2.31.1

bedtools getfasta -tab -fi /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/genome.337.fasta -bed /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0_variants_final_tri.bed -fo /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0_variants_context.txt
