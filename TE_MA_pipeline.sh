#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=80gb

#S paradoxus TE MA pipepline

##ancestors: spike-ins


cd /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA

module load java/jdk1.8.0_20 fastqc
time fastqc -o QCed_ancestors_spike_ins HM_H0_S16_R1_001.fastq.gz HM_D20_S15_R1_001.fastq.gz HMM_D1_S14_R1_001.fastq.gz HMM_D0_S13_R1_001.fastq.gz

#######################################################################################

cd /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA

module load fastx/0.0.14
time fastx_trimmer -f 10 \
-i /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/HM_H0_S16_R1_001.fastq \
 -o /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_H0_trim.fastq

 time fastx_trimmer -f 10 \
 -i /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/HM_D20_S15_R1_001.fastq \
  -o /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20_trim.fastq

  time fastx_trimmer -f 10 \
  -i /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/HMM_D1_S14_R1_001.fastq \
   -o /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D1_trim.fastq

   time fastx_trimmer -f 10 \
   -i /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/HMM_D0_S13_R1_001.fastq \
    -o /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D0_trim.fastq

 #######################################################################################


 module load bwa/0.7.15

 #index the ref genome
 bwa index /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna

 time bwa aln /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
/lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_H0_trim.fastq \
 > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HM_H0_trim.output.fastq.align.sai

 time bwa aln /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
/lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20_trim.fastq \
 > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HM_D20_trim.output.fastq.align.sai

 time bwa aln /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
/lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D1_trim.fastq \
 > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HMM_D1_trim.output.fastq.align.sai

 time bwa aln /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
/lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D0_trim.fastq \
 > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HMM_D0_trim.output.fastq.align.sai
