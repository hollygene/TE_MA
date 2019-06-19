!/bin/bash
PBS -q batch
PBS -l nodes=1:ppn=1:AMD
PBS -l walltime=480:00:00
PBS -l mem=80gb

#S paradoxus TE MA pipepline

##ancestors: spike-ins


cd /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019

#module load java/jdk1.8.0_20 fastqc
#time fastqc -o QCed_ancestors_spike_ins HM_H0_S16_R1_001.fastq.gz HM_D20_S15_R1_001.fastq.gz HMM_D1_S14_R1_001.fastq.gz HMM_D0_S13_R1_001.fastq.gz

#######################################################################################

cd /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019

#module load fastx/0.0.14
#time fastx_trimmer -f 10 \
#-i /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/HM_H0_S16_R1_001.fastq \
# -o /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_H0_trim.fastq

# time fastx_trimmer -f 10 \
# -i /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/HM_D20_S15_R1_001.fastq \
#  -o /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20_trim.fastq

#  time fastx_trimmer -f 10 \
#  -i /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/HMM_D1_S14_R1_001.fastq \
#   -o /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D1_trim.fastq

#   time fastx_trimmer -f 10 \
#   -i /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/HMM_D0_S13_R1_001.fastq \
#    -o /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D0_trim.fastq

 #######################################################################################


 module load bwa/0.7.15

 #index the ref genome
 bwa index /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa

 for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.fastq.gz

 do

 FBASE=$(basename $file .fastq)
 BASE=${FBASE%.fastq}

 time bwa aln /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.fastq \
 > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.fastq.align.sai

 done


 #time bwa aln /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
#/lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20_trim.fastq \
 #> /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HM_D20_trim.output.fastq.align.sai

 #time bwa aln /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
#/lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D1_trim.fastq \
# > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HMM_D1_trim.output.fastq.align.sai

 #time bwa aln /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
#/lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D0_trim.fastq \
# > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HMM_D0_trim.output.fastq.align.sai

 ###########################################################################################

 #create the sam files
 for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.fastq.align.sai

 do

 FBASE=$(basename $file .fastq.align.sai)
 BASE=${FBASE%.fastq.align.sai}

bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.fastq.align.sai \
    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.fastq \
    > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sam

done

#bwa samse /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
#          /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HMM_D0_trim.output.fastq.align.sai \
#          /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D0_trim.fastq \
#          > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D0.sam

#bwa samse /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
#                    /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HMM_D1_trim.output.fastq.align.sai \
#                    /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HMM_D1_trim.fastq \
#                    > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D1.sam

#bwa samse /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus_gen.fna \
#        /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/align/HM_D20_trim.output.fastq.align.sai \
#        /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20_trim.fastq \
#        > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20.sam

#########################################################################################

#samtools
module load samtools/1.6

#index reference genome

samtools faidx /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa

#convert sam files to bam files
for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.sam

do

FBASE=$(basename $file .sam)
BASE=${FBASE%.sam}

samtools view -bt /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa.fai \
/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sam \
  > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.bam

done

#samtools view -bt paradoxus_gen.fna.fai /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D0.sam  > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D0.bam

#samtools view -bt paradoxus_gen.fna.fai /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D1.sam  > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D1.bam

#samtools view -bt paradoxus_gen.fna.fai /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20.sam  > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20.bam

#########################################################################################
