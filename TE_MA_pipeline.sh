#PBS -S /bin/bash
#PBS -q highmem_1
#PBS -N test_Spike_InsJune2019
#PBS -l nodes=1:ppn=10:Intel
#PBS -lwalltime=480:00:00
#PBS -l  mem=200gb
#PBS -M hcm14449@uga.edu
#PBS -m abe


#S paradoxus TE MA pipepline

##ancestors: spike-ins


cd /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019

mkdir /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/QCed_spike_ins

module load FastQC/0.11.8-Java-1.8.0_144

for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}

time fastqc -o /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/QCed_spike_ins /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.fastq

done


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


module load BWA/0.7.17-foss-2016b

 #index the ref genome
 bwa index /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa

 for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.fastq

 do

 FBASE=$(basename $file .fastq)
 BASE=${FBASE%.fastq}

 time bwa aln /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.fastq \
 > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.fastq.align.sai

 done


### for ancestors
for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}

time bwa aln /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/${BASE}.fastq \
> /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/${BASE}.fastq.align.sai

done

 ## have to do this separately for the arabidopsis samples bc different genome
# arabidopsis genome used :
wget ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

bwa index /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas

for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}

time bwa aln /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas \
/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.fastq \
> /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.fastq.align.sai

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

### for ancestors

for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/*.fastq.align.sai

do

FBASE=$(basename $file .fastq.align.sai)
BASE=${FBASE%.fastq.align.sai}

bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
   /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/${BASE}.fastq.align.sai \
   /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/${BASE}.fastq \
   > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/${BASE}.sam

done
##### for arabidopsis samples
for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/*.fastq.align.sai

do

FBASE=$(basename $file .fastq.align.sai)
BASE=${FBASE%.fastq.align.sai}

bwa samse /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas \
   /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.fastq.align.sai \
   /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.fastq \
   > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.sam

done

bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Col_200_S24_R1_001.fastq.align.sai \
    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Col_200_S24_R1_001.fastq \
    > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Col_200_S24_R1_001.sam

bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
        /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Col_500_S18_R1_001.fastq.align.sai \
        /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Col_500_S18_R1_001.fastq \
        > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Col_500_S18_R1_001.sam

bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
                /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_12_S14_R1_001.fastq.align.sai \
                /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_12_S14_R1_001.fastq \
                > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_12_S14_R1_001.sam

bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
              /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_27_New_S21_R1_001.fastq.align.sai \
              /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_27_New_S21_R1_001.fastq \
              > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_27_New_S21_R1_001.sam

bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
                            /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_44_S15_R1_001.fastq.align.sai \
                            /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_44_S15_R1_001.fastq \
                            > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_44_S15_R1_001.sam

bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_5_S13_R1_001.fastq.align.sai \
    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_5_S13_R1_001.fastq \
    > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_5_S13_R1_001.sam

    bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
        /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D1_11_New_S20_R1_001.fastq.align.sai \
        /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D1_11_New_S20_R1_001.fastq \
        > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D1_11_New_S20_R1_001.sam

  bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_1_S16_R1_001.fastq.align.sai \
    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_1_S16_R1_001.fastq \
    > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_1_S16_R1_001.sam

    bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
      /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_2_New_S19_R1_001.fastq.align.sai \
      /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_2_New_S19_R1_001.fastq \
      > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_2_New_S19_R1_001.sam

      bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
        /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_8_S17_R1_001.fastq.align.sai \
        /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_8_S17_R1_001.fastq \
        > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_8_S17_R1_001.sam

        bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
          /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_33_New_S22_R1_001.fastq.align.sai \
          /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_33_New_S22_R1_001.fastq \
          > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_33_New_S22_R1_001.sam


                  bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
                    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_37_New_S23_R1_001.fastq.align.sai \
                    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_37_New_S23_R1_001.fastq \
                    > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_37_New_S23_R1_001.sam
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
module load SAMtools/1.9-foss-2016b

samtools index *.bam
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

### for ancestors
for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/*.sam

do

FBASE=$(basename $file .sam)
BASE=${FBASE%.sam}

samtools view -bt /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa.fai \
/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/${BASE}.sam \
  > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/${BASE}.bam

done
### for arabidopsis samples
samtools faidx /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas

for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/*.sam

do

FBASE=$(basename $file .sam)
BASE=${FBASE%.sam}

samtools view -bt /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas.fai \
/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.sam \
  > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.bam

done

### sort the bam files
for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}

samtools sort -o /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sorted.bam \
   /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.bam

done


#### for ancestors
for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}

samtools sort -o /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/${BASE}.sorted.bam \
   /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/${BASE}.bam

done
#### Arabidopsis

for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}

samtools sort -o /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.sorted.bam \
   /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.bam

done




#samtools view -bt paradoxus_gen.fna.fai /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D0.sam  > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D0.bam

#samtools view -bt paradoxus_gen.fna.fai /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D1.sam  > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D1.bam

#samtools view -bt paradoxus_gen.fna.fai /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20.sam  > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20.bam

#########################################################################################


for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.sorted.bam

do

FBASE=$(basename $file .sorted.bam)
BASE=${FBASE%.sorted.bam}

bedtools bamtobed -i /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sorted.bam > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sorted.bed

done

###########################################################################################

for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.sorted.bed

do

FBASE=$(basename $file .sorted.bed)
BASE=${FBASE%.sorted.bed}

bedtools genomecov -i /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sorted.bed \
-g /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa.fai \
 > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sorted.bedg

done
