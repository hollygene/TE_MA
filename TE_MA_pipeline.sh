#PBS -S /bin/bash
#PBS -q highmem_1
#PBS -N test_Spike_InsJune2019
#PBS -l nodes=1:ppn=10:Intel
#PBS -lwalltime=480:00:00
#PBS -l  mem=200gb
#PBS -M hcm14449@uga.edu
#PBS -m abe

#S paradoxus TE MA Quality Control pipepline

#location of current update of muver
muver_module="muver/0.1.0-foss-2016b-Python-2.7.14-20190318"
#location of trimgalore moedule
trimgalore_module="Trim_Galore/0.4.5-foss-2016b"
#location of fastqc module
fastqc_module="FastQC/0.11.8-Java-1.8.0_144"
#location of BWA module
bwa_module="BWA/0.7.17-foss-2016b"
#location of samtools module
samtools_module="SAMtools/1.9-foss-2016"
#location of bedtools module
bedtools_module="BEDTools/2.28.0-foss-2018a"
#location of python module
python_module="Python/3.5.2-foss-2016b"
#location of bamtoBigWig script and accessories
script_location="/scratch/hcm14449/TE_MA_Paradoxus/jbscripts"
#location of bam to bigwig script
bamToBigWig="/scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py"
#location of data to be analyzed
data_dir="/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein"
anc_dir="/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/"
#location of reference genome to be used
ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/S288C_reference_sequence_R64-2-1_20150113.fa"
#directory reference genome is located in
ref_genome_dir="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/"
#text file listing the fastq files with their full extensions
fastq_list="/home/hcm14449/Github/TE_MA/FASTQ_LIST.txt"
#what sample should all other samples be compared to?
control_sample_name="Ancestor"
#where should the output be sent
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/QCed_spike_ins"
mkdir $experiment_directory
#location of TRIMMED data to be used in the analysis
trimmed_data="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/samples/trimmed"


########################################################################################################################
cd ${data_dir}

mkdir ${output_directory}

module load ${fastq_module}

for file in ${data_dir}/*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}

time fastqc -o ${output_directory} ${data_dir}/${BASE}.fastq

done


#######################################################################################

cd ${data_dir}

mkdir ${trimmed_data}

module load ${trimgalore_module}

# trim all fastq files
for file in $data_dir/*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}

trim_galore --phred33 -q 20 -o $trimmed_data ${BASE}.fastq

done

module unload ${trimgalore_module}

 #######################################################################################


module load ${bwa_module}

 #index the ref genome
 bwa index ${ref_genome}

 for file in ${data_dir}*.fastq

 do

 FBASE=$(basename $file .fastq)
 BASE=${FBASE%.fastq}

 time bwa aln ${ref_genome} \
${data_dir}/${BASE}.fastq \
 > ${data_dir}/${BASE}.fastq.align.sai

 done


### for ancestors
for file in ${anc_dir}*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}

time bwa aln ${ref_genome} \
${anc_dir}/${BASE}.fastq \
> ${anc_dir}/${BASE}.fastq.align.sai

done

#######################################################################################

 ## have to do this separately for the arabidopsis samples bc different genome
# arabidopsis genome used :
# wget ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
#
# bwa index /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas

# for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/*.fastq
#
# do
#
# FBASE=$(basename $file .fastq)
# BASE=${FBASE%.fastq}
#
# time bwa aln /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas \
# /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.fastq \
# > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.fastq.align.sai
#
# done
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
 for file in ${data_dir}/*.fastq.align.sai

 do

 FBASE=$(basename $file .fastq.align.sai)
 BASE=${FBASE%.fastq.align.sai}

bwa samse ${ref_genome} \
    ${data_dir}/${BASE}.fastq.align.sai \
    ${data_dir}/${BASE}.fastq \
    > ${data_dir}/${BASE}.sam

done

### for ancestors

for file in ${anc_dir}/*.fastq.align.sai

do

FBASE=$(basename $file .fastq.align.sai)
BASE=${FBASE%.fastq.align.sai}

bwa samse ${ref_genome}\
   ${anc_dir}/${BASE}.fastq.align.sai \
   ${anc_dir}/${BASE}.fastq \
   > ${anc_dir}/${BASE}.sam

done

##### for arabidopsis samples
# for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/*.fastq.align.sai
#
# do
#
# FBASE=$(basename $file .fastq.align.sai)
# BASE=${FBASE%.fastq.align.sai}
#
# bwa samse /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas \
#    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.fastq.align.sai \
#    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.fastq \
#    > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.sam
#
# done

for file in ${data_dir}/*.fastq.align.sai

do

FBASE=$(basename $file .fastq.align.sai)
BASE=${FBASE%.fastq.align.sai}

bwa samse ${ref_genome} \
    ${data_dir}/${BASE}.fastq.align.sai \
    ${data_dir}/${BASE}.fastq \
    > ${data_dir}/${BASE}.sam

done

# bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
#         /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Col_500_S18_R1_001.fastq.align.sai \
#         /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Col_500_S18_R1_001.fastq \
#         > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Col_500_S18_R1_001.sam
#
# bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
#                 /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_12_S14_R1_001.fastq.align.sai \
#                 /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_12_S14_R1_001.fastq \
#                 > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_12_S14_R1_001.sam
#
# bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
#               /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_27_New_S21_R1_001.fastq.align.sai \
#               /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_27_New_S21_R1_001.fastq \
#               > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_27_New_S21_R1_001.sam
#
# bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
#                             /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_44_S15_R1_001.fastq.align.sai \
#                             /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_44_S15_R1_001.fastq \
#                             > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_44_S15_R1_001.sam
#
# bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
#     /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_5_S13_R1_001.fastq.align.sai \
#     /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_5_S13_R1_001.fastq \
#     > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D0_5_S13_R1_001.sam
#
#     bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
#         /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D1_11_New_S20_R1_001.fastq.align.sai \
#         /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D1_11_New_S20_R1_001.fastq \
#         > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D1_11_New_S20_R1_001.sam
#
#   bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
#     /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_1_S16_R1_001.fastq.align.sai \
#     /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_1_S16_R1_001.fastq \
#     > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_1_S16_R1_001.sam
#
#     bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
      # /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_2_New_S19_R1_001.fastq.align.sai \
      # /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_2_New_S19_R1_001.fastq \
      # > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_2_New_S19_R1_001.sam
      #
      # bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
      #   /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_8_S17_R1_001.fastq.align.sai \
      #   /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_8_S17_R1_001.fastq \
      #   > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/D20_8_S17_R1_001.sam
      #
      #   bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
      #     /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_33_New_S22_R1_001.fastq.align.sai \
      #     /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_33_New_S22_R1_001.fastq \
      #     > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_33_New_S22_R1_001.sam
      #
      #
      #             bwa samse /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
      #               /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_37_New_S23_R1_001.fastq.align.sai \
      #               /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_37_New_S23_R1_001.fastq \
      #               > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/H0_37_New_S23_R1_001.sam
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
module load ${samtools_module}

# samtools index *.bam
#index reference genome

samtools faidx ${ref_genome}

#convert sam files to bam files
for file in ${data_dir}*.sam

do

FBASE=$(basename $file .sam)
BASE=${FBASE%.sam}

samtools view -bt ${ref_genome_dir}/*.fai \
${data_dir}/${BASE}.sam \
  > ${data_dir}/${BASE}.bam

done

### for ancestors
for file in ${anc_dir}/*.sam

do

FBASE=$(basename $file .sam)
BASE=${FBASE%.sam}

samtools view -bt ${ref_genome_dir}/*.fai \
${anc_dir}/${BASE}.sam \
  > ${anc_dir}/${BASE}.bam

done

### for arabidopsis samples
# samtools faidx /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas
#
# for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/*.sam
#
# do
#
# FBASE=$(basename $file .sam)
# BASE=${FBASE%.sam}
#
# samtools view -bt /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas.fai \
# /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.sam \
#   > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.bam
#
# done

### sort the bam files
for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}

samtools sort -o ${data_dir}/${BASE}.sorted.bam \
   ${data_dir}/${BASE}.bam

done


#### for ancestors
for file in ${anc_dir}/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}

samtools sort -o ${anc_dir}/${BASE}.sorted.bam \
   ${anc_dir}/${BASE}.bam

done
#### Arabidopsis
#
# for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/*.bam
#
# do
#
# FBASE=$(basename $file .bam)
# BASE=${FBASE%.bam}
#
# samtools sort -o /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.sorted.bam \
#    /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/${BASE}.bam
#
# done

#samtools view -bt paradoxus_gen.fna.fai /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D0.sam  > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D0.bam

#samtools view -bt paradoxus_gen.fna.fai /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D1.sam  > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D1.bam

#samtools view -bt paradoxus_gen.fna.fai /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20.sam  > /lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/trim/HM_D20.bam

#########################################################################################
# for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.sorted.bam
#
# do
#
# FBASE=$(basename $file .sorted.bam)
# BASE=${FBASE%.sorted.bam}
#
# bedtools bamtobed -i /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sorted.bam > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sorted.bed
#
# done
#
# ###########################################################################################
#
# for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/*.sorted.bed
#
# do
#
# FBASE=$(basename $file .sorted.bed)
# BASE=${FBASE%.sorted.bed}
#
# bedtools genomecov -i /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sorted.bed \
# -g /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa.fai \
#  > /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/${BASE}.sorted.bedg
#
# done

#################################################################################################################
#bam to BigWig
module load ${bedtools_module}
module load ${python_module}
module load ${samtools_module}
export PATH=${PATH}:${script_location}

## Loop
for file in ${data_dir}*.sorted.bam

do

FBASE=$(basename $file .sorted.bam)
BASE=${FBASE%.sorted.bam}

python3 ${bamToBigWig} -sort ${ref_genome_dir}/*.fai ${data_dir}/${BASE}.sorted.bam

done


### for ancestors
for file in ${anc_dir}/*.sorted.bam

do

FBASE=$(basename $file .sorted.bam)
BASE=${FBASE%.sorted.bam}

python3 ${bamToBigWig} -sort ${ref_genome_dir}/*.fai ${anc_dir}/${BASE}.sorted.bam

done
