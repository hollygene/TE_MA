#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N testScriptJuly19
#PBS -l nodes=1:ppn=32:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=512gb
#PBS -M hcm14449@uga.edu
#PBS -m abe

#S paradoxus TE MA Quality Control and Mutation Calling pipeline

#location of current update of muver
# muver_module="muver/0.1.0-foss-2016b-Python-2.7.14-20190318"
#location of trimgalore moedule
trimgalore_module="Trim_Galore/0.4.5-foss-2016b"
#location of fastqc module
# fastqc_module="FastQC/0.11.8-Java-1.8.0_144"
#location of BWA module
bwa_module="BWA/0.7.15-foss-2016b"
#location of samtools module
samtools_module="SAMtools/1.6-foss-2016b"
#location of bedtools module
bedtools_module="BEDTools/2.28.0-foss-2018a"
#location of python module
python_module="Python/3.5.2-foss-2016b"
#location of picard module
picard_module="picard/2.16.0-Java-1.8.0_144"
#location of GATK module
GATK_module="GATK/4.0.3.0-Java-1.8.0_144"
#location of bamtoBigWig script and accessories
script_location="/scratch/hcm14449/TE_MA_Paradoxus/jbscripts"
#location of bam to bigwig script
bamToBigWig="/scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py"
#location of data to be analyzed
data_dir="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq"
#location of reference genome to be used
ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa"
#directory reference genome is located in
ref_genome_dir="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus"
#text file listing the fastq files with their full extensions
# fastq_list="/home/hcm14449/Github/TE_MA/FASTQ_LIST.txt"
#what sample should all other samples be compared to?
# control_sample_name="Ancestor"
#where should the output be sent
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out"
# mkdir $output_directory
#location of TRIMMED data to be used in the analysis
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq"


# ########################################################################################################################
# # works: unzips files and runs FASTQC on all of them
# ########################################################################################################################
# # cd ${data_dir}
# #
# # mkdir ${output_directory}
# #
# # gunzip /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/*.gz
# #
# # module load ${fastqc_module}
# #
# # for file in ${data_dir}/*.fastq
# #
# # do
# #
# # FBASE=$(basename $file .fastq)
# # BASE=${FBASE%.fastq}
# #
# # time fastqc -o ${output_directory} ${data_dir}/${BASE}.fastq
# #
# # done
#
#
# #######################################################################################
# # works: trims files
# #######################################################################################
# cd ${data_dir}
#
# mkdir ${trimmed_data}
#
# module load ${trimgalore_module}
#
# #need to remove TruSeq adapters Index 6
#
#
#
# # trim all fastq files
# for file in $data_dir/*.fastq
#
# do
#
# FBASE=$(basename $file .fastq)
# BASE=${FBASE%.fastq}
#
# trim_galore --phred33 -q 20 -o $trimmed_data ${BASE}.fastq
#
# done
#
# module unload ${trimgalore_module}
#######################################################################################
# convert fastq to unmapped bam (ubam)
# add read group information
#######################################################################################

java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
    FASTQ=${raw_data}/${BASE}_R1_001.fastq \ #first read file of pair
    FASTQ2=${raw_data}/${BASE}_R2_001.fastq  \ #second read file of pair
    OUTPUT=${raw_data}/${BASE}_fastqtosam.bam \
    READ_GROUP_NAME=H0164.2 \ #required; changed from default of A
    SAMPLE_NAME=NA12878 \ #required
    LIBRARY_NAME=Solexa-272222 \ #required
    PLATFORM_UNIT=H0164ALXX140820.2 \
    PLATFORM=illumina \ #recommended
    SEQUENCING_CENTER=BI \
    RUN_DATE=2014-08-20T00:00:00-0400

 #######################################################################################
# works: aligns samples to reference genome. Output is a .sam file
#######################################################################################

module load ${bwa_module}
#
#  #index the ref genome
bwa index ${ref_genome}
#
for file in ${raw_data}/*_R1_001.fastq

 do

 FBASE=$(basename $file _R1_001.fastq)
 BASE=${FBASE%_R1_001.fastq}

bwa mem -M -t 12 ${ref_genome} ${raw_data}/${BASE}_R1_001.fastq ${raw_data}/${BASE}_R2_001.fastq > ${output_directory}/${BASE}_aln.sam

 done


#
# #########################################################################################
# #samtools: converts sam files to bam files and sorts them
# #########################################################################################
#
module load ${samtools_module}
#

#
# #convert sam files to bam files
# for file in ${output_directory}/D0/*_aln.sam
#
# do
#
# FBASE=$(basename $file _aln.sam)
# BASE=${FBASE%_aln.sam}
#
# samtools view -bt -@ 12 ${ref_genome_dir}/*.fai \
# ${output_directory}/D0/${BASE}_aln.sam \
#   > ${output_directory}/D0/${BASE}.bam
#
# done
#
# #convert sam files to bam files
# for file in ${output_directory}/D20/*_aln.sam
#
# do
#
# FBASE=$(basename $file _aln.sam)
# BASE=${FBASE%_aln.sam}
#
# samtools view -bt -@ 12 ${ref_genome_dir}/*.fai \
# ${output_directory}/D20/${BASE}_aln.sam \
#   > ${output_directory}/D20/${BASE}.bam
#
# done
#
# #convert sam files to bam files
# for file in ${output_directory}/H0/*_aln.sam
#
# do
#
# FBASE=$(basename $file _aln.sam)
# BASE=${FBASE%_aln.sam}
#
# samtools view -bt -@ 12 ${ref_genome_dir}/*.fai \
# ${output_directory}/H0/${BASE}_aln.sam \
#   > ${output_directory}/H0/${BASE}.bam
#
# done
#
# ############################
# ### sort the bam files
# ############################
#
for file in ${output_directory}/H0/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}

samtools sort -o -@ 12 ${output_directory}/H0/${BASE}.sorted.bam \
   ${output_directory}/H0/${BASE}.bam

done
#
# for file in ${output_directory}/D0/*.bam
#
# do
#
# FBASE=$(basename $file .bam)
# BASE=${FBASE%.bam}
#
# samtools sort -o -@ 12 ${output_directory}/D0/${BASE}.sorted.bam \
#    ${output_directory}/D0/${BASE}.bam
#
# done
#
# for file in ${output_directory}/D1/*.bam
#
# do
#
# FBASE=$(basename $file .bam)
# BASE=${FBASE%.bam}
#
# samtools sort -o -@ 12 ${output_directory}/D1/${BASE}.sorted.bam \
#    ${output_directory}/D1/${BASE}.bam
#
# done
#
# for file in ${output_directory}/D20/*.bam
#
# do
#
# FBASE=$(basename $file .bam)
# BASE=${FBASE%.bam}
#
# samtools sort -o -@ 12 ${output_directory}/D20/${BASE}.sorted.bam \
#    ${output_directory}/D20/${BASE}.bam
#
# done
#
#
#
#
#
# ###################################################################################################
# ## Picard to mark duplicates
# ###################################################################################################
#
module load ${picard_module}

for file in ${output_directory}/H0/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}

time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
I=${output_directory}/H0/${BASE}.bam \
O=${output_directory}/H0/${BASE}_removedDuplicates.bam \
M=${output_directory}/H0/${BASE}_removedDupsMetrics.txt

done

#
# for file in ${output_directory}/D0/*.bam
#
# do
#
# FBASE=$(basename $file .bam)
# BASE=${FBASE%.bam}
#
# time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/D0/${BASE}.bam \
# O=${output_directory}/D0/${BASE}_removedDuplicates.bam \
# M=${output_directory}/D0/${BASE}_removedDupsMetrics.txt
#
# done


# for file in ${output_directory}/D1/*.bam
#
# do
#
# FBASE=$(basename $file .bam)
# BASE=${FBASE%.bam}
#
# time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/D1/${BASE}.bam \
# O=${output_directory}/D1/${BASE}_removedDuplicates.bam \
# M=${output_directory}/D1/${BASE}_removedDupsMetrics.txt
#
# done


# for file in ${output_directory}/D20/*.bam
#
# do
#
# FBASE=$(basename $file .bam)
# BASE=${FBASE%.bam}
#
# time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/D20/${BASE}.bam \
# O=${output_directory}/D20/${BASE}_removedDuplicates.bam \
# M=${output_directory}/D20/${BASE}_removedDupsMetrics.txt
#
# done

###################################################################################################
# Using GATK HaplotypeCaller in GVCF mode
# apply appropriate ploidy for each sample
# will need to do this separtely for haploid and diploid samples
###################################################################################################

module load ${GATK_module}

### H0 samples
for file in ${output_directory}/H0/*_removedDuplicates.bam

do

FBASE=$(basename $file *_removedDuplicates.bam)
BASE=${FBASE%*_removedDuplicates.bam}


time gatk HaplotypeCaller \
     -R ${ref_genome} \
     -ERC GVCF \
     -I ${output_directory}/H0/${BASE}_removedDuplicates.bam \
     -ploidy 1 \
     -O ${output_directory}/${BASE}_variants.g.vcf

done

#### D0 samples
for file in ${output_directory}/D0/*_removedDuplicates.bam

do

FBASE=$(basename $file *_removedDuplicates.bam)
BASE=${FBASE%*_removedDuplicates.bam}


time gatk HaplotypeCaller \
     -R ${ref_genome} \
     -ERC GVCF \
     -I ${output_directory}/${BASE}_removedDuplicates.bam \
     -ploidy 2 \
     -O ${output_directory}/${BASE}_variants.g.vcf

done

#### D1 samples
for file in ${output_directory}/D1/*_removedDuplicates.bam

do

FBASE=$(basename $file *_removedDuplicates.bam)
BASE=${FBASE%*_removedDuplicates.bam}


time gatk HaplotypeCaller \
     -R ${ref_genome} \
     -ERC GVCF \
     -I ${output_directory}/${BASE}_removedDuplicates.bam \
     -ploidy 2 \
     -O ${output_directory}/${BASE}_variants.g.vcf

done

#### D20 samples
for file in ${output_directory}/*_removedDuplicates.bam

do

FBASE=$(basename $file *_removedDuplicates.bam)
BASE=${FBASE%*_removedDuplicates.bam}


time gatk HaplotypeCaller \
     -R ${ref_genome} \
     -ERC GVCF \
     -I ${output_directory}/${BASE}_removedDuplicates.bam \
     -ploidy 2 \
     -O ${output_directory}/${BASE}_variants.g.vcf

done

###################################################################################################
### Aggregate the GVCF files using GenomicsDBImport
###################################################################################################

# gatk --java-options "-Xmx4g -Xms4g" \
#        GenomicsDBImport \
#        --genomicsdb-workspace-path /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB \
#        --batch-size 50 \
#        --sample-name-map /home/hcm14449/Github/TE_MA/H0_sample_map.txt \
#        --tmp-dir= /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB/tmp \
#        --reader-threads 12

###################################################################################################
### Alternatively, call SNPs with SAMtools
###################################################################################################

# samtools faidx ${reference_genome}
#
# samtools mpileup -g -f ${fasta} ${sorted_bam} > ${raw_bcf}
#
# bcftools view -bvcg ${raw_bcf} > ${var_bcf}
#
# #filter snps
# bcftools view ${var_bcf} | vcfutils.pl varFilter - > ${var_final.vcf}


#################################################################################################################
#bam to BigWig > for quality control purposes
#################################################################################################################
# module load ${bedtools_module}
# module load ${python_module}
# module load ${samtools_module}
# export PATH=${PATH}:${script_location}
#
## Loop
for file in ${output_directory}/H0/*.sorted.bam

do

FBASE=$(basename $file .sorted.bam)
BASE=${FBASE%.sorted.bam}

python3 ${bamToBigWig} -sort ${ref_genome_dir}/*.fai ${output_directory}/H0/${BASE}.sorted.bam

done



python3 /scratch/hcm14449/TE_MA_Paradoxus/bedGraphToBigWigScript -sort /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/*.fai /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D20/HM-D20-37_sorted.bam

#
# for file in ${output_directory}/D0/*.sorted.bam
#
# do
#
# FBASE=$(basename $file .sorted.bam)
# BASE=${FBASE%.sorted.bam}
#
# python3 ${bamToBigWig} -sort ${ref_genome_dir}/*.fai ${output_directory}/D0/${BASE}.sorted.bam
#
# done
#
# for file in ${output_directory}/D1/*.sorted.bam
#
# do
#
# FBASE=$(basename $file .sorted.bam)
# BASE=${FBASE%.sorted.bam}
#
# python3 ${bamToBigWig} -sort ${ref_genome_dir}/*.fai ${output_directory}/D1/${BASE}.sorted.bam
#
# done
#
for file in ${output_directory}/D20/*.sorted.bam

do

FBASE=$(basename $file .sorted.bam)
BASE=${FBASE%.sorted.bam}

python3 ${bamToBigWig} -sort ${ref_genome_dir}/*.fai ${output_directory}/D20/${BASE}.sorted.bam

done






for file in /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/BedGraphs/TELocate/D20/raw/*.bed

do

FBASE=$(basename $file .bed)
BASE=${FBASE%.bed}

awk 'NR==FNR{a[$1]; next} (($1) in a)' ${BASE}.bed D20-A-R_R1_val_telocate_raw.TY1.bed > ${BASE}_A_Common.txt

done
