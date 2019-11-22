#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N H0_scripts
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l walltime=30:00:00
#PBS -l mem=300gb
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
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0"
# mkdir $output_directory
#location of TRIMMED data to be used in the analysis
trimmed_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/H0"
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/H0"
genomicsdb_workspace_path="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB"
sample_name_map="/home/hcm14449/Github/TE_MA/H0_sample_map.txt"
tmp_DIR="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB/tmp"

cd ${output_directory}
rm *

module load ${picard_module}
module load ${bwa_module}
module load ${samtools_module}
module load ${GATK_module}

#######################################################################################
# create a uBAM file
#######################################################################################

for file in ${raw_data}/*_R1_001.fastq

do

FBASE=$(basename $file _R1_001.fastq)
BASE=${FBASE%_R1_001.fastq}
java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
    FASTQ=${raw_data}/${BASE}_R1_001.fastq \
    FASTQ2=${raw_data}/${BASE}_R2_001.fastq  \
    OUTPUT=${raw_data}/${BASE}_fastqtosam.bam \
    READ_GROUP_NAME=${BASE} \
    SAMPLE_NAME=${BASE} \
    LIBRARY_NAME=H0 \
    PLATFORM=illumina \
    SEQUENCING_CENTER=GGBC

done

#######################################################################################
# works: aligns samples to reference genome. Output is a .sam file
#######################################################################################


#
#  #index the ref genome
# bwa index ${ref_genome}
# #
# for file in ${trimmed_data}/*_R1_001_trimmed.fq
#
# do
#
# FBASE=$(basename $file _R1_001_trimmed.fq)
# BASE=${FBASE%_R1_001_trimmed.fq}
#
# bwa mem -M -p -t 12 ${ref_genome} ${trimmed_data}/${BASE}_R1_001_trimmed.fq ${trimmed_data}/${BASE}_R2_001_trimmed.fq > ${output_directory}/${BASE}_aln.sam
#
# done

# #########################################################################################
# #samtools: converts sam files to bam files and sorts them
# #########################################################################################
#


# #index reference genome
#
# samtools faidx ${ref_genome}
# #
# # # #convert sam files to bam files
# for file in ${output_directory}/*_aln.sam
#
# do
#
# FBASE=$(basename $file _aln.sam)
# BASE=${FBASE%_aln.sam}
#
# samtools view -bt ${ref_genome_dir}/*.fai \
# ${output_directory}/${BASE}_aln.sam \
#   > ${output_directory}/${BASE}.bam
#
# done


# ############################
# ### sort the bam files
# ############################

# for file in ${output_directory}/*.bam
#
# do
#
# FBASE=$(basename $file .bam)
# BASE=${FBASE%.bam}
#
# samtools sort -@ 12 -o ${output_directory}/${BASE}.sorted.bam \
#    ${output_directory}/${BASE}.bam
#
# done

# ############################
# ### index the bam files
# ############################

# for file in ${output_directory}/*.sorted.bam
#
# do
#
# FBASE=$(basename $file .sorted.bam)
# BASE=${FBASE%.sorted.bam}
#
# samtools index -@ 12 ${output_directory}/${BASE}.sorted.bam
#
# done

###################################################################################################
# ## Picard to Validate Sam Files and mark duplicates
# ###################################################################################################
# #
#
# #
# for file in ${output_directory}/*.sorted.bam
#
# do
#
# FBASE=$(basename $file .sorted.bam)
# BASE=${FBASE%.sorted.bam}
#
# time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar ValidateSamFile \
#       I=${output_directory}/${BASE}.sorted.bam \
#       IGNORE_WARNINGS=true \
#       MODE=VERBOSE
#
# done
##################################################################################################
# for file in ${output_directory}/*.sorted.bam
#
# do
#
# FBASE=$(basename $file .sorted.bam)
# BASE=${FBASE%.sorted.bam}
#
# time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/${BASE}.sorted.bam \
# O=${output_directory}/${BASE}_removedDuplicates.bam \
# M=${output_directory}/${BASE}_removedDupsMetrics.txt
#
# done

# #### check the bam files again
# for file in ${output_directory}/*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
#
# time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar ValidateSamFile \
#       I=${output_directory}/${BASE}_removedDuplicates.bam \
#       IGNORE_WARNINGS=true \
#       MODE=VERBOSE
#
# done
##################################################################################################
###################################################################################################
# Using GATK HaplotypeCaller in GVCF mode
# apply appropriate ploidy for each sample
# will need to do this separtely for haploid and diploid samples
###################################################################################################



### H0 samples
# for file in ${output_directory}/*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
#
#
# time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -I ${output_directory}/${BASE}_removedDuplicates.bam \
#      -ploidy 1 \
#      -O ${output_directory}/${BASE}_variants.g.vcf
#
# done
# -ERC GVCF \


###################################################################################################
### Aggregate the GVCF files using GenomicsDBImport
###################################################################################################
# mkdir ${genomicsdb_workspace_path}
# mkdir ${tmp_DIR}
#
# gatk  --java-options "-Xmx4g -Xms4g" \
#        GenomicsDBImport \
#        --genomicsdb-workspace-path ${genomicsdb_workspace_path} \
#        --batch-size 50 \
#        --sample-name-map ${sample_name_map} \
#        --TMP_DIR:${tmp_DIR}/tmp \
#        --reader-threads 12
