#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N H0_scripts_base_recal
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l walltime=48:00:00
#PBS -l mem=250gb
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
picard_module="picard/2.4.1-Java-1.8.0_144"
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
do_again="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/H0/do_again"
genomicsdb_workspace_path="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB"
sample_name_map="/home/hcm14449/Github/TE_MA/H0_sample_map.txt"
tmp_DIR="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB/tmp"

# cd ${output_directory}
# rm *

# module load ${picard_module}
# module load ${bwa_module}
# module load ${samtools_module}
module load ${GATK_module}

### Much of the following obtained from https://software.broadinstitute.org/gatk/documentation/article?id=6483#step3
#######################################################################################
# create a uBAM file
#######################################################################################

# for file in ${raw_data}/*_R1_001.fastq
#
# do
#
# FBASE=$(basename $file _R1_001.fastq)
# BASE=${FBASE%_R1_001.fastq}
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
#     FASTQ=${raw_data}/${BASE}_R1_001.fastq \
#     FASTQ2=${raw_data}/${BASE}_R2_001.fastq  \
#     OUTPUT=${raw_data}/${BASE}_fastqtosam.bam \
#     READ_GROUP_NAME=${BASE} \
#     SAMPLE_NAME=${BASE} \
#     LIBRARY_NAME=H0 \
#     PLATFORM=illumina \
#     SEQUENCING_CENTER=GGBC
#
# done

#######################################################################################
# mark Illumina adapters
#######################################################################################

# mkdir ${raw_data}/TMP
#
# for file in ${raw_data}/*_fastqtosam.bam
#
# do
#
# FBASE=$(basename $file _fastqtosam.bam)
# BASE=${FBASE%_fastqtosam.bam}
#
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
# I=${raw_data}/${BASE}_fastqtosam.bam \
# O=${raw_data}/${BASE}_markilluminaadapters.bam \
# M=${raw_data}/${BASE}_markilluminaadapters_metrics.txt \
# TMP_DIR=${raw_data}/TMP \
# USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
#
# done

#######################################################################################
# convert BAM to FASTQ and discount adapter sequences using SamToFastq
#######################################################################################

# for file in ${raw_data}/*_markilluminaadapters.bam
#
# do
#
# FBASE=$(basename $file _markilluminaadapters.bam)
# BASE=${FBASE%_markilluminaadapters.bam}
#
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar SamToFastq \
# I=${raw_data}/${BASE}_markilluminaadapters.bam \
# FASTQ=${raw_data}/${BASE}_samtofastq_interleaved.fq \
# CLIPPING_ATTRIBUTE=XT \
# CLIPPING_ACTION=2 \
# INTERLEAVE=true \
# NON_PF=true \
# TMP_DIR=${raw_data}/TMP
#
# done

#######################################################################################
# works: aligns samples to reference genome. Output is a .sam file
#######################################################################################

#  #index the ref genome
# bwa index ${ref_genome}
# #
# for file in ${raw_data}/*_samtofastq_interleaved.fq
#
# do
#
# FBASE=$(basename $file _samtofastq_interleaved.fq)
# BASE=${FBASE%_samtofastq_interleaved.fq}
#
# bwa mem -M -p -t 12 ${ref_genome} ${raw_data}/${BASE}_samtofastq_interleaved.fq > ${output_directory}/${BASE}_bwa_mem.sam
#
#
#
# done


### Piped command: SamToFastq, then bwa mem, then MergeBamAlignment
# for file in ${raw_data}/*_markilluminaadapters.bam
#
# do
#
# FBASE=$(basename $file _markilluminaadapters.bam)
# BASE=${FBASE%_markilluminaadapters.bam}
#
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
# I=${raw_data}/${BASE}_markilluminaadapters.bam \
# FASTQ=/dev/stdout \
# CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
# TMP_DIR=${raw_data}/TMP | \
# bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
# ALIGNED_BAM=/dev/stdin \
# UNMAPPED_BAM=${raw_data}/${BASE}_fastqtosam.bam \
# OUTPUT=${raw_data}/${BASE}_piped.bam \
# R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
# CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
# INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
# PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
# TMP_DIR=${raw_data}/TMP
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
## Sort the
# java -jar picard.jar SortSam \
#     INPUT=aligned_reads.sam \
#     OUTPUT=sorted_reads.bam \
#     SORT_ORDER=coordinate
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
#
#
# java -jar picard.jar BuildBamIndex \
#     INPUT=dedup_reads.bam
##################################################################################################
###################################################################################################
# Using GATK HaplotypeCaller in GVCF mode
# apply appropriate ploidy for each sample
# will need to do this separtely for haploid and diploid samples
###################################################################################################
# ###################################################################################################
# #
module load ${GATK_module}

#### H0 samples
# for file in ${raw_data}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
#
# time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I ${raw_data}/${BASE}_piped.bam \
#      -ploidy 1 \
#      -O ${output_directory}/${BASE}_variants.g.vcf
#
# done


# ###################################################################################################
### Combine gVCFs before joint genotyping
# ###################################################################################################


# time gatk CombineGVCFs \
#  -O ${output_directory}/H0_cohort.g.vcf \
#  -R ${ref_genome} \
#  --variant ${output_directory}/HM-H0-A_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-10_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-11_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-12_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-13_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-14_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-15_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-16_variants.g.vcf


# ###################################################################################################
# ### Jointly genotype 8 random samples to identify consensus sequences
# ###################################################################################################

# time gatk GenotypeGVCFs \
#         -R ${ref_genome} \
#         --variant ${output_directory}/H0_cohort.g.vcf \
#         -O ${output_directory}/H0_variants_8Samples.vcf

# ###################################################################################################
# ## Recalibrate base quality scores in all samples to mask any likely consensus variants
# ###################################################################################################

for file in ${raw_data}/${BASE}*_piped.bam

do

FBASE=$(basename $file _piped.bam)
BASE=${FBASE%_piped.bam}


time gatk BaseRecalibrator \
   -I ${raw_data}/${BASE}_piped.bam \
   --known-sites ${output_directory}/H0_variants_8Samples.vcf \
   -O ${output_directory}/${BASE}_recal_data.table \
   -R ${ref_genome}

done



# ###################################################################################################
# ### Aggregate the GVCF files using GenomicsDBImport
# ###################################################################################################
# mkdir ${genomicsdb_workspace_path}
# mkdir ${tmp_DIR}

# gatk --java-options "-Xmx4g -Xms4g" \
#        GenomicsDBImport \
#        --genomicsdb-workspace-path ${genomicsdb_workspace_path} \
#        --batch-size 50 \
#        --sample-name-map ${sample_name_map} \
#        --TMP_DIR: ${tmp_DIR} \
#        --reader-threads 12
