#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N D20_scripts
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l walltime=48:00:00
#PBS -l mem=300gb
#PBS -M hcm14449@uga.edu
#PBS -m abe

#S paradoxus TE MA Quality Control and Mutation Calling pipeline

#location of current update of muver
# muver_module="muver/0.1.0-foss-2016b-Python-2.7.14-20190318"
#location of trimgalore moedule
# trimgalore_module="Trim_Galore/0.4.5-foss-2016b"
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
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D20"
# mkdir $output_directory
#location of TRIMMED data to be used in the analysis
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/D20"
trimmed_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/D20"
# raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq"
genomicsdb_workspace_path="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D20/GenDB"
sample_name_map="/home/hcm14449/Github/TE_MA/D20_sample_map.txt"
tmp_DIR="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D20/GenDB/tmp"

# cd ${output_directory}
# rm *

module load ${picard_module}
module load ${bwa_module}
module load ${samtools_module}
module load ${GATK_module}

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
#     LIBRARY_NAME=D20 \
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
# USE_JDK_DEFLATER=true \
# USE_JDK_INFLATER=true
#
# done


### Piped command: SamToFastq, then bwa mem, then MergeBamAlignment
for file in ${raw_data}/*_markilluminaadapters.bam

do

FBASE=$(basename $file _markilluminaadapters.bam)
BASE=${FBASE%_markilluminaadapters.bam}

java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
I=${raw_data}/${BASE}_markilluminaadapters.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=/path/shlee | \
bwa mem -M -t 7 -p ${ref_genome} /dev/stdin | \
java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=${raw_data}/${BASE}_fastqtosam.bam \
OUTPUT=${raw_data}/${BASE}_piped.bam \
R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=/path/shlee

done

# #######################################################################################
# # works: aligns samples to reference genome. Output is a .sam file
# #######################################################################################
#
# module load ${bwa_module}
# #
# #  #index the ref genome
# bwa index ${ref_genome}
# #
# for file in ${trimmed_data}/*_R1_001_trimmed.fq
#
# do
#
# FBASE=$(basename $file _R1_001_trimmed.fq)
# BASE=${FBASE%_R1_001_trimmed.fq}
#
# bwa mem -p -M -t 12 ${ref_genome} ${trimmed_data}/${BASE}_R1_001_trimmed.fq ${trimmed_data}/${BASE}_R2_001_trimmed.fq > ${output_directory}/${BASE}_aln.sam
#
# done
# # #########################################################################################
# # #samtools: converts sam files to bam files and sorts them
# # #########################################################################################
# #
# module load ${samtools_module}
#
#
# #convert sam files to bam files
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
#
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
#
# # ############################
# # ### index the bam files
# # ############################
#
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
#
# # ###################################################################################################
# # ## Picard to mark duplicates
# # ###################################################################################################
# #
# module load ${picard_module}
#
#
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
#       MODE=VERBOSE
#
# done
# ###################################################################################################
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
# ###################################################################################################
# # Using GATK HaplotypeCaller in GVCF mode
# # apply appropriate ploidy for each sample
# # will need to do this separtely for haploid and diploid samples
# ###################################################################################################
# #
# module load ${GATK_module}
#
# #### D20 samples
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
#      -ploidy 2 \
#      -O ${output_directory}/${BASE}_variants.g.vcf
#
# done
# # -ERC GVCF \
#
# ###################################################################################################
# ### Aggregate the GVCF files using GenomicsDBImport
# ###################################################################################################
# # mkdir ${genomicsdb_workspace_path}
# # mkdir ${tmp_DIR}
# #
# # gatk --java-options "-Xmx4g -Xms4g" \
# #        GenomicsDBImport \
# #        --genomicsdb-workspace-path ${genomicsdb_workspace_path} \
# #        --batch-size 50 \
# #        --sample-name-map ${sample_name_map} \
# #        --TMP_DIR: ${tmp_DIR} \
# #        --reader-threads 12
