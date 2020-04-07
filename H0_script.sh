#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N genotype
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=72:00:00
#PBS -l mem=200gb
#PBS -M hcm14449@uga.edu
#PBS -m abe

#S paradoxus TE MA Quality Control and Mutation Calling pipeline

#location of trimgalore moedule
trimgalore_module="Trim_Galore/0.4.5-foss-2016b"
#location of fastqc module
fastqc_module="FastQC/0.11.8-Java-1.8.0_144"
#location of BWA module
bwa_module="BWA/0.7.15-foss-2016b"
#location of samtools module
samtools_module="SAMtools/1.6-foss-2016b"
#location of bedtools module
bedtools_module="BEDTools/2.28.0-foss-2018a"
#location of python module
python_module="Python/3.5.2-foss-2016b"
#location of picard module
picard_module="picard/2.21.6-Java-11"
#location of GATK module
GATK_module="GATK/4.0.3.0-Java-1.8.0_144"
#deeptools location
deeptools_module="deepTools/3.2.1-foss-2018a-Python-3.6.4"
#location of bam to bigwig script
bamToBigWig="/scratch/hcm14449/TE_MA_Paradoxus/bedGraphToBigWigScript/file_to_bigwig_pe.py"
#location of reference genome to be used
ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/genome.337.fasta"
#directory reference genome is located in
ref_genome_dir="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/"
#where should the output be sent
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0"
# mkdir $output_directory
#location of data to be used in the analysis
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas"
mcc_bams="/scratch/jc33471/paradoxusHolly/run0217all"
mcc_bam_indiv="/scratch/jc33471/paradoxusHolly/run0217all/out/Spar"


# cd ${raw_data}
cd ${output_directory}

module load ${picard_module}
module load ${bwa_module}
module load ${samtools_module}
module load ${GATK_module}
#
# ### Much of the following obtained from https://software.broadinstitute.org/gatk/documentation/article?id=6483#step3
# #######################################################################################
# # create a uBAM file
# #######################################################################################
# #
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

# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
#     FASTQ=${raw_data}/HM-H0-11_R1_001.fastq.gz \
#     FASTQ2=${raw_data}/HM-H0-11_R2_001.fastq.gz  \
#     OUTPUT=${output_directory}/GATK_workflow_files/HM-H0-11_fastqtosam.bam \
#     READ_GROUP_NAME=HM-H0-11 \
#     SAMPLE_NAME=HM-H0-11 \
#     LIBRARY_NAME=H0 \
#     PLATFORM=illumina \
#     SEQUENCING_CENTER=GGBC

# for file in ${raw_data}/*_R1_001.fastq.gz
#
# do
#   FBASE=$(basename $file _R1_001.fastq.gz)
#   BASE=${FBASE%_R1_001.fastq.gz}
# 	OUT="${BASE}_FastqToSam.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_FastqToSam" >> ${OUT}
# 	echo "#PBS -l walltime=12:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# 	echo "#PBS -q batch" >> ${OUT}
# 	echo "#PBS -l mem=40gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${raw_data}" >> ${OUT}
# 	echo "module load ${picard_module}" >> ${OUT}
#   echo "module load ${bwa_module}" >> ${OUT}
#   echo "module load ${samtools_module}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
#       FASTQ=${raw_data}/${BASE}_R1_001.fastq.gz \
#       FASTQ2=${raw_data}/${BASE}_R2_001.fastq.gz  \
#       OUTPUT=${raw_data}/${BASE}_fastqtosam.bam \
#       READ_GROUP_NAME=${BASE} \
#       SAMPLE_NAME=${BASE} \
#       PLATFORM=illumina \
#       SEQUENCING_CENTER=GGBC" >> ${OUT}
# 	qsub ${OUT}
# done
#
# #
# #
# #
# #
# #
#
# for file in ${raw_data}/*R1.fq.gz
#
# do
#   FBASE=$(basename $file R1.fq.gz)
#   BASE=${FBASE%R1.fq.gz}
# 	OUT="${BASE}_FastqToSam.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_FastqToSam" >> ${OUT}
# 	echo "#PBS -l walltime=12:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# 	echo "#PBS -q batch" >> ${OUT}
# 	echo "#PBS -l mem=40gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${raw_data}" >> ${OUT}
# 	echo "module load ${picard_module}" >> ${OUT}
#   echo "module load ${bwa_module}" >> ${OUT}
#   echo "module load ${samtools_module}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
#       FASTQ=${raw_data}/${BASE}R1.fq.gz \
#       FASTQ2=${raw_data}/${BASE}R2.fq.gz \
#       OUTPUT=${output_directory}/${BASE}_fastqtosam.bam \
#       READ_GROUP_NAME=${BASE} \
#       SAMPLE_NAME=${BASE} \
#       PLATFORM=illumina \
#       SEQUENCING_CENTER=GGBC" >> ${OUT}
# 	qsub ${OUT}
# done



# for file in ${mcc_bam_indiv}/*_val/bam/*_val.bam;
#
# do
#
# FBASE=$(basename $file _val.bam)
# BASE=${FBASE%_val.bam}
# OUT="${BASE}_RevertSam.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_RevertSam" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=80gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
# /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar RevertSam \
# I=${mcc_bam_indiv}/${BASE}_val/bam/${BASE}_val.bam \
# O=${output_directory}/mcc_bams_out/${BASE}_Reverted.bam \
# SANITIZE=true \
# MAX_DISCARD_FRACTION=0.005 \
# ATTRIBUTE_TO_CLEAR=XT \
# ATTRIBUTE_TO_CLEAR=XN \
# ATTRIBUTE_TO_CLEAR=AS \
# ATTRIBUTE_TO_CLEAR=OC \
# ATTRIBUTE_TO_CLEAR=OP \
# SORT_ORDER=queryname \
# RESTORE_ORIGINAL_QUALITIES=true \
# REMOVE_DUPLICATE_INFORMATION=true \
# REMOVE_ALIGNMENT_INFORMATION=true" >> ${OUT}
#
# qsub ${OUT}
#
# done





# # #######################################################################################
# # # mark Illumina adapters
# # #######################################################################################
# #
mkdir ${raw_data}/TMP

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
#
#
#

# for file in ${output_directory}/*_fastqtosam.bam
#
# do
#   FBASE=$(basename $file _fastqtosam.bam)
#   BASE=${FBASE%_fastqtosam.bam}
# 	OUT="${BASE}_MarkIlluminaAdapters.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_MarkIlluminaAdapters" >> ${OUT}
# 	echo "#PBS -l walltime=12:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# 	echo "#PBS -q batch" >> ${OUT}
# 	echo "#PBS -l mem=50gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}/D0/" >> ${OUT}
# 	echo "module load ${picard_module}" >> ${OUT}
#   echo "module load ${bwa_module}" >> ${OUT}
#   echo "module load ${samtools_module}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "mkdir ${output_directory}/D0/TMP" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
#   I=${output_directory}/${BASE}_fastqtosam.bam \
#   O=${output_directory}/${BASE}_markilluminaadapters.bam \
#   M=${output_directory}/${BASE}_markilluminaadapters_metrics.txt \
#   TMP_DIR=${output_directory}/TMP" >> ${OUT}
# 	qsub ${OUT}
# done

############################################################

# for file in ${output_directory}/mcc_bams_out/*_Reverted.bam
#
# do
#   FBASE=$(basename $file _Reverted.bam)
#   BASE=${FBASE%_Reverted.bam}
# 	OUT="${BASE}_MarkIlluminaAdapters.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_MarkIlluminaAdapters" >> ${OUT}
# 	echo "#PBS -l walltime=12:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# 	echo "#PBS -q batch" >> ${OUT}
# 	echo "#PBS -l mem=50gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}/mcc_bams_out/" >> ${OUT}
# 	echo "module load ${picard_module}" >> ${OUT}
#   echo "module load ${bwa_module}" >> ${OUT}
#   echo "module load ${samtools_module}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "mkdir ${output_directory}/mcc_bams_out/TMP" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
#   I=${output_directory}/mcc_bams_out/${BASE}_Reverted.bam \
#   O=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters.bam \
#   M=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters_metrics.txt \
#   TMP_DIR=${output_directory}/mcc_bams_out/TMP" >> ${OUT}
# 	qsub ${OUT}
# done


# # #######################################################################################
# # # convert BAM to FASTQ and discount adapter sequences using SamToFastq **** THIS IS IN PIPED COMMAND *****
# # #######################################################################################
# #
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
# #
# #######################################################################################
# # works: aligns samples to reference genome. Output is a .sam file ** THIS IS ALSO IN PIPED COMMAND **
# #######################################################################################
#
#  # index the ref genome
# bwa index ${ref_genome}
#
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
# done
#
#
# #
# for file in ${output_directory}/*_markilluminaadapters.bam
#
# do
#   FBASE=$(basename $file _markilluminaadapters.bam)
#   BASE=${FBASE%_markilluminaadapters.bam}
# 	OUT="${BASE}_SamToFastq.sh"
# 	echo "#!/bin/bash" >> ${OUT}
# 	echo "#PBS -N ${BASE}_SamToFastq" >> ${OUT}
# 	echo "#PBS -l walltime=12:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# 	echo "#PBS -q batch" >> ${OUT}
# 	echo "#PBS -l mem=10gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}" >> ${OUT}
# 	echo "module load ${picard_module}" >> ${OUT}
#   echo "module load ${bwa_module}" >> ${OUT}
#   echo "module load ${samtools_module}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
#   I=${output_directory}/${BASE}_markilluminaadapters.bam \
#   FASTQ=${output_directory}/${BASE}_samtofastq_interleaved.fq \
#   CLIPPING_ATTRIBUTE=XT \
#   CLIPPING_ACTION=2 \
#   INTERLEAVE=true \
#   NON_PF=true \
#   TMP_DIR=${output_directory}/TMP" >> ${OUT}
# 	qsub ${OUT}
# done
#
#
# #######################################################################################
# # Piped command: SamToFastq, then bwa mem, then MergeBamAlignment
## #######################################################################################
#
# # java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# # /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar CreateSequenceDictionary \
# #       R=${ref_genome} \
# #       O=${ref_genome_dir}/genome.337.dict
#


#ask unix to stop piped command if any of it fails and report errors
# set -o pipefail
#
# for file in ${output_directory}/D0/*_markilluminaadapters.bam
#
# do
#
# FBASE=$(basename $file _markilluminaadapters.bam)
# BASE=${FBASE%_markilluminaadapters.bam}
# OUT="${BASE}_pipedNewRef.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_pipedNewRef" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=150gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D0" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
#   I=${output_directory}/D0/${BASE}_markilluminaadapters.bam \
#   FASTQ=/dev/stdout \
#   CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
#   TMP_DIR=${output_directory}/D0/TMP | \
#   bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
#   java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
#   ALIGNED_BAM=/dev/stdin \
#   UNMAPPED_BAM=${output_directory}/D0/${BASE}_fastqtosam.bam \
#   OUTPUT=${output_directory}/D0/${BASE}_pipedNewRef.bam \
#   R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
#   CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
#   INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
#   PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
#   TMP_DIR=${output_directory}/D0/TMP" >> ${OUT}
# qsub ${OUT}
#
# done

######################################################################################
# for file in ${output_directory}/mcc_bams_out/*_markilluminaadapters.bam
#
# do
#
# FBASE=$(basename $file _markilluminaadapters.bam)
# BASE=${FBASE%_markilluminaadapters.bam}
# OUT="${BASE}_pipedNewRef.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_pipedNewRef" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=150gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
#   I=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters.bam \
#   FASTQ=/dev/stdout \
#   CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
#   TMP_DIR=${output_directory}/mcc_bams_out/TMP | \
#   bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
#   java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
#   ALIGNED_BAM=/dev/stdin \
#   UNMAPPED_BAM=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters.bam \
#   OUTPUT=${output_directory}/mcc_bams_out/${BASE}_pipedNewRef.bam \
#   R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
#   CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
#   INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
#   PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
#   TMP_DIR=${output_directory}/mcc_bams_out/TMP" >> ${OUT}
# qsub ${OUT}
#
# done





# Sort the piped command output

# for file in ${output_directory}/D0/*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
# OUT="${BASE}_sortSam.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_sortSam" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D0" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#    /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
#     INPUT=${output_directory}/D0/${BASE}_pipedNewRef.bam \
#     OUTPUT=${output_directory}/D0/${BASE}_sorted.bam \
#     SORT_ORDER=coordinate" >> ${OUT}
#
# qsub ${OUT}
#
# done

# # ############################


# for file in ${output_directory}/H0/*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
# OUT="${BASE}_sortSam.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_sortSam" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/H0" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#    /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
#     INPUT=${output_directory}/H0/${BASE}_pipedNewRef.bam \
#     OUTPUT=${output_directory}/H0/${BASE}_sorted.bam \
#     SORT_ORDER=coordinate" >> ${OUT}
#
# qsub ${OUT}
#
# done


# # ############################


# for file in ${output_directory}/D1/*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
# OUT="${BASE}_sortSam.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_sortSam" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D1" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#    /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
#     INPUT=${output_directory}/D1/${BASE}_pipedNewRef.bam \
#     OUTPUT=${output_directory}/D1/${BASE}_sorted.bam \
#     SORT_ORDER=coordinate" >> ${OUT}
#
# qsub ${OUT}
#
# done

# # ############################


# for file in ${output_directory}/D20/*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
# OUT="${BASE}_sortSam.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_sortSam" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D20" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#    /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
#     INPUT=${output_directory}/D20/${BASE}_pipedNewRef.bam \
#     OUTPUT=${output_directory}/D20/${BASE}_sorted.bam \
#     SORT_ORDER=coordinate" >> ${OUT}
#
# qsub ${OUT}
#
# done

# # ############################


# for file in ${output_directory}/mcc_bams/*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
# OUT="${BASE}_sortSam.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_sortSam" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/mcc_bams" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#    /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
#     INPUT=${output_directory}/mcc_bams/${BASE}_pipedNewRef.bam \
#     OUTPUT=${output_directory}/mcc_bams/${BASE}_sorted.bam \
#     SORT_ORDER=coordinate" >> ${OUT}
#
# qsub ${OUT}
#
# done
# # ############################
# # ### index the bam files
# # ############################
#

for file in /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/mcc_bams_out/*_sorted.bam

do

FBASE=$(basename $file _sorted.bam)
BASE=${FBASE%_sorted.bam}

FBASE=$(basename $file _sorted.bam)
BASE=${FBASE%_sorted.bam}
OUT="${BASE}_index.sh"
echo "#!/bin/bash" >> ${OUT}
echo "#PBS -N ${BASE}_index" >> ${OUT}
echo "#PBS -l walltime=12:00:00" >> ${OUT}
echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
echo "#PBS -q batch" >> ${OUT}
echo "#PBS -l mem=50gb" >> ${OUT}
echo "" >> ${OUT}
echo "cd /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/mcc_bams_out" >> ${OUT}
echo "module load ${samtools_module}" >> ${OUT}
echo "" >> ${OUT}
echo "samtools index -@ 12/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/mcc_bams_out/${BASE}_sorted.bam" >> ${OUT}
qsub ${OUT}
done


#
# ###################################################################################################
# # ## Picard to Validate Sam Files and mark duplicates
# # ###################################################################################################
# # #
# #
# Sort the
# java -jar picard.jar SortSam \
#     INPUT=aligned_reads.sam \
#     OUTPUT=sorted_reads.bam \
#     SORT_ORDER=coordinate
# # #
# # for file in ${output_directory}/*.sorted.bam
# #
# # do
# #
# # FBASE=$(basename $file .sorted.bam)
# # BASE=${FBASE%.sorted.bam}
# #
# # time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# # /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar ValidateSamFile \
# #       I=${output_directory}/${BASE}.sorted.bam \
# #       IGNORE_WARNINGS=true \
# #       MODE=VERBOSE
# #
# # done
# ##################################################################################################
# for file in ${output_directory}/D0/${BASE}*_sorted.bam
#
# do
#
# FBASE=$(basename $file _sorted.bam)
# BASE=${FBASE%_sorted.bam}
# OUT="${BASE}_sorted.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_sorted" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D0" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
# /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/D0/${BASE}_sorted.bam \
# O=${output_directory}/D0/${BASE}_removedDuplicates.bam \
# M=${output_directory}/D0/${BASE}_removedDupsMetrics.txt" >> ${OUT}
# qsub ${OUT}
#
# done


# ##################################################################################################
# for file in ${output_directory}/D1/${BASE}*_sorted.bam
#
# do
#
# FBASE=$(basename $file _sorted.bam)
# BASE=${FBASE%_sorted.bam}
# OUT="${BASE}_markduplicates.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_markduplicates" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D1" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
# /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/D1/${BASE}_sorted.bam \
# O=${output_directory}/D1/${BASE}_removedDuplicates.bam \
# M=${output_directory}/D1/${BASE}_removedDupsMetrics.txt" >> ${OUT}
# qsub ${OUT}
#
# done

# ##################################################################################################
# for file in ${output_directory}/D20/${BASE}*_sorted.bam
#
# do
#
# FBASE=$(basename $file _sorted.bam)
# BASE=${FBASE%_sorted.bam}
# OUT="${BASE}_markduplicates.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_markduplicates" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D20" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
# /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/D20/${BASE}_sorted.bam \
# O=${output_directory}/D20/${BASE}_removedDuplicates.bam \
# M=${output_directory}/D20/${BASE}_removedDupsMetrics.txt" >> ${OUT}
# qsub ${OUT}
#
# done

#####################################################################
# for file in ${output_directory}/H0/${BASE}*_sorted.bam
#
# do
#
# FBASE=$(basename $file _sorted.bam)
# BASE=${FBASE%_sorted.bam}
# OUT="${BASE}_sorted.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_sorted" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=40gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/H0" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
# /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/H0/${BASE}_sorted.bam \
# O=${output_directory}/H0/${BASE}_removedDuplicates.bam \
# M=${output_directory}/H0/${BASE}_removedDupsMetrics.txt" >> ${OUT}
# qsub ${OUT}
#
# done

#####################################################################
# for file in ${output_directory}/mcc_bams_out/${BASE}*_sorted.bam
#
# do
#
# FBASE=$(basename $file _sorted.bam)
# BASE=${FBASE%_sorted.bam}
# OUT="${BASE}_markduplicates.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_markduplicates" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=40gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/H0" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
# /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/mcc_bams_out/${BASE}_sorted.bam \
# O=${output_directory}/mcc_bams_out/${BASE}_removedDuplicates.bam \
# M=${output_directory}/mcc_bams_out/${BASE}_removedDupsMetrics.txt" >> ${OUT}
# qsub ${OUT}
#
# done

# for file in ${output_directory}/D0/${BASE}*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
# OUT="${BASE}_HC.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_HC" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D0" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/D0/${BASE}_pipedNewRef.bam \
# O=${output_directory}/D0/${BASE}_removedDuplicates.bam \
# M=${output_directory}/D0/${BASE}_removedDupsMetrics.txt" >> ${OUT}
# qsub ${OUT}
#
# done


# for file in ${output_directory}/D1/${BASE}*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
# OUT="${BASE}_removeDuplicates.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_removeDuplicates" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D1" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/D1/${BASE}_pipedNewRef.bam \
# O=${output_directory}/D1/${BASE}_removedDuplicates.bam \
# M=${output_directory}/D1/${BASE}_removedDupsMetrics.txt" >> ${OUT}
# qsub ${OUT}
#
# done


# for file in ${output_directory}/D20/${BASE}*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
# OUT="${BASE}_removedDuplicates.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_removeDuplicates" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D20" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/D20/${BASE}_pipedNewRef.bam \
# O=${output_directory}/D20/${BASE}_removedDuplicates.bam \
# M=${output_directory}/D20/${BASE}_removedDupsMetrics.txt" >> ${OUT}
# qsub ${OUT}
#
# done
# #
# #
# # java -jar picard.jar BuildBamIndex \
# #     INPUT=dedup_reads.bam



##############################################################################################
############## ** Picard to BuildBamIndex              #################################
#################################################################################################

# for file in ${output_directory}/mcc_bams_out/*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_index.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_index" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
#    /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
#     INPUT=${output_directory}/mcc_bams_out/${BASE}_removedDuplicates.bam " >> ${OUT}
#
# qsub ${OUT}
#
# done

# for file in ${output_directory}/H0/*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_index.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_index" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/H0" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
#    /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
#     INPUT=${output_directory}/H0/${BASE}_removedDuplicates.bam " >> ${OUT}
#
# qsub ${OUT}
#
# done

# for file in ${output_directory}/D20/*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_index.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_index" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D20" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
#    /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
#     INPUT=${output_directory}/D20/${BASE}_removedDuplicates.bam " >> ${OUT}
#
# qsub ${OUT}
#
# done



# for file in ${output_directory}/D1/*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_index.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_index" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D1" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
#    /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
#     INPUT=${output_directory}/D1/${BASE}_removedDuplicates.bam " >> ${OUT}
#
# qsub ${OUT}
#
# done


# for file in ${output_directory}/D0/*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_index.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_index" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D0" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
#    /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
#     INPUT=${output_directory}/D0/${BASE}_removedDuplicates.bam " >> ${OUT}
#
# qsub ${OUT}
#
# done
# ##################################################################################################
# ###################################################################################################
# # Using GATK HaplotypeCaller in GVCF mode
# # apply appropriate ploidy for each sample
# # will need to do this separtely for haploid and diploid samples
# ###################################################################################################
# # ###################################################################################################
# # #
module load ${GATK_module}



# need to make folder for H0 in the mcc_bams_out because different ploidy

# for file in ${output_directory}/mcc_bams_out/H0/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_HC.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_HC" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/mcc_bams_out/H0" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I /${output_directory}/mcc_bams_out/H0/${BASE}_removedDuplicates.bam \
#      -ploidy 1 \
#      -O ${output_directory}/mcc_bams_out/H0/${BASE}_variantsNewRef.g.vcf" >> ${OUT}
# qsub ${OUT}
#
# done


# for file in ${output_directory}/mcc_bams_out/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_HC.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_HC" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I /${output_directory}/mcc_bams_out/${BASE}_removedDuplicates.bam \
#      -ploidy 2 \
#      -O ${output_directory}/mcc_bams_out/${BASE}_variantsNewRef.g.vcf" >> ${OUT}
# qsub ${OUT}
#
# done




# for file in ${output_directory}/D20/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_HC.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_HC" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D20" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I /${output_directory}/D20/${BASE}_removedDuplicates.bam \
#      -ploidy 2 \
#      -O ${output_directory}/D20/${BASE}_variantsNewRef.g.vcf" >> ${OUT}
# qsub ${OUT}
#
# done



#
# for file in ${output_directory}/D0/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_samToolsIndex.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_samToolsIndex" >> ${OUT}
# echo "#PBS -l walltime=24:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D0" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "samtools index -@ 12 ${output_directory}/D0/${BASE}*_removedDuplicates.bam" >> ${OUT}
# qsub ${OUT}
#
# done



# for file in ${output_directory}/D0/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_Hap.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_Hap" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D0" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I /${output_directory}/D0/${BASE}_removedDuplicates.bam \
#      -ploidy 2 \
#      -O ${output_directory}/D0/${BASE}_variants.g.vcf" >> ${OUT}
# qsub ${OUT}
#
# done





# for file in ${output_directory}/D1/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_HC.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_HC" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/D1" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I /${output_directory}/D1/${BASE}_removedDuplicates.bam \
#      -ploidy 2 \
#      -O ${output_directory}/D1/${BASE}_variantsNewRef.g.vcf" >> ${OUT}
# qsub ${OUT}
#
# done

# for file in ${output_directory}/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_HC.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_HC" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I ${output_directory}/${BASE}_removedDuplicates.bam \
#      -ploidy 1 \
#      -O ${output_directory}/${BASE}_variantsNewRef.g.vcf" >> ${OUT}
# qsub ${OUT}
#
# done


















# time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I ${output_directory}/H0/HM-H0-10_pipedNewRef.bam \
#      -ploidy 1 \
#      -O ${output_directory}/H0/HM-H0-10_variants.g.vcf

# # # time gatk HaplotypeCaller \
# #           -R ${ref_genome} \
# #           -ERC GVCF \
# #           -I ${raw_data}/HM-H0-45_piped.bam \
# #           -ploidy 1 \
# #           -O ${output_directory}/HM-H0-45_variants.g.vcf
# # time gatk HaplotypeCaller \
# #                -R ${ref_genome} \
# #                -ERC GVCF \
# #                -I ${raw_data}/HM-H0-46_piped.bam \
# #                -ploidy 1 \
# #                -O ${output_directory}/HM-H0-46_variants.g.vcf
#
#                # time gatk HaplotypeCaller \
#                #      -R ${ref_genome} \
#                #      -ERC GVCF \
#                #      -I ${raw_data}/HM-H0-47_piped.bam \
#                #      -ploidy 1 \
#                #      -O ${output_directory}/HM-H0-47_variants.g.vcf
# #
#                     # time gatk HaplotypeCaller \
#                     #      -R ${ref_genome} \
#                     #      -ERC GVCF \
#                     #      -I ${raw_data}/HM-H0-48_piped.bam \
#                     #      -ploidy 1 \
#                     #      -O ${output_directory}/HM-H0-48_variants.g.vcf
#
# ###################################################################################################
### Combine gVCFs before joint genotyping
# ###################################################################################################


# time gatk CombineGVCFs \
#  -O ${output_directory}/H0_variants_8SamplesNewRef.g.vcf \
#  -R ${ref_genome} \
#  --variant ${output_directory}/H0-A__variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-H0-10_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-H0-11_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-H0-12_variantsNewRef.g.vcf \
#  --variant ${output_directory}/H0-13__variantsNewRef.g.vcf \
#  --variant ${output_directory}/H0-14__variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-H0-15_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-H0-16_variantsNewRef.g.vcf




###################################################################################################
### Jointly genotype 8 random samples to identify consensus sequences
###################################################################################################

time gatk GenotypeGVCFs \
        -ploidy 1 \
        -R ${ref_genome} \
        --variant ${output_directory}/H0_variants_8SamplesNewRef.g.vcf \
        -O ${output_directory}/H0_variants_8SamplesNewRef_ploidy1.vcf


# ###################################################################################################
# ## Recalibrate base quality scores in all samples to mask any likely consensus variants
# ###################################################################################################
#
# for file in ${output_directory}/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_BR.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_BR" >> ${OUT}
# echo "#PBS -l walltime=24:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk BaseRecalibrator \
# -I ${output_directory}/${BASE}_removedDuplicates.bam \
# --known-sites ${output_directory}/H0_variants_8SamplesNewRef_ploidy1.vcf \
# -O ${output_directory}/${BASE}_recal_dataP1.table \
# -R ${ref_genome}" >> ${OUT}
# qsub ${OUT}
#
# done


# ###################################################################################################
# ## Apply BQSR to bam files
# ###################################################################################################
#
module load ${GATK_module}

# for file in ${output_directory}/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# OUT="${BASE}_HC.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_HC" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "gatk ApplyBQSR \
# -R ${ref_genome} \
# -I ${output_directory}/${BASE}_removedDuplicates.bam \
# -bqsr ${output_directory}/${BASE}_recal_dataP1.table \
# -O ${output_directory}/${BASE}_recalibratedNewRefP1.bam" >> ${OUT}
# qsub ${OUT}
#
# done



  # ###################################################################################################
  ### Run HaplotypeCaller again on recalibrated samples
  # ###################################################################################################
  # ###################################################################################################
  # #

        # # D1 samples
        # for file in ${output_directory}/${BASE}*_recalibratedNewRef.bam
        #
        # do
        #
        # FBASE=$(basename $file _recalibratedNewRef.bam)
        # BASE=${FBASE%_recalibratedNewRef.bam}
        #
        # time gatk HaplotypeCaller \
        # -R ${ref_genome} \
        # -ERC GVCF \
        # -I ${output_directory}/${BASE}_recalibratedNewRef.bam \
        # -ploidy 1 \
        # -O ${output_directory}/${BASE}_variants.Recal.g.vcf
        # done
for file in ${output_directory}/${BASE}*_recalibratedNewRefP1.bam

do

FBASE=$(basename $file _recalibratedNewRefP1.bam)
BASE=${FBASE%_recalibratedNewRefP1.bam}
OUT="${BASE}_HaplotypeCaller.sh"
echo "#!/bin/bash" > ${OUT}
echo "#PBS -N ${BASE}_HaplotypeCaller" >> ${OUT}
echo "#PBS -l walltime=72:00:00" >> ${OUT}
echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
echo "#PBS -q highmem_q" >> ${OUT}
echo "#PBS -l mem=200gb" >> ${OUT}
echo "" >> ${OUT}
echo "cd ${output_directory}" >> ${OUT}
echo "module load ${GATK_module}" >> ${OUT}
echo "" >> ${OUT}
echo "time gatk HaplotypeCaller \
-R ${ref_genome} \
-ERC GVCF \
-I ${output_directory}/${BASE}_recalibratedNewRefP1.bam \
-ploidy 1 \
-O ${output_directory}/${BASE}_variants.RecalP1.g.vcf" >> ${OUT}
qsub ${OUT}

done


#
#           ###################################################################################################
#           ## Combine gvcfs
#           ###################################################################################################
#           ###################################################################################################
#           #
#           #
module load ${GATK_module}
#
time gatk CombineGVCFs \
-R ${ref_genome} \
-O ${output_directory}/H0_FullCohort.g.vcf \
-V ${output_directory}/H0-A__variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-1_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-2_variants.Recal.g.vcf \
-V ${output_directory}/H0-3__variants.Recal.g.vcf \
-V ${output_directory}/H0-4__variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-5_variants.Recal.g.vcf \
-V ${output_directory}/H0-7__variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-8_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-9_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-10_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-11_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-12_variants.Recal.g.vcf \
-V ${output_directory}/H0-13__variants.Recal.g.vcf \
-V ${output_directory}/H0-14__variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-15_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-16_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-17_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-18_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-19_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-20_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-21_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-22_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-23_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-24_variants.Recal.g.vcf \
-V ${output_directory}/H0-25__variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-26_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-27_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-28_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-29_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-30_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-31_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-32_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-33_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-34_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-35_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-36_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-37_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-38_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-39_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-40_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-41_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-42_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-43_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-44_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-45_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-46_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-47_variants.Recal.g.vcf \
-V ${output_directory}/HM-H0-48_variants.Recal.g.vcf
#
#              ###################################################################################################
#              ## Genotype gVCFs (jointly)
#              ###################################################################################################
#              ###################################################################################################
#
#
time gatk GenotypeGVCFs \
-R ${ref_genome} \
-ploidy 1 \
--variant ${output_directory}/H0_FullCohort.g.vcf \
-O ${output_directory}/H0_FullCohort.vcf




#### Remove sites with mappability < 0.9
low_mappability="/scratch/jc33471/pilon/337/mappability/337_lowmappability.bed"
module load ${bedtools_module}

# bedtools sort -i ${low_mappability} > ${output_directory}/337_lowmappability_sorted.bed
bedtools intersect -v -a ${output_directory}/H0_FullCohort.vcf -b ${low_mappability} -header > ${output_directory}/H0_reducedGEM.vcf

#count number of lines between original vcf and reduced vcf
wc -l ${output_directory}/D0_FullCohort.vcf
#1595
wc -l ${output_directory}/reducedTest.vcf
#675

# ###################################################################################################
# ### Find coverage and put into 10k chunks
# ###################################################################################################

# module load ${deeptools_module}
#
#
# for file in ${raw_data}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
# OUT="${BASE}_bamCoverage.sh"
# echo "#!/bin/bash" > ${OUT}
# echo "#PBS -N ${BASE}_bamCoverage" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=20gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${raw_data}" >> ${OUT}
# echo "module load ${deeptools_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "bamCoverage -b ${raw_data}/${BASE}_piped.bam -o ${output_directory}/${BASE}.bedgraph -of bedgraph -bs 10000" >> ${OUT}
# qsub ${OUT}
#
# done


# for file in ${raw_data}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
# samtools sort ${raw_data}/${BASE}_piped.bam \
# -o ${raw_data}/${BASE}.sorted.bam
#
# samtools depth \
# ${raw_data}/${BASE}.sorted.bam \
# |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${raw_data}/${BASE}.txt
#
# done
#
# #combine depths with filenames into the same file
# find . -type f -name "*.txt" -exec awk '{s=$0};END{if(s)print FILENAME,s}' {} \; > D0_depth.txt


#                  # ################
# # ###################################################################################################
# # ### Filter variants
# # Can easily run these interactively
# # ###################################################################################################
#
# Get only those lines where there is actually a genotype call in the ancestor
gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/H0_FullCohort.vcf \
-O ${output_directory}/H0_FullCohort_AncCalls.vcf \
-select 'vc.getGenotype("HM-D0-A").isCalled()'

#
# remove all lines in the ancestor that have a heterozygous genotype
gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/H0_FullCohort_AncCalls.vcf \
-O ${output_directory}/H0_FullCohort_AncCalls_NoHets.vcf \
-select '!vc.getGenotype("HM-D0-A").isHet()'

# filter out sites with low read depth
gatk VariantFiltration \
   -R ${ref_genome} \
   -V ${output_directory}/H0_FullCohort_AncCalls_NoHets.vcf \
   -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBias.vcf \
   --set-filtered-genotype-to-no-call TRUE \
   -G-filter "DP < 10"  -G-filter-name "depthGr10" \
   -filter "MQ < 50.0" -filter-name "MQ50" \
   -filter "SOR < 0.01" -filter-name "strandBias"

#
  # remove filtered sites (these were set to no calls ./.)
gatk SelectVariants \
   -R ${ref_genome} \
   -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBias.vcf \
   -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil.vcf \
   --exclude-filtered TRUE

gatk SelectVariants \
   -R ${ref_genome} \
   -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil.vcf \
   -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
   -select 'vc.getGenotype("HM-D0-A").isCalled()'

# cd ${output_directory}
#
gatk SelectVariants \
   -R ${ref_genome} \
   -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
   -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_SNPs.vcf \
   --max-nocall-fraction 0 \
   --exclude-non-variants TRUE \
   -select-type SNP

gatk SelectVariants \
   -R ${ref_genome} \
   -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
   -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_Indels.vcf \
   --max-nocall-fraction 0.001 \
   -select-type INDEL

#
#    -select-type INDEL \
#    -select-type MIXED \
#    -select-type MNP \
#    -select-type SYMBOLIC
#
# #
# # #gives a final dataset with only called sites in the Ancestor, no heterozygous sites in the ancestor,
# # # depth > 10, mapping quality > 50, and strand bias (SOR) > 0.01 (not significant)
# #
# # #Variants to table
gatk VariantsToTable \
     -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_vars.vcf \
     -F CHROM -F POS -F REF -F ALT -F QUAL \
     -GF AD -GF DP -GF GQ -GF GT \
     -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_vars.txt

gatk VariantsToTable \
      -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_SNPs.vcf \
      -F CHROM -F POS -F REF -F ALT  \
      -GF AD -GF GT \
      -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_SNPs.txt

gatk VariantsToTable \
      -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_Indels.vcf \
      -F CHROM -F POS -F REF -F ALT  \
      -GF AD -GF GT \
      -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_Indels.txt

# # ##################################################################################################
# # #### UNIFIED GENOTYPER
# # ##################################################################################################






# # #########################################################################################
# # #samtools: converts sam files to bam files and sorts them
# # #########################################################################################
# #
#
#
# # #index reference genome
# #
# samtools faidx ${ref_genome}
# # #
# # # # #convert sam files to bam files
# # for file in ${output_directory}/*_aln.sam
# #
# # do
# #
# # FBASE=$(basename $file _aln.sam)
# # BASE=${FBASE%_aln.sam}
# #
# # samtools view -bt ${ref_genome_dir}/*.fai \
# # ${output_directory}/${BASE}_aln.sam \
# #   > ${output_directory}/${BASE}.bam
# #
# # done
#
#
# # ############################
# # ### sort the bam files
# # ############################
#
# # for file in ${output_directory}/*.bam
# #
# # do
# #
# # FBASE=$(basename $file .bam)
# # BASE=${FBASE%.bam}
# #
# # samtools sort -@ 12 -o ${output_directory}/${BASE}.sorted.bam \
# #    ${output_directory}/${BASE}.bam
# #
# # done
#
### MCC Bam Files
# sort the bam files first
#    for file in ${mcc_bam_indiv}/*_val/bam/*_val.bam;
#
#    do
#
#    FBASE=$(basename $file _val.bam)
#    BASE=${FBASE%_val.bam}
#    OUT="${BASE}_sort.sh"
#    echo "#!/bin/bash" >> ${OUT}
#    echo "#PBS -N ${BASE}_sort" >> ${OUT}
#    echo "#PBS -l walltime=12:00:00" >> ${OUT}
#    echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
#    echo "#PBS -q batch" >> ${OUT}
#    echo "#PBS -l mem=20gb" >> ${OUT}
#    echo "" >> ${OUT}
#    echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
#    echo "module load ${samtools_module}" >> ${OUT}
#    echo "" >> ${OUT}
#    echo "samtools sort -@ 12 -o ${output_directory}/mcc_bams_out/${BASE}_sorted.bam \
#       ${mcc_bam_indiv}/${BASE}_val/bam/${BASE}_val.bam" >> ${OUT}
#    qsub ${OUT}
#
# done
