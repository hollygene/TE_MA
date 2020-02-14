#PBS -S /bin/bash
#PBS -q batch
#PBS -N genotype
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=72:00:00
#PBS -l mem=50gb
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
picard_module="picard/2.4.1-Java-1.8.0_144"
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
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/H0"


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
#
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
#
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
#   echo "mkdir ${output_directory}/GATK_workflow_files" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
#       FASTQ=${raw_data}/${BASE}_R1_001.fastq.gz \
#       FASTQ2=${raw_data}/${BASE}_R2_001.fastq.gz  \
#       OUTPUT=${output_directory}/GATK_workflow_files/${BASE}_fastqtosam.bam \
#       READ_GROUP_NAME=${BASE} \
#       SAMPLE_NAME=${BASE} \
#       LIBRARY_NAME=H0 \
#       PLATFORM=illumina \
#       SEQUENCING_CENTER=GGBC" >> ${OUT}
# 	qsub ${OUT}
# done
#
#
#
#
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
#       FASTQ2=${raw_data}/${BASE}R1.fq.gz \
#       OUTPUT=${output_directory}/GATK_workflow_files/${BASE}_fastqtosam.bam \
#       READ_GROUP_NAME=${BASE} \
#       SAMPLE_NAME=${BASE} \
#       LIBRARY_NAME=H0 \
#       PLATFORM=illumina \
#       SEQUENCING_CENTER=GGBC" >> ${OUT}
# 	qsub ${OUT}
# done
# # #######################################################################################
# # # mark Illumina adapters
# # #######################################################################################
# #
# # mkdir ${raw_data}/TMP
# #
# # for file in ${raw_data}/*_fastqtosam.bam
# #
# # do
# #
# # FBASE=$(basename $file _fastqtosam.bam)
# # BASE=${FBASE%_fastqtosam.bam}
# #
# # java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# # /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
# # I=${raw_data}/${BASE}_fastqtosam.bam \
# # O=${raw_data}/${BASE}_markilluminaadapters.bam \
# # M=${raw_data}/${BASE}_markilluminaadapters_metrics.txt \
# # TMP_DIR=${raw_data}/TMP \
# # USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
# #
# # done
# #
#
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
# I=${output_directory}/GATK_workflow_files/HM-H0-11_fastqtosam.bam \
# O=${output_directory}/GATK_workflow_files/HM-H0-11_markilluminaadapters.bam \
# M=${output_directory}/GATK_workflow_files/HM-H0-11_markilluminaadapters_metrics.txt \
# TMP_DIR=${output_directory}/GATK_workflow_files/TMP \
# USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
#
#
#
# for file in ${output_directory}/GATK_workflow_files/*_fastqtosam.bam
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
# 	echo "#PBS -l mem=20gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}GATK_workflow_files/" >> ${OUT}
# 	echo "module load ${picard_module}" >> ${OUT}
#   echo "module load ${bwa_module}" >> ${OUT}
#   echo "module load ${samtools_module}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "mkdir ${output_directory}/GATK_workflow_files/TMP" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
#   I=${output_directory}/GATK_workflow_files/${BASE}_fastqtosam.bam \
#   O=${output_directory}/GATK_workflow_files/${BASE}_markilluminaadapters.bam \
#   M=${output_directory}/GATK_workflow_files/${BASE}_markilluminaadapters_metrics.txt \
#   TMP_DIR=${output_directory}/GATK_workflow_files/TMP" >> ${OUT}
# 	qsub ${OUT}
# done
# # #######################################################################################
# # # convert BAM to FASTQ and discount adapter sequences using SamToFastq
# # #######################################################################################
# #
# # # for file in ${raw_data}/*_markilluminaadapters.bam
# # #
# # # do
# # #
# # # FBASE=$(basename $file _markilluminaadapters.bam)
# # # BASE=${FBASE%_markilluminaadapters.bam}
# # #
# # # java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# # # /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar SamToFastq \
# # # I=${raw_data}/${BASE}_markilluminaadapters.bam \
# # # FASTQ=${raw_data}/${BASE}_samtofastq_interleaved.fq \
# # # CLIPPING_ATTRIBUTE=XT \
# # # CLIPPING_ACTION=2 \
# # # INTERLEAVE=true \
# # # NON_PF=true \
# # # TMP_DIR=${raw_data}/TMP
# # #
# # # done
# #
# #######################################################################################
# # works: aligns samples to reference genome. Output is a .sam file
# #######################################################################################
#
#  # index the ref genome
# bwa index ${ref_genome}
#
# # for file in ${raw_data}/*_samtofastq_interleaved.fq
# #
# # do
# #
# # FBASE=$(basename $file _samtofastq_interleaved.fq)
# # BASE=${FBASE%_samtofastq_interleaved.fq}
# #
# # bwa mem -M -p -t 12 ${ref_genome} ${raw_data}/${BASE}_samtofastq_interleaved.fq > ${output_directory}/${BASE}_bwa_mem.sam
# #
# #
# #
# # done
#
java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
I=${output_directory}/GATK_workflow_files/HM-H0-11_markilluminaadapters.bam \
FASTQ=${output_directory}/GATK_workflow_files/HM-H0-11_samtofastq_interleaved.fq \
CLIPPING_ATTRIBUTE=XT \
CLIPPING_ACTION=2 \
INTERLEAVE=true \
NON_PF=true \
TMP_DIR=${output_directory}/GATK_workflow_files/TMP
#
#
#
#
# for file in ${output_directory}/GATK_workflow_files/*_markilluminaadapters.bam
#
# do
#   FBASE=$(basename $file _markilluminaadapters.bam)
#   BASE=${FBASE%_markilluminaadapters.bam}
# 	OUT="${BASE}_SamToFastq.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_SamToFastq" >> ${OUT}
# 	echo "#PBS -l walltime=12:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# 	echo "#PBS -q batch" >> ${OUT}
# 	echo "#PBS -l mem=10gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}/GATK_workflow_files" >> ${OUT}
# 	echo "module load ${picard_module}" >> ${OUT}
#   echo "module load ${bwa_module}" >> ${OUT}
#   echo "module load ${samtools_module}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
#   I=${output_directory}/GATK_workflow_files/${BASE}_markilluminaadapters.bam \
#   FASTQ=${output_directory}/GATK_workflow_files/${BASE}_samtofastq_interleaved.fq \
#   CLIPPING_ATTRIBUTE=XT \
#   CLIPPING_ACTION=2 \
#   INTERLEAVE=true \
#   NON_PF=true \
#   TMP_DIR=${output_directory}/GATK_workflow_files/TMP" >> ${OUT}
# 	qsub ${OUT}
# done
#
#
# # Piped command: SamToFastq, then bwa mem, then MergeBamAlignment
#
# # java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# # /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar CreateSequenceDictionary \
# #       R=${ref_genome} \
# #       O=${ref_genome_dir}/genome.337.dict
#
#
#
java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
I=${output_directory}/GATK_workflow_files/HM-H0-11_markilluminaadapters.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=${output_directory}/GATK_workflow_files/TMP | \
bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=${output_directory}/GATK_workflow_files/HM-H0-11_fastqtosam.bam \
OUTPUT=${output_directory}/GATK_workflow_files/HM-H0-11_pipedNewRef.bam \
R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=${output_directory}/GATK_workflow_files/TMP

# for file in ${output_directory}/GATK_workflow_files/*_markilluminaadapters.bam
#
# do
#
# FBASE=$(basename $file _markilluminaadapters.bam)
# BASE=${FBASE%_markilluminaadapters.bam}
# OUT="${BASE}_pipedNewRef.sh"
# echo "#!/bin/bash" > ${OUT}
# echo "#PBS -N ${BASE}_pipedNewRef" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=150gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}/GATK_workflow_files/" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
# I=${output_directory}/GATK_workflow_files/${BASE}_markilluminaadapters.bam \
# FASTQ=/dev/stdout \
# CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
# TMP_DIR=${output_directory}GATK_workflow_files//TMP | \
# bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
# ALIGNED_BAM=/dev/stdin \
# UNMAPPED_BAM=${output_directory}/GATK_workflow_files/${BASE}_fastqtosam.bam \
# OUTPUT=${output_directory}/GATK_workflow_files/${BASE}_pipedNewRef.bam \
# R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
# CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
# INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
# PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
# TMP_DIR=${output_directory}/GATK_workflow_files/TMP" >> ${OUT}
# qsub ${OUT}
#
# done
# # #########################################################################################
# # #samtools: converts sam files to bam files and sorts them
# # #########################################################################################
# #
#
#
# # #index reference genome
# #
samtools faidx ${ref_genome}
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
# # ############################
# # ### index the bam files
# # ############################
#
# # for file in ${output_directory}/*.sorted.bam
# #
# # do
# #
# # FBASE=$(basename $file .sorted.bam)
# # BASE=${FBASE%.sorted.bam}
# #
# # samtools index -@ 12 ${output_directory}/${BASE}.sorted.bam
# #
# # done
#
# ###################################################################################################
# # ## Picard to Validate Sam Files and mark duplicates
# # ###################################################################################################
# # #
# #
# ## Sort the
# # java -jar picard.jar SortSam \
# #     INPUT=aligned_reads.sam \
# #     OUTPUT=sorted_reads.bam \
# #     SORT_ORDER=coordinate
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
# # for file in ${output_directory}/*.sorted.bam
# #
# # do
# #
# # FBASE=$(basename $file .sorted.bam)
# # BASE=${FBASE%.sorted.bam}
# #
# # time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# # /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
# # REMOVE_DUPLICATES=TRUE \
# # I=${output_directory}/${BASE}.sorted.bam \
# # O=${output_directory}/${BASE}_removedDuplicates.bam \
# # M=${output_directory}/${BASE}_removedDupsMetrics.txt
# #
# # done
# #
# #
# # java -jar picard.jar BuildBamIndex \
# #     INPUT=dedup_reads.bam
# ##################################################################################################
# ###################################################################################################
# # Using GATK HaplotypeCaller in GVCF mode
# # apply appropriate ploidy for each sample
# # will need to do this separtely for haploid and diploid samples
# ###################################################################################################
# # ###################################################################################################
# # #
module load ${GATK_module}

# H0 samples
# for file in ${output_directory}/GATK_workflow_files/${BASE}*_pipedNewRef.bam
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
# echo "cd ${output_directory}/GATK_workflow_files/" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I /${output_directory}/GATK_workflow_files/${BASE}_pipedNewRef.bam \
#      -ploidy 1 \
#      -O ${output_directory}/GATK_workflow_files/${BASE}_variantsNewRef.g.vcf" >> ${OUT}
# qsub ${OUT}
#
# done

time gatk HaplotypeCaller \
     -R ${ref_genome} \
     -ERC GVCF \
     -I ${output_directory}/GATK_workflow_files/HM-H0-10_pipedNewRef.bam \
     -ploidy 1 \
     -O ${output_directory}/GATK_workflow_files/HM-H0-10_variants.g.vcf

# # time gatk HaplotypeCaller \
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


time gatk CombineGVCFs \
 -O ${output_directory}/GATK_workflow_files/H0_smallCohortNewRef.g.vcf \
 -R ${ref_genome} \
 --variant ${output_directory}/GATK_workflow_files/H0-A_variantsNewRef.g.vcf \
 --variant ${output_directory}/GATK_workflow_files/HM-H0-10_variants.g.vcf \
 --variant ${output_directory}/GATK_workflow_files/HM-H0-11_variantsNewRef.g.vcf \
 --variant ${output_directory}/GATK_workflow_files/HM-H0-12_variantsNewRef.g.vcf \
 --variant ${output_directory}/GATK_workflow_files/H0-13_variantsNewRef.g.vcf \
 --variant ${output_directory}/GATK_workflow_files/H0-14_variantsNewRef.g.vcf \
 --variant ${output_directory}/GATK_workflow_files/HM-H0-15_variantsNewRef.g.vcf \
 --variant ${output_directory}/GATK_workflow_files/HM-H0-16_variantsNewRef.g.vcf


# ###################################################################################################
# ### Jointly genotype 8 random samples to identify consensus sequences
# ###################################################################################################

time gatk GenotypeGVCFs \
        -R ${ref_genome} \
        -ploidy 1 \
        --variant ${output_directory}/H0_smallCohortNewRef.g.vcf \
        -O ${output_directory}/H0_variants_8SamplesNewRef.vcf

# # ###################################################################################################
# # ## Recalibrate base quality scores in all samples to mask any likely consensus variants
# # ###################################################################################################
#
for file in ${output_directory}/${BASE}*_piped.bam

do

FBASE=$(basename $file _piped.bam)
BASE=${FBASE%_piped.bam}
OUT="${BASE}_BR.sh"
echo "#!/bin/bash" > ${OUT}
echo "#PBS -N ${BASE}_BR" >> ${OUT}
echo "#PBS -l walltime=12:00:00" >> ${OUT}
echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
echo "#PBS -q batch" >> ${OUT}
echo "#PBS -l mem=50gb" >> ${OUT}
echo "" >> ${OUT}
echo "cd ${output_directory}" >> ${OUT}
echo "module load ${GATK_module}" >> ${OUT}
echo "" >> ${OUT}
echo "time gatk BaseRecalibrator \
   -I ${output_directory}/${BASE}_piped.bam \
   --known-sites ${output_directory}/H0_variants_8SamplesNewRef.vcf \
   -O ${output_directory}/${BASE}_recal_dataNewRef.table \
   -R ${ref_genome}" >> ${OUT}
qsub ${OUT}

done

# ###################################################################################################
# ## Apply BQSR to bam files
# ###################################################################################################
#

# for file in ${raw_data}/${BASE}*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
#
#
# gatk ApplyBQSR \
#    -R ${ref_genome} \
#    -I ${raw_data}/${BASE}_pipedNewRef.bam \
#    -bqsr ${output_directory}/${BASE}_recal_dataNewRef.table \
#    -O ${output_directory}/${BASE}_recalibratedNewRef.bam
#
# done

# for file in ${output_directory}/${BASE}*_pipedNewRef.bam
#
# do
#
# FBASE=$(basename $file _pipedNewRef.bam)
# BASE=${FBASE%_pipedNewRef.bam}
# OUT="${BASE}_BQSR.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_BQSR" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${output_directory}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "gatk ApplyBQSR \
#    -R ${ref_genome} \
#    -I ${output_directory}/${BASE}_pipedNewRef.bam \
#    -bqsr ${output_directory}/${BASE}_recal_dataNewRef.table \
#    -O ${output_directory}/${BASE}_recalibratedNewRef.bam" >> ${OUT}
# qsub ${OUT}
#
# done


# ###################################################################################################
# ## Run HaplotypeCaller again on recalibrated samples
# ###################################################################################################
# ###################################################################################################
# #
# module load ${GATK_module}
#
### D1 samples
# for file in ${output_directory}/${BASE}*_recalibratedNewRef.bam
#
# do
#
# FBASE=$(basename $file _recalibratedNewRef.bam)
# BASE=${FBASE%_recalibratedNewRef.bam}
# OUT="${BASE}_recalHC2.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_recalHC2" >> ${OUT}
# echo "#PBS -l walltime=72:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=200gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
# -R ${ref_genome} \
# -ERC GVCF \
# -I ${output_directory}/${BASE}_recalibratedNewRef.bam \
# -ploidy 1 \
# -O ${output_directory}/${BASE}_variants.RecalNewRef.g.vcf" >> ${OUT}
# qsub ${OUT}
#
# done
#
# OUT="HM-H0-12_recalhapCall.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N HM-H0-12_recalHC" >> ${OUT}
# echo "#PBS -l walltime=30:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=100gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
# -R ${ref_genome} \
# -ERC GVCF \
# -I ${output_directory}/HM-H0-12_recalibrated.bam \
# -ploidy 1 \
# -O ${output_directory}/HM-H0-12_variants.Recal.g.vcf" >> ${OUT}
# qsub ${OUT}




###################################################################################################
## Combine gvcfs
###################################################################################################
###################################################################################################
#
#
# module load ${GATK_module}
#
# time gatk CombineGVCFs \
#    -R ${ref_genome} \
#    -O /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0_fullCohortNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-A_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-10_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-20_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-31_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-42_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-21_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-32_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-44_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-11_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-22_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-33_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-45_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-12_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-24_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-35_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-46_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-13_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-34_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-36_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-4_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-14_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-26_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-37_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-5_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-15_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-27_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-38_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-23_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-16_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-28_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-39_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-7_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-17_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-29_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-1_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-8_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-18_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-2_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-9_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-19_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-30_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-40_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-43_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-47_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-48_variants.RecalNewRef.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-41_variants.RecalNewRef.g.vcf
#
# # #
# #
# #    ###################################################################################################
# #    ## Genotype gVCFs (jointly)
# #    ###################################################################################################
# #    ###################################################################################################
# #    #
# #
#    time gatk GenotypeGVCFs \
#         -R ${ref_genome} \
#         -ploidy 1 \
#         --variant ${output_directory}/H0_fullCohort.g.vcf \
#         -O ${output_directory}/H0_fullCohort_int.vcf
#
# # ###################################################################################################
# # ## Genotype gVCFs (individually)
# # ###################################################################################################
# # ###################################################################################################
# # #
#
# # for file in ${output_directory}/${BASE}*_variants.Recal.g.vcf
# #
# # do
# #
# #   FBASE=$(basename $file _variants.Recal.g.vcf)
# #   BASE=${FBASE%_variants.Recal.g.vcf}
# #   OUT="${BASE}_recalHC.sh"
# #   echo "#!/bin/bash" >> ${OUT}
# #   echo "#PBS -N ${BASE}_recalHC" >> ${OUT}
# #   echo "#PBS -l walltime=12:00:00" >> ${OUT}
# #   echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# #   echo "#PBS -q batch" >> ${OUT}
# #   echo "#PBS -l mem=50gb" >> ${OUT}
# #   echo "" >> ${OUT}
# #   echo "module load ${GATK_module}" >> ${OUT}
# #   echo "" >> ${OUT}
# #   echo "time gatk GenotypeGVCFs \
# #      -R ${ref_genome} \
# #      --variant ${output_directory}/${BASE}_variants.Recal.g.vcf \
# #      -O ${output_directory}/${BASE}.recal.vcf" >> ${OUT}
# #   qsub ${OUT}
# #
# # done
# #
# # # ###################################################################################################
# # # ### Find coverage and put into 10k chunks
# # # ###################################################################################################
# #
# module load ${deeptools_module}
#
#
# for file in ${output_directory}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
#
# bamCoverage -b ${output_directory}/${BASE}_piped.bam -o ${output_directory}/${BASE}.bedgraph -of bedgraph -bs 10000
#
# done
#
# # for file in ${raw_data}/${BASE}*_piped.bam
# #
# # do
# #
# # FBASE=$(basename $file _piped.bam)
# # BASE=${FBASE%_piped.bam}
# # OUT="${BASE}_bamCoverage.sh"
# # echo "#!/bin/bash" >> ${OUT}
# # echo "#PBS -N ${BASE}_bamCoverage" >> ${OUT}
# # echo "#PBS -l walltime=12:00:00" >> ${OUT}
# # echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# # echo "#PBS -q batch" >> ${OUT}
# # echo "#PBS -l mem=20gb" >> ${OUT}
# # echo "" >> ${OUT}
# # echo "cd ${raw_data}" >> ${OUT}
# # echo "module load ${deeptools_module}" >> ${OUT}
# # echo "" >> ${OUT}
# # echo "bamCoverage -b ${raw_data}/${BASE}_piped.bam -o ${output_directory}/${BASE}.bedgraph -of bedgraph -bs 10000" >> ${OUT}
# # qsub ${OUT}
# #
# # done
#
#
# for file in ${output_directory}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
# samtools sort ${output_directory}/${BASE}_piped.bam \
# -o ${output_directory}/${BASE}.sorted.bam
#
# samtools depth \
# ${output_directory}/${BASE}.sorted.bam \
# |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${raw_data}/${BASE}.txt
#
# done
#
# # ###################################################################################################
# # ### Filter variants
# # Can easily run these interactively
# # ###################################################################################################
#
# # Get only those lines where there is actually a genotype call in the ancestor
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/H0_fullCohort_int.vcf \
# -O ${output_directory}/H0_FullCohort_AncCalls.vcf \
# -select 'vc.getGenotype("H0-A").isCalled()'
#
#
# # remove all lines in the ancestor that have a heterozygous genotype
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/H0_FullCohort_AncCalls.vcf \
# -O ${output_directory}/H0_FullCohort_AncCalls_NoHets.vcf \
# -select '!vc.getGenotype("H0-A").isHet()'
#
# # # filter out sites with low read depth
# gatk VariantFiltration \
#    -R ${ref_genome} \
#    -V ${output_directory}/H0_FullCohort_AncCalls_NoHets.vcf \
#    -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBias.vcf \
#    --set-filtered-genotype-to-no-call TRUE \
#    -G-filter "DP < 10"  -G-filter-name "depthGr10" \
#    -filter "MQ < 50.0" -filter-name "MQ50" \
#    -filter "SOR < 0.01" -filter-name "strandBias"
#
# #
# #   # remove filtered sites (these were set to no calls ./.)
#    gatk SelectVariants \
#    -R ${ref_genome} \
#    -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBias.vcf \
#    -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil.vcf \
#    --exclude-filtered TRUE
#
#    gatk SelectVariants \
#    -R ${ref_genome} \
#    -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil.vcf \
#    -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
#    -select 'vc.getGenotype("H0-A").isCalled()'
#
#    gatk SelectVariants \
#    -R ${ref_genome} \
#    -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
#    -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_vars.vcf \
#    -select-type SNP \
#    -select-type INDEL \
#    -select-type MIXED \
#    -select-type MNP \
#    -select-type SYMBOLIC
#
# # #
# # # #gives a final dataset with only called sites in the Ancestor, no heterozygous sites in the ancestor,
# # # # depth > 10, mapping quality > 50, and strand bias (SOR) > 0.01 (not significant)
# # #
# # # #Variants to table
# gatk VariantsToTable \
#      -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_vars.vcf \
#      -F CHROM -F POS -F REF -F ALT -F QUAL \
#      -GF AD -GF DP -GF GQ -GF GT \
#      -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_vars.txt
#
#      gatk VariantsToTable \
#           -V ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_vars.vcf \
#           -F CHROM -F POS -F REF -F ALT  \
#           -GF AD -GF GT \
#           -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_varsGT.txt
#
# # ##################################################################################################
# # #### UNIFIED GENOTYPER
# # ##################################################################################################
