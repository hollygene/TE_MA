#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N genotype
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=48:00:00
#PBS -l mem=150gb
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
#deeptools location
deeptools_module="deepTools/3.2.1-foss-2018a-Python-3.6.4"
#location of bamtoBigWig script and accessories
script_location="/scratch/hcm14449/TE_MA_Paradoxus/jbscripts"
#location of bam to bigwig script
bamToBigWig="/scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py"
#location of data to be analyzed
data_dir="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq"
#location of reference genome to be used
# ref_genome="/scratch/jc33471/pilon/337/annotation/genome.337.fasta"
ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa"
#directory reference genome is located in
# ref_genome_dir="/scratch/jc33471/pilon/337/"
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
# do_again="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/H0/do_again"
# genomicsdb_workspace_path="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB"
# sample_name_map="/home/hcm14449/Github/TE_MA/H0_sample_map.txt"
# tmp_DIR="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB/tmp"

# cd ${output_directory}
# rm *

# # module load ${picard_module}
# # module load ${bwa_module}
# module load ${samtools_module}
# # module load ${GATK_module}
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


# for file in ${raw_data}/*_R1_001.fastq
#
# do
#   FBASE=$(basename $file _R1_001.fastq)
#   BASE=${FBASE%_R1_001.fastq}
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
#       FASTQ=${raw_data}/${BASE}_R1_001.fastq \
#       FASTQ2=${raw_data}/${BASE}_R2_001.fastq  \
#       OUTPUT=${raw_data}/${BASE}_fastqtosam.bam \
#       READ_GROUP_NAME=${BASE} \
#       SAMPLE_NAME=${BASE} \
#       LIBRARY_NAME=D0 \
#       PLATFORM=illumina \
#       SEQUENCING_CENTER=GGBC" >> ${OUT}
# 	qsub ${OUT}
# done





# for file in ${raw_data}/*R1.fq
#
# do
#   FBASE=$(basename $file R1.fq)
#   BASE=${FBASE%R1.fq}
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
#       FASTQ=${raw_data}/${BASE}R1.fq \
#       FASTQ2=${raw_data}/${BASE}R1.fq  \
#       OUTPUT=${raw_data}/${BASE}_fastqtosam.bam \
#       READ_GROUP_NAME=${BASE} \
#       SAMPLE_NAME=${BASE} \
#       LIBRARY_NAME=D0 \
#       PLATFORM=illumina \
#       SEQUENCING_CENTER=GGBC" >> ${OUT}
# 	qsub ${OUT}
# done
# #######################################################################################
# # mark Illumina adapters
# #######################################################################################
#
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
#
# for file in ${raw_data}/*_fastqtosam.bam
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
# 	echo "cd ${raw_data}" >> ${OUT}
# 	echo "module load ${picard_module}" >> ${OUT}
#   echo "module load ${bwa_module}" >> ${OUT}
#   echo "module load ${samtools_module}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "mkdir ${raw_data}/TMP" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
#   I=${raw_data}/${BASE}_fastqtosam.bam \
#   O=${raw_data}/${BASE}_markilluminaadapters.bam \
#   M=${raw_data}/${BASE}_markilluminaadapters_metrics.txt \
#   TMP_DIR=${raw_data}/TMP" >> ${OUT}
# 	qsub ${OUT}
# done
# #######################################################################################
# # convert BAM to FASTQ and discount adapter sequences using SamToFastq
# #######################################################################################
#
# # for file in ${raw_data}/*_markilluminaadapters.bam
# #
# # do
# #
# # FBASE=$(basename $file _markilluminaadapters.bam)
# # BASE=${FBASE%_markilluminaadapters.bam}
# #
# # java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# # /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar SamToFastq \
# # I=${raw_data}/${BASE}_markilluminaadapters.bam \
# # FASTQ=${raw_data}/${BASE}_samtofastq_interleaved.fq \
# # CLIPPING_ATTRIBUTE=XT \
# # CLIPPING_ACTION=2 \
# # INTERLEAVE=true \
# # NON_PF=true \
# # TMP_DIR=${raw_data}/TMP
# #
# # done
#
# #######################################################################################
# # works: aligns samples to reference genome. Output is a .sam file
# #######################################################################################
#
# #  #index the ref genome
# # bwa index ${ref_genome}
# # #
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
# for file in ${raw_data}/*_markilluminaadapters.bam
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
# 	echo "cd ${raw_data}" >> ${OUT}
# 	echo "module load ${picard_module}" >> ${OUT}
#   echo "module load ${bwa_module}" >> ${OUT}
#   echo "module load ${samtools_module}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
#   I=${raw_data}/${BASE}_markilluminaadapters.bam \
#   FASTQ=${raw_data}/${BASE}_samtofastq_interleaved.fq \
#   CLIPPING_ATTRIBUTE=XT \
#   CLIPPING_ACTION=2 \
#   INTERLEAVE=true \
#   NON_PF=true \
#   TMP_DIR=${raw_data}/TMP" >> ${OUT}
# 	qsub ${OUT}
# done
## Piped command: SamToFastq, then bwa mem, then MergeBamAlignment
# for file in ${raw_data}/*_markilluminaadapters.bam
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
# echo "cd ${raw_data}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
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
# OUTPUT=${raw_data}/${BASE}_pipedNewRef.bam \
# R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
# CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
# INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
# PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
# TMP_DIR=${raw_data}/TMP" >> ${OUT}
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
# # samtools faidx ${ref_genome}
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
# module load ${GATK_module}
#
# ## H0 samples
# for file in /${raw_data}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
# OUT="${BASE}_HC.sh"
# echo "#!/bin/bash" > ${OUT}
# echo "#PBS -N ${BASE}_HC" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${raw_data}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "mkdir ${raw_data}/TMP" >> ${OUT}
# echo "time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -ERC GVCF \
#      -I /${raw_data}/${BASE}_piped.bam \
#      -ploidy 1 \
#      -O ${output_directory}/${BASE}_variants.g.vcf" >> ${OUT}
# qsub ${OUT}
# done
#
# # time gatk HaplotypeCaller \
# #      -R ${ref_genome} \
# #      -ERC GVCF \
# #      -I ${raw_data}/HM-H0-44_piped.bam \
# #      -ploidy 1 \
# #      -O ${output_directory}/HM-H0-44_variants.g.vcf
#
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
# # ###################################################################################################
# ### Combine gVCFs before joint genotyping
# # ###################################################################################################
#
#
# time gatk CombineGVCFs \
#  -O ${output_directory}/H0_smallCohort.g.vcf \
#  -R ${ref_genome} \
#  --variant ${output_directory}/H0-A_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-10_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-11_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-12_variants.g.vcf \
#  --variant ${output_directory}/H0-13_variants.g.vcf \
#  --variant ${output_directory}/H0-14_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-15_variants.g.vcf \
#  --variant ${output_directory}/HM-H0-16_variants.g.vcf
#
#
# # ###################################################################################################
# # ### Jointly genotype 8 random samples to identify consensus sequences
# # ###################################################################################################
#
# time gatk GenotypeGVCFs \
#         -R ${ref_genome} \
#         -ploidy 1 \
#         --variant ${output_directory}/H0_smallCohort.g.vcf \
#         -O ${output_directory}/H0_variants_7Samples.vcf
#
# # ###################################################################################################
# # ## Recalibrate base quality scores in all samples to mask any likely consensus variants
# # ###################################################################################################
#
# for file in ${raw_data}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
# OUT="${BASE}_BR.sh"
# echo "#!/bin/bash" > ${OUT}
# echo "#PBS -N ${BASE}_BR" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${raw_data}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk BaseRecalibrator \
#    -I ${raw_data}/${BASE}_piped.bam \
#    --known-sites ${output_directory}/H0_variants_8Samples.vcf \
#    -O ${output_directory}/${BASE}_recal_data.table \
#    -R ${ref_genome}" >> ${OUT}
# qsub ${OUT}
#
# done
#
# ###################################################################################################
# ## Apply BQSR to bam files
# ###################################################################################################
#
# for file in ${raw_data}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
# OUT="${BASE}_BQSR.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_BQSR" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${raw_data}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "gatk ApplyBQSR \
#    -R ${ref_genome} \
#    -I ${raw_data}/${BASE}_piped.bam \
#    -bqsr ${output_directory}/${BASE}_recal_data.table \
#    -O ${output_directory}/${BASE}_recalibrated.bam" >> ${OUT}
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
# for file in ${output_directory}/${BASE}*_recalibrated.bam
#
# do
#
# FBASE=$(basename $file _recalibrated.bam)
# BASE=${FBASE%_recalibrated.bam}
# OUT="${BASE}_recalHC.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_recalHC" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# echo "#PBS -q batch" >> ${OUT}
# echo "#PBS -l mem=50gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time gatk HaplotypeCaller \
# -R ${ref_genome} \
# -ERC GVCF \
# -I ${output_directory}/${BASE}_recalibrated.bam \
# -ploidy 1 \
# -O ${output_directory}/${BASE}_variants.Recal.g.vcf" >> ${OUT}
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
module load ${GATK_module}

# time gatk CombineGVCFs \
#    -R ${ref_genome} \
#    -O /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0_fullCohort.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-A_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-10_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-20_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-31_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-42_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-21_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-32_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-44_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-11_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-22_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-33_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-45_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-12_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-24_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-35_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-46_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-13_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-34_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-36_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-4_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-14_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-26_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-37_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-5_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-15_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-27_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-38_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-23_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-16_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-28_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-39_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0-7_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-17_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-29_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-1_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-8_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-18_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-2_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-9_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-19_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-30_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-40_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-43_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-47_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-48_variants.Recal.g.vcf \
#    -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-41_variants.Recal.g.vcf
#
# #
#
#    ###################################################################################################
#    ## Genotype gVCFs (jointly)
#    ###################################################################################################
#    ###################################################################################################
#    #
#
   time gatk GenotypeGVCFs \
        -R ${ref_genome} \
        -ploidy 1 \
        --variant ${output_directory}/H0_fullCohort.g.vcf \
        -O ${output_directory}/H0_fullCohort_int.vcf

# ###################################################################################################
# ## Genotype gVCFs (individually)
# ###################################################################################################
# ###################################################################################################
# #

# for file in ${output_directory}/${BASE}*_variants.Recal.g.vcf
#
# do
#
#   FBASE=$(basename $file _variants.Recal.g.vcf)
#   BASE=${FBASE%_variants.Recal.g.vcf}
#   OUT="${BASE}_recalHC.sh"
#   echo "#!/bin/bash" >> ${OUT}
#   echo "#PBS -N ${BASE}_recalHC" >> ${OUT}
#   echo "#PBS -l walltime=12:00:00" >> ${OUT}
#   echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
#   echo "#PBS -q batch" >> ${OUT}
#   echo "#PBS -l mem=50gb" >> ${OUT}
#   echo "" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
#   echo "" >> ${OUT}
#   echo "time gatk GenotypeGVCFs \
#      -R ${ref_genome} \
#      --variant ${output_directory}/${BASE}_variants.Recal.g.vcf \
#      -O ${output_directory}/${BASE}.recal.vcf" >> ${OUT}
#   qsub ${OUT}
#
# done
#
# # ###################################################################################################
# # ### Find coverage and put into 10k chunks
# # ###################################################################################################
#
# module load ${deeptools_module}
#
#
# for file in ${raw_data}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
#
# bamCoverage -b ${raw_data}/${BASE}_piped.bam -o ${output_directory}/${BASE}.bedgraph -of bedgraph -bs 10000
#
# done
# for file in ${raw_data}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
# OUT="${BASE}_bamCoverage.sh"
# echo "#!/bin/bash" >> ${OUT}
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

# ###################################################################################################
# ### Filter variants
# Can easily run these interactively
# ###################################################################################################

# Get only those lines where there is actually a genotype call in the ancestor
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/H0_fullCohort_int.vcf \
# -O ${output_directory}/H0_FullCohort_AncCalls.vcf \
# -select 'vc.getGenotype("H0-A").isCalled()'
#

# remove all lines in the ancestor that have a heterozygous genotype
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/H0_FullCohort_AncCalls.vcf \
# -O ${output_directory}/H0_FullCohort_AncCalls_NoHets.vcf \
# -select '!vc.getGenotype("H0-A").isHet()'
#
# # filter out sites with low read depth
# gatk VariantFiltration \
#    -R ${ref_genome} \
#    -V ${output_directory}/H0_FullCohort_AncCalls_NoHets.vcf \
#    -O ${output_directory}/H0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBias.vcf \
#    --set-filtered-genotype-to-no-call TRUE \
#    -G-filter "DP < 10"  -G-filter-name "depthGr10" \
#    -filter "MQ < 50.0" -filter-name "MQ50" \
#    -filter "SOR < 0.01" -filter-name "strandBias"
#
#
#   # remove filtered sites (these were set to no calls ./.)
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
# #
# # #gives a final dataset with only called sites in the Ancestor, no heterozygous sites in the ancestor,
# # # depth > 10, mapping quality > 50, and strand bias (SOR) > 0.01 (not significant)
# #
# # #Variants to table
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
# ##################################################################################################
# #### UNIFIED GENOTYPER
# ##################################################################################################
