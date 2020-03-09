#PBS -S /bin/bash
#PBS -q batch
#PBS -N test
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=4:00:00
#PBS -l mem=5gb
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
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/"
# mkdir $output_directory
#location of data to be used in the analysis
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas"
mcc_bams="/scratch/jc33471/paradoxusHolly/run0217all"
mcc_bam_indiv="/scratch/jc33471/paradoxusHolly/run0217all/out/Spar"
test_dir="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out"


# cd ${output_directory}
# rm *

module load ${picard_module}
module load ${bwa_module}
module load ${samtools_module}
module load ${GATK_module}

#######################################################################################
# create a uBAM file
#######################################################################################

# mcclintock bams

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
# O=${output_directory}/mcc_bams_out/${BASE}_unmapped.bam \
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

#######################################################################################
# mark Illumina adapters
#######################################################################################

# mv ${raw_data}/*_unmapped.bam ${output_directory}


for file in ${raw_data}/*_R1_001.fastq.gz

do
  FBASE=$(basename $file _R1_001.fastq.gz)
  BASE=${FBASE%_R1_001.fastq.gz}
	OUT="${BASE}_fqToBam.sh"
	echo "#!/bin/bash" > ${OUT}
	echo "#PBS -N ${BASE}_fqToBam" >> ${OUT}
	echo "#PBS -l walltime=72:00:00" >> ${OUT}
	echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
	echo "#PBS -q highmem_q" >> ${OUT}
	echo "#PBS -l mem=300gb" >> ${OUT}
	echo "" >> ${OUT}
	echo "cd ${output_directory}" >> ${OUT}
	echo "module load ${picard_module}" >> ${OUT}
  echo "module load ${bwa_module}" >> ${OUT}
  echo "module load ${samtools_module}" >> ${OUT}
  echo "module load ${GATK_module}" >> ${OUT}
	echo "" >> ${OUT}
  echo "mkdir ${output_directory}/TMP" >> ${OUT}
  echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
  FASTQ=${raw_data}/${BASE}_R1_001.fastq.gz \
  FASTQ2=${raw_data}/${BASE}_R2_001.fastq.gz  \
  OUTPUT=${raw_data}/${BASE}_unmapped.bam \
  READ_GROUP_NAME=${BASE} \
  SAMPLE_NAME=${BASE} \
  PLATFORM=illumina \
  SEQUENCING_CENTER=GGBC

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
  I=${output_directory}/${BASE}_unmapped.bam \
  O=${output_directory}/${BASE}_markilluminaadapters.bam \
  M=${output_directory}/${BASE}_markilluminaadapters_metrics.txt \
  TMP_DIR=${output_directory}/TMP

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
  I=${output_directory}/${BASE}_markilluminaadapters.bam \
  FASTQ=/dev/stdout \
  CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
  TMP_DIR=${output_directory}/TMP | \
  bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
  ALIGNED_BAM=/dev/stdin \
  UNMAPPED_BAM=${output_directory}/${BASE}_unmapped.bam \
  OUTPUT=${output_directory}/${BASE}_pipedNewRef.bam \
  R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
  CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
  INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
  PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
  TMP_DIR=${output_directory}/TMP

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
  INPUT=${output_directory}/${BASE}_pipedNewRef.bam \
  OUTPUT=${output_directory}/${BASE}_sorted.bam \
  SORT_ORDER=coordinate

  time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
  /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
  REMOVE_DUPLICATES=TRUE \
  I=${output_directory}/${BASE}_sorted.bam \
  O=${output_directory}/${BASE}_removedDuplicates.bam \
  M=${output_directory}/${BASE}_removedDupsMetrics.txt

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
  /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
  INPUT=${output_directory}/${BASE}_removedDuplicates.bam" >> ${OUT}

	qsub ${OUT}

done

#####################################################
# seq run 2

for file in ${raw_data}/*R1.fq.gz

do
  FBASE=$(basename $file R1.fq.gz)
  BASE=${FBASE%R1.fq.gz}
	OUT="${BASE}_fqToBam.sh"
	echo "#!/bin/bash" > ${OUT}
	echo "#PBS -N ${BASE}_fqToBam" >> ${OUT}
	echo "#PBS -l walltime=72:00:00" >> ${OUT}
	echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
	echo "#PBS -q highmem_q" >> ${OUT}
	echo "#PBS -l mem=300gb" >> ${OUT}
	echo "" >> ${OUT}
	echo "cd ${output_directory}" >> ${OUT}
	echo "module load ${picard_module}" >> ${OUT}
  echo "module load ${bwa_module}" >> ${OUT}
  echo "module load ${samtools_module}" >> ${OUT}
  echo "module load ${GATK_module}" >> ${OUT}
	echo "" >> ${OUT}
  echo "mkdir ${output_directory}/TMP" >> ${OUT}
  echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
  FASTQ=${raw_data}/${BASE}R1.fq.gz \
  FASTQ2=${raw_data}/${BASE}R2.fq.gz \
  OUTPUT=${output_directory}/${BASE}_unmapped.bam \
  READ_GROUP_NAME=${BASE} \
  SAMPLE_NAME=${BASE} \
  PLATFORM=illumina \
  SEQUENCING_CENTER=GGBC

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
  I=${output_directory}/${BASE}_unmapped.bam \
  O=${output_directory}/${BASE}_markilluminaadapters.bam \
  M=${output_directory}/${BASE}_markilluminaadapters_metrics.txt \
  TMP_DIR=${output_directory}/TMP

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
  I=${output_directory}/${BASE}_markilluminaadapters.bam \
  FASTQ=/dev/stdout \
  CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
  TMP_DIR=${output_directory}/TMP | \
  bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
  ALIGNED_BAM=/dev/stdin \
  UNMAPPED_BAM=${output_directory}/${BASE}_unmapped.bam \
  OUTPUT=${output_directory}/${BASE}_pipedNewRef.bam \
  R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
  CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
  INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
  PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
  TMP_DIR=${output_directory}/TMP

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
  INPUT=${output_directory}/${BASE}_pipedNewRef.bam \
  OUTPUT=${output_directory}/${BASE}_sorted.bam \
  SORT_ORDER=coordinate

  time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
  /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
  REMOVE_DUPLICATES=TRUE \
  I=${output_directory}/${BASE}_sorted.bam \
  O=${output_directory}/${BASE}_removedDuplicates.bam \
  M=${output_directory}/${BASE}_removedDupsMetrics.txt

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
  /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
  INPUT=${output_directory}/${BASE}_removedDuplicates.bam" >> ${OUT}

	qsub ${OUT}

done
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
# for file in ${output_directory}/*_markilluminaadapters.bam
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
# echo "cd ${output_directory}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
#   I=${output_directory}/${BASE}_markilluminaadapters.bam \
#   FASTQ=/dev/stdout \
#   CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
#   TMP_DIR=${output_directory}/TMP | \
#   bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
#   java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
#   ALIGNED_BAM=/dev/stdin \
#   UNMAPPED_BAM=${output_directory}/${BASE}_fastqtosam.bam \
#   OUTPUT=${output_directory}/${BASE}_pipedNewRef.bam \
#   R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
#   CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
#   INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
#   PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
#   TMP_DIR=${output_directory}/TMP" >> ${OUT}
# qsub ${OUT}
#
# done
## #######################################################################################

# Sort the piped command output
## #######################################################################################

# for file in ${output_directory}/*_pipedNewRef.bam
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
# echo "cd ${output_directory}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
#    /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
#     INPUT=${output_directory}/${BASE}_pipedNewRef.bam \
#     OUTPUT=${output_directory}/${BASE}_sorted.bam \
#     SORT_ORDER=coordinate" >> ${OUT}
#
# qsub ${OUT}
#
# done

###################################################################################################
## Picard to mark and remove duplicates
###################################################################################################

# for file in ${output_directory}/${BASE}*_sorted.bam
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
# echo "cd ${output_directory}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
# /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/${BASE}_sorted.bam \
# O=${output_directory}/${BASE}_removedDuplicates.bam \
# M=${output_directory}/${BASE}_removedDupsMetrics.txt" >> ${OUT}
# qsub ${OUT}
#
# done

##############################################################################################
############## ** Picard to BuildBamIndex              #################################
#################################################################################################

# for file in ${output_directory}/*_removedDuplicates.bam
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
# echo "cd ${output_directory}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
#    /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
#     INPUT=${output_directory}/${BASE}_removedDuplicates.bam " >> ${OUT}
#
# qsub ${OUT}
#
# done




###################################################################################################
# Using GATK HaplotypeCaller in GVCF mode
# apply appropriate ploidy for each sample
# will need to do this separtely for haploid and diploid samples
###################################################################################################
#
# module load ${GATK_module}

## D0 samples
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
#      -ploidy 2 \
#      -O ${output_directory}/${BASE}_variants.g.vcf
#
# done

# for file in ${raw_data}/${BASE}*_piped.bam
#
# do
#
# FBASE=$(basename $file _piped.bam)
# BASE=${FBASE%_piped.bam}
# 	OUT="${BASE}_HaplotypeCaller.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_HaplotypeCaller" >> ${OUT}
# 	echo "#PBS -l walltime=12:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
# 	echo "#PBS -q batch" >> ${OUT}
# 	echo "#PBS -l mem=50gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${raw_data}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "time gatk HaplotypeCaller \
#        -R ${ref_genome} \
#        -ERC GVCF \
#        -I ${raw_data}/${BASE}_piped.bam \
#        -ploidy 2 \
#        -O ${output_directory}/${BASE}_variants.g.vcf" >> ${OUT}
# 	qsub ${OUT}
# done

# for file in ${raw_data}/${BASE}*_piped.bam
#
# do
#   FBASE=$(basename $file _piped.bam)
#   BASE=${FBASE%_piped.bam}
# 	OUT="${BASE}_HaplotypeCaller.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_HaplotypeCaller" >> ${OUT}
# 	echo "#PBS -l walltime=72:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# 	echo "#PBS -q highmem_q" >> ${OUT}
# 	echo "#PBS -l mem=200gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${raw_data}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "time gatk HaplotypeCaller \
#        -R ${ref_genome} \
#        -ERC GVCF \
#        -I ${raw_data}/${BASE}_piped.bam \
#        -ploidy 2 \
#        -O ${output_directory}/${BASE}_variants.g.vcf" >> ${OUT}
# 	qsub ${OUT}
# done

# ###################################################################################################
### Combine gVCFs before joint genotyping
# ###################################################################################################


# time gatk CombineGVCFs \
#  -O ${output_directory}/D0_cohortNewRef.g.vcf \
#  -R ${ref_genome} \
#  --variant ${output_directory}/HM-D0-A_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-D0-10_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-D0-11_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-D0-12_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-D0-13_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-D0-14_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-D0-15_variantsNewRef.g.vcf \
#  --variant ${output_directory}/HM-D0-16_variantsNewRef.g.vcf


###################################################################################################
### Jointly genotype 8 random samples to identify consensus sequences
###################################################################################################

# time gatk GenotypeGVCFs \
#         -R ${ref_genome} \
#         --variant ${output_directory}/D0_cohortNewRef.g.vcf \
#         -O ${output_directory}/D0_variants_8SamplesNewRef.vcf


# ###################################################################################################
# ## Recalibrate base quality scores in all samples to mask any likely consensus variants
# ###################################################################################################
#
# for file in ${raw_data}/${BASE}*_piped.bam
#
#         do
#
#         FBASE=$(basename $file _piped.bam)
#         BASE=${FBASE%_piped.bam}
#
#
#         time gatk BaseRecalibrator \
#            -I ${raw_data}/${BASE}_piped.bam \
#            --known-sites ${output_directory}/D0_variants_8Samples.vcf \
#            -O ${output_directory}/${BASE}_recal_data.table \
#            -R ${ref_genome}
#
#         done

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
#    --known-sites ${output_directory}/D0_variants_8SamplesNewRef.vcf \
#    -O ${output_directory}/${BASE}_recal_dataNewRef.table \
#    -R ${ref_genome}" >> ${OUT}
# qsub ${OUT}
#
# done
# ###################################################################################################
  # ## Apply BQSR to bam files
  # ###################################################################################################
#
# for file in ${raw_data}/${BASE}*_pipedNewRef.bam
#
#       do
#         FBASE=$(basename $file _pipedNewRef.bam)
#         BASE=${FBASE%_pipedNewRef.bam}
#
#
#         gatk ApplyBQSR \
#            -R ${ref_genome} \
#            -I ${raw_data}/${BASE}_pipedNewRef.bam \
#            -bqsr ${output_directory}/${BASE}_recal_dataNewRef.table \
#            -O ${output_directory}/${BASE}_recalibratedNewRef.bam
#
#         done

  # ###################################################################################################
  ### Run HaplotypeCaller again on recalibrated samples
  # ###################################################################################################
  # ###################################################################################################
  # #
        # module load ${GATK_module}
        #
        # # D1 samples
        # for file in ${output_directory}/do_again/${BASE}*_recalibrated.bam
        #
        # do
        #
        # FBASE=$(basename $file _recalibrated.bam)
        # BASE=${FBASE%_recalibrated.bam}
        #
        # time gatk HaplotypeCaller \
        # -R ${ref_genome} \
        # -ERC GVCF \
        # -I ${output_directory}/do_again/${BASE}_recalibrated.bam \
        # -ploidy 2 \
        # -O ${output_directory}/do_again/${BASE}_variants.Recal.g.vcf
        # done
        #
        # for file in ${output_directory}/${BASE}*_recalibratedNewRef.bam
        #
        # do
        #   FBASE=$(basename $file _recalibratedNewRef.bam)
        #   BASE=${FBASE%_recalibratedNewRef.bam}
        # 	OUT="${BASE}_HaplotypeCaller.sh"
        # 	echo "#!/bin/bash" > ${OUT}
        # 	echo "#PBS -N ${BASE}_HaplotypeCaller" >> ${OUT}
        # 	echo "#PBS -l walltime=12:00:00" >> ${OUT}
        # 	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
        # 	echo "#PBS -q batch" >> ${OUT}
        # 	echo "#PBS -l mem=10gb" >> ${OUT}
        # 	echo "" >> ${OUT}
        # 	echo "cd ${raw_data}" >> ${OUT}
        #   echo "module load ${GATK_module}" >> ${OUT}
        # 	echo "" >> ${OUT}
        #   echo "time gatk HaplotypeCaller \
        #   -R ${ref_genome} \
        #   -ERC GVCF \
        #   -I ${output_directory}/${BASE}_recalibratedNewRef.bam \
        #   -ploidy 2 \
        #   -O ${output_directory}/${BASE}_variants.RecalNewRef.g.vcf" >> ${OUT}
        # 	qsub ${OUT}
        # done

#
#           ###################################################################################################
#           ## Combine gvcfs
#           ###################################################################################################
#           ###################################################################################################
#           #
#           #
  # module load ${GATK_module}
#
#           time gatk CombineGVCFs \
#              -R ${ref_genome} \
#              -O /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/D0_FullCohort.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-A_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-10_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-20_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-31_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-42_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-21_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-32_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-44_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-11_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-22_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-33_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-12_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-24_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-35_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-46_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-13_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-34_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-36_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-4_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-14_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-26_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-37_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-5_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-15_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-27_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-38_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-23_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-16_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-3_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-39_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-7_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-17_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-25_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-1_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-6_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-18_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-2_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-9_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-19_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-30_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-40_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-43_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-47_variants.Recal.g.vcf \
#              -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-48_variants.Recal.g.vcf
#
# #
#              ###################################################################################################
#              ## Genotype gVCFs (jointly)
#              ###################################################################################################
#              ###################################################################################################
#
#
#              time gatk GenotypeGVCFs \
#                   -R ${ref_genome} \
#                   -ploidy 2 \
#                   --variant ${output_directory}/D0_FullCohort.g.vcf \
#                   -O ${output_directory}/D0_FullCohort.vcf
#
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
# # Get only those lines where there is actually a genotype call in the ancestor
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/D0_FullCohort.vcf \
# -O ${output_directory}/D0_FullCohort_AncCalls.vcf \
# -select 'vc.getGenotype("HM-D0-A").isCalled()'
#
#
# # remove all lines in the ancestor that have a heterozygous genotype
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/D0_FullCohort_AncCalls.vcf \
# -O ${output_directory}/D0_FullCohort_AncCalls_NoHets.vcf \
# -select '!vc.getGenotype("HM-D0-A").isHet()'
#
# # filter out sites with low read depth
# gatk VariantFiltration \
#    -R ${ref_genome} \
#    -V ${output_directory}/D0_FullCohort_AncCalls_NoHets.vcf \
#    -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBias.vcf \
#    --set-filtered-genotype-to-no-call TRUE \
#    -G-filter "DP < 10"  -G-filter-name "depthGr10" \
#    -filter "MQ < 50.0" -filter-name "MQ50" \
#    -filter "SOR < 0.01" -filter-name "strandBias"
#
#
#   # remove filtered sites (these were set to no calls ./.)
#    gatk SelectVariants \
#    -R ${ref_genome} \
#    -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBias.vcf \
#    -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil.vcf \
#    --exclude-filtered TRUE
#
#    gatk SelectVariants \
#    -R ${ref_genome} \
#    -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil.vcf \
#    -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
#    -select 'vc.getGenotype("HM-D0-A").isCalled()'
#
# cd ${output_directory}
#
#    gatk SelectVariants \
#    -R ${ref_genome} \
#    -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
#    -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_SNPs.vcf \
#    --max-nocall-fraction 0 \
#    --exclude-non-variants TRUE \
#    -select-type SNP
#
#    gatk SelectVariants \
#    -R ${ref_genome} \
#    -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
#    -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_Indels.vcf \
#       --max-nocall-fraction 0.001 \
#       -select-type INDEL
#
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
# gatk VariantsToTable \
#      -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_vars.vcf \
#      -F CHROM -F POS -F REF -F ALT -F QUAL \
#      -GF AD -GF DP -GF GQ -GF GT \
#      -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_vars.txt
#
#      gatk VariantsToTable \
#           -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_SNPs.vcf \
#           -F CHROM -F POS -F REF -F ALT  \
#           -GF AD -GF GT \
#           -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_SNPs.txt
#
#           gatk VariantsToTable \
#                -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_Indels.vcf \
#                -F CHROM -F POS -F REF -F ALT  \
#                -GF AD -GF GT \
#                -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_Indels.txt

# # ###################################################################################################
### This stuff is in Excel
# # ###################################################################################################
#
# 3 columns after last genotype: # samples (count # samples you have), #same (meaning number of samples that have the same genotype as the ancestor),
# and # variant (# samples different genotype than the ancestor).
# Codes are:
#same: =COUNTIFS(E2:CL2,CN2)
#variant: =CO2-CP2 ( just subtract the # same from the # samples)
