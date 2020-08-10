#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N D0_pipeline
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=96:00:00
#PBS -l mem=250gb
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
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/gVCFs"
# mkdir $output_directory
#location of data to be used in the analysis
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas"
mcc_bams="/scratch/jc33471/paradoxusHolly/run0217all"
mcc_bam_indiv="/scratch/jc33471/paradoxusHolly/run0217all/out/Spar"
test_dir="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out"


cd ${output_directory}
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

expr $(samtools view -f 4 H0-A_.sam -c) + $(samtools view -F 2308 H0-A_.sam -c)

echo $(zcat D0-A_R2.fq.gz|wc -l)/4|bc

samtools view -F260 H0-A_.sam | cut -f3 | datamash -s -g 1 count 1 > counts.txt

grep -c "^>" genome.337.fasta


#######################################################################################
# mark Illumina adapters
#######################################################################################

# mv ${raw_data}/*_unmapped.bam ${output_directory}


for file in ${raw_data}/*_R1_001.fastq.gz

do
  FBASE=$(basename $file _R1_001.fastq.gz)
  BASE=${FBASE%_R1_001.fastq.gz}
	OUT="${BASE}_BWA.sh"
	echo "#!/bin/bash" > ${OUT}
	echo "#PBS -N ${BASE}_BWA" >> ${OUT}
	echo "#PBS -l walltime=12:00:00" >> ${OUT}
	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
	echo "#PBS -q batch" >> ${OUT}
	echo "#PBS -l mem=30gb" >> ${OUT}
	echo "" >> ${OUT}
	echo "cd ${output_directory}/D0" >> ${OUT}
  echo "module load ${bwa_module}" >> ${OUT}
	echo "" >> ${OUT}
  echo "mkdir ${output_directory}/D0/TMP" >> ${OUT}
	echo "mkdir ${output_directory}/D0/fromFASTQ" >> ${OUT}
	echo "cd ${output_directory}/D0/fromFASTQ" >> ${OUT}
  echo "bwa mem -M -t 7 -p ${ref_genome} ${raw_data}/${BASE}_R1_001.fastq.gz > ${output_directory}/D0/fromFASTQ/${BASE}.sam" >> ${OUT}
	qsub ${OUT}

done

#####################################################
# seq run 2

for file in ${raw_data}/*R1.fq.gz

do
  FBASE=$(basename $file R1.fq.gz)
  BASE=${FBASE%R1.fq.gz}
	OUT="${BASE}_BWA.sh"
	echo "#!/bin/bash" > ${OUT}
	echo "#PBS -N ${BASE}_BWA" >> ${OUT}
	echo "#PBS -l walltime=12:00:00" >> ${OUT}
	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
	echo "#PBS -q batch" >> ${OUT}
	echo "#PBS -l mem=30gb" >> ${OUT}
	echo "" >> ${OUT}
	echo "cd ${output_directory}" >> ${OUT}
	echo "module load ${picard_module}" >> ${OUT}
  echo "module load ${bwa_module}" >> ${OUT}
  echo "module load ${samtools_module}" >> ${OUT}
  echo "module load ${GATK_module}" >> ${OUT}
	echo "" >> ${OUT}
	echo "mkdir ${output_directory}/D0/TMP" >> ${OUT}
	echo "mkdir ${output_directory}/D0/fromFASTQ" >> ${OUT}
	echo "cd ${output_directory}/D0/fromFASTQ" >> ${OUT}
  echo "bwa mem -M -t 7 -p ${ref_genome} ${raw_data}/${BASE}R1.fq.gz > ${output_directory}/D0/fromFASTQ/${BASE}.sam" >> ${OUT}

	qsub ${OUT}

done



for file in ${mcc_bam_indiv}/*_val/bam/*_val.bam;

do

FBASE=$(basename $file _val.bam)
BASE=${FBASE%_val.bam}
	OUT="${BASE}_fqToIndexed.sh"
	echo "#!/bin/bash" > ${OUT}
	echo "#PBS -N ${BASE}_fqToIndexed" >> ${OUT}
	echo "#PBS -l walltime=72:00:00" >> ${OUT}
	echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
	echo "#PBS -q highmem_q" >> ${OUT}
	echo "#PBS -l mem=300gb" >> ${OUT}
	echo "" >> ${OUT}
	echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
	echo "module load ${picard_module}" >> ${OUT}
  echo "module load ${bwa_module}" >> ${OUT}
  echo "module load ${samtools_module}" >> ${OUT}
  echo "module load ${GATK_module}" >> ${OUT}
	echo "" >> ${OUT}
  echo "mkdir ${output_directory}/mcc_bams_out/TMP" >> ${OUT}
  echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
  /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar RevertSam \
  I=${mcc_bam_indiv}/${BASE}_val/bam/${BASE}_val.bam \
  O=${output_directory}/mcc_bams_out/${BASE}_unmapped.bam \
  SANITIZE=true \
  MAX_DISCARD_FRACTION=0.005 \
  ATTRIBUTE_TO_CLEAR=XT \
  ATTRIBUTE_TO_CLEAR=XN \
  ATTRIBUTE_TO_CLEAR=AS \
  ATTRIBUTE_TO_CLEAR=OC \
  ATTRIBUTE_TO_CLEAR=OP \
  SORT_ORDER=queryname \
  RESTORE_ORIGINAL_QUALITIES=true \
  REMOVE_DUPLICATE_INFORMATION=true \
  REMOVE_ALIGNMENT_INFORMATION=true" >> ${OUT}
	qsub ${OUT}

done


######################################################################################
##################### ############ ########### ########### ########### ########## ####

for file in ${output_directory}/mcc_bams_out/*_unmapped.bam;
do

FBASE=$(basename $file _unmapped.bam)
BASE=${FBASE%_unmapped.bam}
	OUT="${BASE}_mark.sh"
	echo "#!/bin/bash" > ${OUT}
	echo "#PBS -N ${BASE}_mark" >> ${OUT}
	echo "#PBS -l walltime=72:00:00" >> ${OUT}
	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
	echo "#PBS -q batch" >> ${OUT}
	echo "#PBS -l mem=50gb" >> ${OUT}
	echo "" >> ${OUT}
	echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
	echo "module load ${picard_module}" >> ${OUT}
  echo "module load ${bwa_module}" >> ${OUT}
  echo "module load ${samtools_module}" >> ${OUT}
  echo "module load ${GATK_module}" >> ${OUT}
	echo "" >> ${OUT}
  echo "mkdir ${output_directory}/mcc_bams_out/TMP" >> ${OUT}
  echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
  I=${output_directory}/mcc_bams_out/${BASE}_unmapped.bam \
  O=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters.bam \
  M=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters_metrics.txt \
  TMP_DIR=${output_directory}/mcc_bams_out/TMP" >> ${OUT}
	qsub ${OUT}

done


for file in ${output_directory}/mcc_bams_out/*_markilluminaadapters.bam;
do

FBASE=$(basename $file _markilluminaadapters.bam)
BASE=${FBASE%_markilluminaadapters.bam}
	OUT="${BASE}_piped.sh"
	echo "#!/bin/bash" > ${OUT}
	echo "#PBS -N ${BASE}_piped" >> ${OUT}
	echo "#PBS -l walltime=72:00:00" >> ${OUT}
	echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
	echo "#PBS -q batch" >> ${OUT}
	echo "#PBS -l mem=50gb" >> ${OUT}
	echo "" >> ${OUT}
	echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
	echo "module load ${picard_module}" >> ${OUT}
  echo "module load ${bwa_module}" >> ${OUT}
  echo "module load ${samtools_module}" >> ${OUT}
  echo "module load ${GATK_module}" >> ${OUT}
	echo "" >> ${OUT}
  echo "mkdir ${output_directory}/mcc_bams_out/TMP" >> ${OUT}
  echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
	  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
	  I=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters.bam \
	  FASTQ=/dev/stdout \
	  CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
	  TMP_DIR=${output_directory}/mcc_bams_out/TMP | \
	  bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
	  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
	  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
	  ALIGNED_BAM=/dev/stdin \
	  UNMAPPED_BAM=${output_directory}/mcc_bams_out/${BASE}_unmapped.bam \
	  OUTPUT=${output_directory}/mcc_bams_out/${BASE}_pipedNewRef.bam \
	  R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
	  CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
	  INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
	  PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
	  TMP_DIR=${output_directory}/mcc_bams_out/TMP" >> ${OUT}
	qsub ${OUT}

done

  # java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  # /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
  # I=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters.bam \
  # FASTQ=/dev/stdout \
  # CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
  # TMP_DIR=${output_directory}/mcc_bams_out/TMP | \
  # bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
  # java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  # /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
  # ALIGNED_BAM=/dev/stdin \
  # UNMAPPED_BAM=${output_directory}/mcc_bams_out/${BASE}_unmapped.bam \
  # OUTPUT=${output_directory}/mcc_bams_out/${BASE}_pipedNewRef.bam \
  # R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
  # CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
  # INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
  # PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
  # TMP_DIR=${output_directory}/mcc_bams_out/TMP
	#
  # java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  # /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
  # INPUT=${output_directory}/mcc_bams_out/${BASE}_pipedNewRef.bam \
  # OUTPUT=${output_directory}/mcc_bams_out/${BASE}_sorted.bam \
  # SORT_ORDER=coordinate
	#
  # time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
  # /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
  # REMOVE_DUPLICATES=TRUE \
  # I=${output_directory}/mcc_bams_out/${BASE}_sorted.bam \
  # O=${output_directory}/mcc_bams_out/${BASE}_removedDuplicates.bam \
  # M=${output_directory}/mcc_bams_out/${BASE}_removedDupsMetrics.txt
	#
  # java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
  # /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
  # INPUT=${output_directory}/mcc_bams_out/${BASE}_removedDuplicates.bam" >> ${OUT}
	#



####################################################
## MCC Bams

for file in ${mcc_bam_indiv}/*_val/bam/*_val.bam;

do

FBASE=$(basename $file _val.bam)
BASE=${FBASE%_val.bam}
	OUT="${BASE}_fqToIndexed.sh"
	echo "#!/bin/bash" > ${OUT}
	echo "#PBS -N ${BASE}_fqToIndexed" >> ${OUT}
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
	echo "mkdir ${output_directory}/mcc_bams" >> ${OUT}
  echo "mkdir ${output_directory}/mcc_bams/TMP" >> ${OUT}
  echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
	/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
	INPUT=${mcc_bam_indiv}/${BASE}_val/bam/${BASE}_val.bam \
	OUTPUT=${output_directory}/mcc_bams_out/${BASE}_sorted.bam \
	SORT_ORDER=queryname

	java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
  I=${output_directory}/mcc_bams_out/${BASE}_sorted.bam \
  O=${output_directory}/mcc_bams/${BASE}_markilluminaadapters.bam \
  M=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters_metrics.txt \
  TMP_DIR=${output_directory}/mcc_bams_out/TMP

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
  /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
  INPUT=${output_directory}/mcc_bams_out/${BASE}_markilluminaadapters.bam \
  OUTPUT=${output_directory}/mcc_bams_out/${BASE}_sorted2.bam \
  SORT_ORDER=coordinate

  time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
  /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
  REMOVE_DUPLICATES=TRUE \
  I=${output_directory}/mcc_bams_out/${BASE}_sorted2.bam \
  O=${output_directory}/mcc_bams_out/${BASE}_removedDuplicates.bam \
  M=${output_directory}/mcc_bams_out/${BASE}_removedDupsMetrics.txt

  java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar" -jar  \
  /usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar BuildBamIndex \
  INPUT=${output_directory}/mcc_bams_out/${BASE}_removedDuplicates.bam" >> ${OUT}
	qsub ${OUT}
done





## Running some things interactively
${output_directory}/mcc_bams

module load ${picard_module}
module load ${bwa_module}
module load ${samtools_module}
module load ${GATK_module}

java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
INPUT=${mcc_bam_indiv}/${BASE}_val/bam/${BASE}_val.bam \
OUTPUT=${output_directory}/mcc_bams_out/${BASE}_sorted.bam \
SORT_ORDER=queryname

# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar FixMateInformation \
# I=${output_directory}/mcc_bams_out/HM-H0-A_R1_001_sorted.bam \
# O=${output_directory}/mcc_bams_out/HM-H0-A_R1_001_sortedFixedMate.bam \
# ADD_MATE_CIGAR=TRUE

java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MarkIlluminaAdapters \
I=${output_directory}/mcc_bams_out/HM-H0-10_R1_001_sorted.bam \
O=${output_directory}/mcc_bams_out/HM-H0-10_R1_001_markilluminaadapters.bam \
M=${output_directory}/mcc_bams_out/HM-H0-10_R1_001_markilluminaadapters_metrics.txt \
TMP_DIR=${output_directory}/mcc_bams_out/TMP






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

for file in ${output_directory}/mcc_bams_out/*_pipedNewRef.bam

do

FBASE=$(basename $file _pipedNewRef.bam)
BASE=${FBASE%_pipedNewRef.bam}
OUT="${BASE}_sortSam.sh"
echo "#!/bin/bash" >> ${OUT}
echo "#PBS -N ${BASE}_sortSam" >> ${OUT}
echo "#PBS -l walltime=12:00:00" >> ${OUT}
echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
echo "#PBS -q batch" >> ${OUT}
echo "#PBS -l mem=50gb" >> ${OUT}
echo "" >> ${OUT}
echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
echo "module load ${picard_module}" >> ${OUT}
echo "module load ${bwa_module}" >> ${OUT}
echo "module load ${samtools_module}" >> ${OUT}
echo "module load ${GATK_module}" >> ${OUT}
echo "" >> ${OUT}
echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
   /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SortSam \
    INPUT=${output_directory}/mcc_bams_out/${BASE}_pipedNewRef.bam \
    OUTPUT=${output_directory}/mcc_bams_out/${BASE}_sorted.bam \
    SORT_ORDER=coordinate" >> ${OUT}

qsub ${OUT}

done

###################################################################################################
## Picard to mark and remove duplicates
###################################################################################################

for file in ${output_directory}/mcc_bams_out/${BASE}*_sorted.bam

do

FBASE=$(basename $file _sorted.bam)
BASE=${FBASE%_sorted.bam}
OUT="${BASE}_sorted.sh"
echo "#!/bin/bash" >> ${OUT}
echo "#PBS -N ${BASE}_sorted" >> ${OUT}
echo "#PBS -l walltime=72:00:00" >> ${OUT}
echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
echo "#PBS -q highmem_q" >> ${OUT}
echo "#PBS -l mem=200gb" >> ${OUT}
echo "" >> ${OUT}
echo "cd ${output_directory}/mcc_bams_out" >> ${OUT}
echo "module load ${GATK_module}" >> ${OUT}
echo "module load ${picard_module}" >> ${OUT}
echo "" >> ${OUT}
echo "time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.21.6-Java-11" -jar  \
/usr/local/apps/eb/picard/2.21.6-Java-11/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
I=${output_directory}/mcc_bams_out/${BASE}_sorted.bam \
O=${output_directory}/mcc_bams_out/${BASE}_removedDuplicates.bam \
M=${output_directory}/mcc_bams_out/${BASE}_removedDupsMetrics.txt" >> ${OUT}
qsub ${OUT}

done

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
### Diploids

# for file in ${output_directory}/D0/${BASE}*_removedDuplicates.bam
#
# do
#
#   FBASE=$(basename $file _removedDuplicates.bam)
#   BASE=${FBASE%_removedDuplicates.bam}
# 	OUT="${BASE}_HaplotypeCaller.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_HaplotypeCaller" >> ${OUT}
# 	echo "#PBS -l walltime=72:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# 	echo "#PBS -q highmem_q" >> ${OUT}
# 	echo "#PBS -l mem=200gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "time gatk HaplotypeCaller \
#        -R ${ref_genome} \
#        -ERC GVCF \
#        -I ${output_directory}/D0/${BASE}_removedDuplicates.bam \
#        -ploidy 2 \
#        -O ${output_directory}/D0/${BASE}_variants.g.vcf" >> ${OUT}
#
# 	qsub ${OUT}
#
# done


###################################################################################################
#


# for file in ${output_directory}/D1/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# 	OUT="${BASE}_HaplotypeCaller.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_HaplotypeCaller" >> ${OUT}
# 	echo "#PBS -l walltime=72:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# 	echo "#PBS -q highmem_q" >> ${OUT}
# 	echo "#PBS -l mem=200gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}/D1" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "time gatk HaplotypeCaller \
#        -R ${ref_genome} \
#        -ERC GVCF \
#        -I ${output_directory}/D1/${BASE}_removedDuplicates.bam \
#        -ploidy 2 \
#        -O ${output_directory}/D1/${BASE}_variants.g.vcf" >> ${OUT}
# 	qsub ${OUT}
# done

###################################################################################################
#

# for file in ${output_directory}/D20/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# 	OUT="${BASE}_HaplotypeCaller.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_HaplotypeCaller" >> ${OUT}
# 	echo "#PBS -l walltime=72:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# 	echo "#PBS -q highmem_q" >> ${OUT}
# 	echo "#PBS -l mem=200gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}/D20" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "time gatk HaplotypeCaller \
#        -R ${ref_genome} \
#        -ERC GVCF \
#        -I ${output_directory}/D20/${BASE}_removedDuplicates.bam \
#        -ploidy 2 \
#        -O ${output_directory}/D20/${BASE}_variants.g.vcf" >> ${OUT}
# 	qsub ${OUT}
# done

###################################################################################################
#
### Haploids

# for file in ${output_directory}/H0/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
# 	OUT="${BASE}_HaplotypeCaller.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_HaplotypeCaller" >> ${OUT}
# 	echo "#PBS -l walltime=72:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# 	echo "#PBS -q highmem_q" >> ${OUT}
# 	echo "#PBS -l mem=200gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}/H0" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "time gatk HaplotypeCaller \
#        -R ${ref_genome} \
#        -ERC GVCF \
#        -I ${output_directory}/H0/${BASE}_removedDuplicates.bam \
#        -ploidy 1 \
#        -O ${output_directory}/H0/${BASE}_variants.g.vcf" >> ${OUT}
# 	qsub ${OUT}
# done

###################################################################################################
#
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
#  --variant ${output_directory}/D0-A__variants.g.vcf \
#  --variant ${output_directory}/HM-D0-10_variants.g.vcf \
#  --variant ${output_directory}/HM-D0-11_variants.g.vcf \
#  --variant ${output_directory}/HM-D0-12_variants.g.vcf \
#  --variant ${output_directory}/HM-D0-13_variants.g.vcf \
#  --variant ${output_directory}/HM-D0-14_variants.g.vcf \
#  --variant ${output_directory}/HM-D0-15_variants.g.vcf \
#  --variant ${output_directory}/HM-D0-16_variants.g.vcf




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
# for file in ${output_directory}/${BASE}*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
#
# time gatk BaseRecalibrator \
# -I ${output_directory}/${BASE}_removedDuplicates.bam \
# --known-sites ${output_directory}/D0_variants_8SamplesNewRef.vcf \
# -O ${output_directory}/${BASE}_recal_data.table \
# -R ${ref_genome}
#
# done


# ###################################################################################################
# ## Apply BQSR to bam files
# ###################################################################################################
#
# module load ${GATK_module}

# for file in ${output_directory}/${BASE}*_removedDuplicates.bam
#
# do
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
#
#
# gatk ApplyBQSR \
# -R ${ref_genome} \
# -I ${output_directory}/${BASE}_removedDuplicates.bam \
# -bqsr ${output_directory}/${BASE}_recal_data.table \
# -O ${output_directory}/${BASE}_recalibratedNewRef.bam
#
# done


  # ###################################################################################################
  ### Run HaplotypeCaller again on recalibrated samples
  # ###################################################################################################
  # ###################################################################################################
  # #


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
# -ploidy 2 \
# -O ${output_directory}/${BASE}_variants.Recal.g.vcf
#
# done



# for file in ${output_directory}/${BASE}*_recalibratedNewRef.bam
#
# do
#
#   FBASE=$(basename $file _recalibratedNewRef.bam)
#   BASE=${FBASE%_recalibratedNewRef.bam}
# 	OUT="${BASE}_HaplotypeCaller.sh"
# 	echo "#!/bin/bash" > ${OUT}
# 	echo "#PBS -N ${BASE}_HaplotypeCaller" >> ${OUT}
# 	echo "#PBS -l walltime=72:00:00" >> ${OUT}
# 	echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# 	echo "#PBS -q highmem_q" >> ${OUT}
# 	echo "#PBS -l mem=200gb" >> ${OUT}
# 	echo "" >> ${OUT}
# 	echo "cd ${output_directory}" >> ${OUT}
#   echo "module load ${GATK_module}" >> ${OUT}
# 	echo "" >> ${OUT}
#   echo "time gatk HaplotypeCaller \
#     -R ${ref_genome} \
#     -ERC GVCF \
#     -I ${output_directory}/${BASE}_recalibratedNewRef.bam \
#     -ploidy 2 \
#     -O ${output_directory}/${BASE}_variants.Recal.g.vcf" >> ${OUT}
# 	qsub ${OUT}
#
# done

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
-O ${output_directory}/D0_FullCohort.g.vcf \
-V ${output_directory}/D0-A__variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-1_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-2_variants.Recal.g.vcf \
-V ${output_directory}/D0-3__variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-4_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-5_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-6_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-7_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-9_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-10_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-11_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-12_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-13_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-14_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-15_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-16_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-17_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-18_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-19_variants.Recal.g.vcf \
-V ${output_directory}/D0-20__variants.Recal.g.vcf \
-V ${output_directory}/D0-21__variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-22_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-23_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-24_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-25_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-26_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-27_variants.Recal.g.vcf \
-V ${output_directory}/D0-28__variants.Recal.g.vcf \
-V ${output_directory}/D0-29__variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-30_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-31_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-32_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-33_variants.Recal.g.vcf \
-V ${output_directory}/D0-34__variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-35_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-36_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-37_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-38_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-39_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-40_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-42_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-43_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-44_variants.Recal.g.vcf \
-V ${output_directory}/D0-45__variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-46_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-47_variants.Recal.g.vcf \
-V ${output_directory}/HM-D0-48_variants.Recal.g.vcf
#
#              ###################################################################################################
#              ## Genotype gVCFs (jointly)
#              ###################################################################################################
#              ###################################################################################################
#
# #
time gatk GenotypeGVCFs \
-R ${ref_genome} \
-ploidy 2 \
--variant ${output_directory}/D0_FullCohort.g.vcf \
-O ${output_directory}/D0_FullCohort.vcf

# ###################################################################################################
# ### Find coverage and put into 10k chunks
# ###################################################################################################

module load ${deeptools_module}


for file in ${raw_data}/${BASE}*_piped.bam

do

FBASE=$(basename $file _piped.bam)
BASE=${FBASE%_piped.bam}
OUT="${BASE}_bamCoverage.sh"
echo "#!/bin/bash" > ${OUT}
echo "#PBS -N ${BASE}_bamCoverage" >> ${OUT}
echo "#PBS -l walltime=12:00:00" >> ${OUT}
echo "#PBS -l nodes=1:ppn=1:AMD" >> ${OUT}
echo "#PBS -q batch" >> ${OUT}
echo "#PBS -l mem=20gb" >> ${OUT}
echo "" >> ${OUT}
echo "cd ${raw_data}" >> ${OUT}
echo "module load ${deeptools_module}" >> ${OUT}
echo "" >> ${OUT}
echo "bamCoverage -b ${raw_data}/${BASE}_piped.bam -o ${output_directory}/${BASE}.bedgraph -of bedgraph -bs 10000" >> ${OUT}
qsub ${OUT}

done


for file in ${raw_data}/${BASE}*_piped.bam

do

FBASE=$(basename $file _piped.bam)
BASE=${FBASE%_piped.bam}
samtools sort ${raw_data}/${BASE}_piped.bam \
-o ${raw_data}/${BASE}.sorted.bam

samtools depth \
${raw_data}/${BASE}.sorted.bam \
|  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${raw_data}/${BASE}.txt

done
#
#combine depths with filenames into the same file
find . -type f -name "*.txt" -exec awk '{s=$0};END{if(s)print FILENAME,s}' {} \; > D0_depth.txt



#### bedtools genomecov

module load ${bedtools_module}
# report gives per-base depth across entire genome

bedtools genomecov -d -ibam ${output_directory}/D0/D0-A__recalibratedNewRef.bam > ${output_directory}/D0/D0-A_depth.txt






#                  # ################
# # ###################################################################################################
# # ### Filter variants
# # Can easily run these interactively
# # ###################################################################################################
#

# module load MultiQC/1.5-foss-2016b-Python-2.7.14
# multiqc ${output_directory}/D0

# module load ${samtools_module}
# samtools view -c -F 260 D0-A__pipedNewRef.bam

#### Remove sites with mappability < 0.9
low_mappability="/scratch/jc33471/pilon/337/mappability/337_lowmappability.bed"
module load ${bedtools_module}

# bedtools sort -i ${low_mappability} > ${output_directory}/337_lowmappability_sorted.bed
# bedtools intersect -v -a ${output_directory}/D0/D0_FullCohort_DpGr10_MQGr50_AncCalls_NoHets_Fil.vcf -b ${low_mappability} -header > ${output_directory}/D0/D0_reducedGEM.vcf

# Get only those lines where there is actually a genotype call in the ancestor
# # filter out sites with low read depth
# gatk VariantFiltration \
#    -R ${ref_genome} \
#    -V ${output_directory}/D0/D0_FullCohort.vcf \
#    -O ${output_directory}/D0/D0_FullCohort_DpGr10.vcf \
# 	 --set-filtered-genotype-to-no-call TRUE \
#    -G-filter "DP < 10"  -G-filter-name "depthGr10"


# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/D0/D0_FullCohort_DpGr10.vcf \
# -O ${output_directory}/D0/D0_FullCohort_DpGr10_Fil.vcf --max-filtered-genotypes 0 --exclude-filtered TRUE
#
#
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/D0/D0_FullCohort_DpGr10_Fil.vcf \
# -O ${output_directory}/D0/D0_FullCohort_DpGr10_Fil_AncCalls.vcf \
# -select 'vc.getGenotype("D0-A_").isCalled()'
# #
# # #
# # remove all lines in the ancestor that have a heterozygous genotype
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/D0/D0_FullCohort_DpGr10_Fil_AncCalls.vcf \
# -O ${output_directory}/D0/D0_FullCohort_DpGr10_Fil_AncCalls_NoHets.vcf \
# -select '!vc.getGenotype("D0-A_").isHet()'
#
# #merge filtered with GEM
# bedtools intersect -v -a ${output_directory}/D0/D0_FullCohort_DpGr10_Fil_AncCalls_NoHets.vcf -b ${low_mappability} -header > ${output_directory}/D0/D0_FullCohort_DpGr10_Fil_AncCalls_NoHets_GEM.vcf
#

## Make sure you get the same dataset when you filter the sites from the reduced Gem dataset with my filters and when you combine my filtered datsset with the Gem one
# # First remove the low mappability sites from the D0_FullCohort using bedtools intersect
# bedtools intersect -v -a ${output_directory}/D0/D0_FullCohort.vcf -b ${low_mappability} -header > ${output_directory}/D0/D0_reducedGEM.vcf
#
#
# gatk VariantFiltration \
#    -R ${ref_genome} \
#    -V ${output_directory}/D0/D0_reducedGEM.vcf \
#    -O ${output_directory}/D0/D0_reducedGEM_DpGr10.vcf \
#    --set-filtered-genotype-to-no-call TRUE \
#    -G-filter "DP < 10"  -G-filter-name "depthGr10"
#
# gatk SelectVariants \
# 	    -R ${ref_genome} \
# 	    -V ${output_directory}/D0/D0_reducedGEM_DpGr10.vcf \
# 	    -O ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil.vcf --max-filtered-genotypes 0 \
# 	    --exclude-filtered TRUE
#
# #prints rows of file2 that are not in file1
# awk 'NR==FNR{a[$1,$2]; next} !(($1,$2) in a)' ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil.vcf ${output_directory}/D0/D0_reducedGEM.vcf > ${output_directory}/D0/filteredOut.vcf
#









########################################################################
#### Remove low and high read depth first
gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/D0_FullCohort.vcf \
-O ${output_directory}/D0_noLow.vcf \
-select 'vc.getGenotype("D0-A_").getDP() > 82'

gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/D0_noLow.vcf \
-O ${output_directory}/D0_noLow_noHigh.vcf \
-select 'vc.getGenotype("D0-A_").getDP() < 206'

low_mappability="/scratch/jc33471/pilon/337/mappability/337_lowmappability.bed"
module load ${bedtools_module}

# bedtools sort -i ${low_mappability} > ${output_directory}/337_lowmappability_sorted.bed
bedtools intersect -v -a ${output_directory}/D0_noLow_noHigh.vcf -b ${low_mappability} -header > ${output_directory}/D0_noLow_noHigh_redGem.vcf

awk 'NR==FNR{a[$1,$2]; next} !(($1,$2) in a)' ${output_directory}/D0/D0_noLow_noHigh_redGem.vcf ${output_directory}/D0/D0_noLow_noHigh.vcf > ${output_directory}/D0/GEMremoved.txt


gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/D0_noLow_noHigh_redGem.vcf \
-O ${output_directory}/D0_noLow_noHigh_redGem_AncCalls.vcf \
-select 'vc.getGenotype("D0-A_").isCalled()'


gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/D0_noLow_noHigh_redGem_AncCalls.vcf \
-O ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets.vcf \
-select '!vc.getGenotype("D0-A_").isHet()'



### Ancestor hets only
gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/D0_noLow_noHigh_redGem_AncCalls.vcf \
-O ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_Hets.vcf \
-select 'vc.getGenotype("D0-A_").isHet()'




### select snps and indels and make into tables
gatk SelectVariants \
   -R ${ref_genome} \
   -V ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets.vcf \
   -O ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets_SNPs.vcf \
   --max-nocall-number 0 \
   --exclude-non-variants TRUE \
	 --restrict-alleles-to BIALLELIC \
   -select-type SNP
#
gatk SelectVariants \
   -R ${ref_genome} \
   -V ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets.vcf \
   -O ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets_Indels.vcf \
   --max-nocall-number 0 \
	 --exclude-non-variants TRUE \
	 --restrict-alleles-to BIALLELIC \
   -select-type INDEL

gatk SelectVariants \
	 -R ${ref_genome} \
   -V ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets.vcf \
   -O ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHetsVars.vcf \
   --max-nocall-number 0 \
	 --exclude-non-variants TRUE \
	 --restrict-alleles-to BIALLELIC


gatk VariantsToTable \
	 -V ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHetsVars.vcf \
	 -F CHROM -F POS -F REF -F ALT -F QUAL \
	 -GF AD -GF DP -GF GQ -GF GT \
	 -O ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets_vars.txt

gatk SelectVariants \
	 	    -R ${ref_genome} \
	 	    -V ${output_directory}/D0_FullCohort.vcf \
	 	    -O ${output_directory}/D0_FullCohortCalls.vcf \
	 	    --max-nocall-number 0

	 gatk VariantsToTable \
	 	 -V ${output_directory}/D0_FullCohortCalls.vcf \
	 	 -F CHROM -F POS -F REF -F ALT -F QUAL \
	 	 -GF AD -GF DP -GF GQ -GF GT \
	 	 -O ${output_directory}/D0_FullCohortCalls.txt

gatk VariantsToTable \
	-V ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets_SNPs.vcf \
	-F CHROM -F POS -F REF -F ALT -F QUAL \
	-GF AD -GF DP -GF GQ -GF GT \
	-O ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets_SNPs.txt

gatk VariantsToTable \
-V ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets_Indels.vcf \
-F CHROM -F POS -F REF -F ALT -F QUAL \
-GF AD -GF DP -GF GQ -GF GT \
-O ${output_directory}/D0_noLow_noHigh_redGem_AncCalls_NoHets_Indels.txt










awk 'NR==FNR{a[$1,$2]; next} !(($1,$2) in a)' ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil_AncDpTest2_NoHets.vcf ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil_AncCalls_NoHets.vcf > ${output_directory}/D0/filtered_out.vcf

#
# #
# remove all lines in the ancestor that have a heterozygous genotype
gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil_AncCalls.vcf \
-O ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil_AncCalls_NoHets.vcf \
-select '!vc.getGenotype("D0-A_").isHet()'

gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil_AncCalls_NoHets.vcf \
-O ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil_AncCalls_NoHets_AncLowDepthFilt.vcf \
-select 'vc.getGenotype("D0-A_").getDP() > 114'

gatk SelectVariants \
-R ${ref_genome} \
-V ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil_AncCalls_NoHets_AncLowDepthFilt.vcf \
-O ${output_directory}/D0/D0_reducedGEM_DpGr10_Fil_AncCalls_NoHets_AncLowDepthFilt_HighDepthFilt.vcf \
-select 'vc.getGenotype("D0-A_").getDP() < 254'


bedtools intersect -v -a ${output_directory}/D0/D0_FullCohort_DpGr10_MQGr50_AncCalls_NoHets_Fil.vcf -b ${low_mappability} -header > ${output_directory}/D0/D0_reducedGEM.vcf



## Find any differences between D0_reducedGEM_DpGr10_Fil_AncCalls_NoHets.vcf and D0_FullCohort_DpGr10_MQGr50_StrBiasFil_AncCalls_NoHets_GEM.vcf
# awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' D0_reducedGEM_DpGr10_Fil_AncCalls_NoHets.vcf D0_FullCohort_DpGr10_Fil_AncCalls_NoHets_GEM.vcf > common.txt



### Find how many sites are variable in the ancestor
## Heterozygous sites
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/D0/reducedTest.vcf \
# -O ${output_directory}/D0/reducedTestAncHets.vcf \
# -select 'vc.getGenotype("D0-A_").isHet()'

## Homozygous variant sites
# gatk SelectVariants \
# -R ${ref_genome} \
# -V ${output_directory}/D0/reducedTest.vcf \
# -O ${output_directory}/D0/reducedTestAncHomVars.vcf \
# -select 'vc.getGenotype("D0-A_").isHomVar()'

#

#
# #
#
# # cd ${output_directory}
# #
gatk SelectVariants \
   -R ${ref_genome} \
   -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
   -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_SNPs.vcf \
   --max-nocall-fraction 0 \
   --exclude-non-variants TRUE \
   -select-type SNP
#
gatk SelectVariants \
   -R ${ref_genome} \
   -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
   -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_Indels.vcf \
      --max-nocall-fraction 0.001 \
      -select-type INDEL

# #
# # #gives a final dataset with only called sites in the Ancestor, no heterozygous sites in the ancestor,
# # # depth > 10, mapping quality > 50, and strand bias (SOR) > 0.01 (not significant)
# #
# # #Variants to table
gatk VariantsToTable \
     -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls.vcf \
     -F CHROM -F POS -F REF -F ALT -F QUAL \
     -GF AD -GF DP -GF GQ -GF GT \
     -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_vars.txt

gatk VariantsToTable \
          -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_SNPs.vcf \
          -F CHROM -F POS -F REF -F ALT -F QUAL \
          -GF AD -GF DP -GF GQ -GF GT \
          -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_SNPs.txt

gatk VariantsToTable \
               -V ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_Indels.vcf \
               -F CHROM -F POS -F REF -F ALT -F QUAL \
               -GF AD -GF DP -GF GQ -GF GT \
               -O ${output_directory}/D0_FullCohort_AnCalls_NoHets_DpGr10_MQGr50_StrBiasFil_Calls_Indels.txt



awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' reducedTest.vcf D0_FullCohort_DpGr10_MQGr50_StrBiasFil_AncCalls_NoHets.vcf > common.txt

#prints rows of file2 that are not in file1
awk 'NR==FNR{a[$1,$2]; next} !(($1,$2) in a)' file1 file2

awk 'NR==FNR{a[$1,$2]; next} !(($1,$2) in a)' reducedTest.vcf D0_FullCohort_DpGr10_MQGr50_AncCalls_NoHets_Fil.vcf > uniqueMine.txt

awk 'NR==FNR{a[$1,$2]; next} !(($1,$2) in a)' D0_FullCohort_DpGr10_MQGr50_AncCalls_NoHets_Fil.vcf reducedTest.vcf > uniqueToGem.txt




### Find common heterozygous sites between the three diploid ancestors
awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' D0_noLow_noHigh_redGem_AncCalls_Hets.vcf D1_noLow_noHigh_redGem_AncCalls_Hets.vcf > D0_D1_common.txt

awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' D20_noLow_noHigh_redGem_AncCalls_Hets.vcf D0_D1_common.txt > D0_D1_D20_common.txt




awk 'NR==FNR{a[$1]; next} (($1) in a)' D20-44_R1_val_telocate_raw.TY1.bed D20-A-R_R1_val_telocate_raw.TY1.bed > common.txt



##################################################################################################
### This stuff is in Excel
# # ###################################################################################################
#
# 3 columns after last genotype: # samples (count # samples you have), #same (meaning number of samples that have the same genotype as the ancestor),
# and # variant (# samples different genotype than the ancestor).
# Codes are:
#same: =COUNTIFS(E2:CL2,CN2)
#variant: =CO2-CP2 ( just subtract the # same from the # samples)
