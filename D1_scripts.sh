#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N D1_recal_gvcfs
#PBS -l nodes=2:ppn=1:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=600gb
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
deeptools_module="deepTools/3.2.1-foss-2018a-Python-3.6.4"
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
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1"
# mkdir $output_directory
#location of TRIMMED data to be used in the analysis
# raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq"
# trimmed_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/D1"
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/D1"
# do_again="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/D1/do_again"
genomicsdb_workspace_path="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/GenDB"
sample_name_map="/home/hcm14449/Github/TE_MA/D1_sample_map.txt"
tmp_DIR="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/GenDB/tmp"

cd ${output_directory}
rm *

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
#     LIBRARY_NAME=D1 \
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

#######################################################################################
# #
# #
# module load ${picard_module}
#
#
# for file in ${raw_data}/*_markilluminaadapters.bam
#
# do
#
# FBASE=$(basename $file _markilluminaadapters.bam)
# BASE=${FBASE%_markilluminaadapters.bam}
#
# time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar ValidateSamFile \
#       I=${raw_data}/${BASE}_markilluminaadapters.bam \
#       MODE=VERBOSE
#
# done

#######################################################################################
# convert BAM to FASTQ and discount adapter sequences using SamToFastq
#######################################################################################

# for file in ${output_directory}/*_markilluminaadapters.bam
#
# do
#
# FBASE=$(basename $file _markilluminaadapters.bam)
# BASE=${FBASE%_markilluminaadapters.bam}
#
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/picard/2.16.0-Java-1.8.0_144 SamToFastq \
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

### Piped command: SamToFastq, then bwa mem, then MergeBamAlignment

# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
# I=${raw_data}/HM-D1-3_markilluminaadapters.bam \
# FASTQ=/dev/stdout \
# CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
# TMP_DIR=${raw_data}/TMP | \
# bwa mem -M -t 7 -p ${ref_genome} /dev/stdin| \
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
# ALIGNED_BAM=/dev/stdin \
# UNMAPPED_BAM=${raw_data}/HM-D1-3_fastqtosam.bam \
# OUTPUT=${do_again}/HM-D1-3_piped.bam \
# R=${ref_genome} CREATE_INDEX=true ADD_MATE_CIGAR=true \
# CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
# INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
# PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
# TMP_DIR=${raw_data}/TMP

#######################################################################################
# works: aligns samples to reference genome. Output is a .sam file
#######################################################################################

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
#
# # #########################################################################################
# # #samtools: converts sam files to bam files and sorts them
# # #########################################################################################
# #
# # #convert sam files to bam files
# module load ${samtools_module}
#
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
# # ############################
# # ### sort the bam files
# # ############################
# #
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
# samtools index -@ 12 -o ${output_directory}/${BASE}.sorted.bam
#
# done
# # ###################################################################################################
# # ## Picard to mark duplicates
# # ###################################################################################################

#
# # ###################################################################################################
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
# ###################################################################################################
# #
module load ${GATK_module}

## D1 samples
for file in ${raw_data}/${BASE}*_piped.bam

do

FBASE=$(basename $file _piped.bam)
BASE=${FBASE%_piped.bam}

time gatk HaplotypeCaller \
     -R ${ref_genome} \
     -ERC GVCF \
     -I ${raw_data}/${BASE}_piped.bam \
     -ploidy 2 \
     -O ${output_directory}/${BASE}_variants.g.vcf

done


# module load GATK/4.0.3.0-Java-1.8.0_144
#
# time gatk HaplotypeCaller \
#      -R /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa \
#      -ERC GVCF \
#      -I /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/D0/HM-D0-10_piped.bam \
#      -ploidy 2 \
#      -O /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/HM-D0-10_variants.g.vcf
#
# ###################################################################################################
### Combine gVCFs before joint genotyping
# ###################################################################################################


time gatk CombineGVCFs \
 -O ${output_directory}/D1_cohort.g.vcf \
 -R ${ref_genome} \
 --variant ${output_directory}/HM-D1-A_variants.g.vcf \
 --variant ${output_directory}/HM-D1-10_variants.g.vcf \
 --variant ${output_directory}/HM-D1-11_variants.g.vcf \
 --variant ${output_directory}/HM-D1-12_variants.g.vcf \
 --variant ${output_directory}/HM-D1-13_variants.g.vcf \
 --variant ${output_directory}/HM-D1-14_variants.g.vcf \
 --variant ${output_directory}/HM-D1-15_variants.g.vcf \
 --variant ${output_directory}/HM-D1-16_variants.g.vcf

#
# # ###################################################################################################
# # ### Jointly genotype 8 random samples to identify consensus sequences
# # ###################################################################################################
#
time gatk GenotypeGVCFs \
        -R ${ref_genome} \
        --variant ${output_directory}/D1_cohort.g.vcf \
        -O ${output_directory}/D1_variants_8Samples.vcf

# ###################################################################################################
# ## Recalibrate base quality scores in all samples to mask any likely consensus variants
# ###################################################################################################

for file in ${raw_data}/${BASE}*_piped.bam

do

FBASE=$(basename $file _piped.bam)
BASE=${FBASE%_piped.bam}


time gatk BaseRecalibrator \
   -I ${raw_data}/${BASE}_piped.bam \
   --known-sites ${output_directory}/D1_variants_8Samples.vcf \
   -O ${output_directory}/${BASE}_recal_data.table \
   -R ${ref_genome}
done

# ###################################################################################################
# ## Apply BQSR to bam files
# ###################################################################################################

for file in ${raw_data}/${BASE}*_piped.bam

do

FBASE=$(basename $file _piped.bam)
BASE=${FBASE%_piped.bam}


gatk ApplyBQSR \
   -R ${ref_genome} \
   -I ${raw_data}/${BASE}_piped.bam \
   -bqsr ${output_directory}/${BASE}_recal_data.table \
   -O ${output_directory}/${BASE}_recalibrated.bam

done

# ###################################################################################################
### Run HaplotypeCaller again on recalibrated samples
# ###################################################################################################
# ###################################################################################################
# #
# module load ${GATK_module}

## D1 samples
for file in ${output_directory}/${BASE}*_recalibrated.bam

do

FBASE=$(basename $file _recalibrated.bam)
BASE=${FBASE%_recalibrated.bam}

time gatk HaplotypeCaller \
     -R ${ref_genome} \
     -ERC GVCF \
     -I ${output_directory}/${BASE}_recalibrated.bam \
     -ploidy 2 \
     -O ${output_directory}/${BASE}_variants.Recal.g.vcf

done
# ###################################################################################################
### Run GenotypeGVCFs on recalibrated samples
# ###################################################################################################
# ###################################################################################################
# #
for file in ${output_directory}/${BASE}*_variants.g.vcf

do

FBASE=$(basename $file _variants.g.vcf)
  BASE=${FBASE%_variants.g.vcf}

time gatk GenotypeGVCFs \
     -R ${ref_genome} \
     --variant ${output_directory}/${BASE}_variants.g.vcf \
     -O ${output_directory}/${BASE}.vcf

done

# filter based on genotype in the Ancestor only
# time gatk GenotypeGVCFs \
#      -R ${ref_genome} \
#      --variant ${output_directory}/HM-D1-A_variants.g.vcf \
#      -O ${output_directory}/HM-D1-A.vcf
#
#
# gatk VariantFiltration \
# -V ${output_directory}/HM-D1-A.vcf \
# -O ${output_directory}/HM-D1-A_hetsFil.vcf \
# --genotype-filter-expression "isHet == 1" \
# --genotype-filter-name "isHetFilter"
#
#
# gatk SelectVariants \
# -V ${output_directory}/HM-D1-A_hetsFil.vcf \
# --set-filtered-gt-to-nocall \
# -O ${output_directory}/HM-D1-A_hetsFil_SV.vcf
#
#
# time gatk SelectVariants \
#        -R ${ref_genome} \
#        -V ${output_directory}/HM-D1-A_hetsFil_SV.vcf \
#        -O ${output_directory}/HM-D1-A_hetsFil_SV_noNoCalls.vcf \
#        --max-nocall-number 0
# # ###################################################################################################
# ### Combine VCFs together
# # ###################################################################################################
# module load GATK/3.8-1-Java-1.8.0_144
#
# java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#    -T CombineVariants \
#    -R ${ref_genome} \
#    -V ${output_directory}/HM-D1-A_hetsFil_SV_noNoCalls.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-10.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-20.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-31.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-42.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-21.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-32.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-44.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-11.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-22.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-33.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-45.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-12.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-24.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-35.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-46.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-13.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-25.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-36.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-4.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-14.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-26.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-37.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-5.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-15.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-27.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-38.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-6.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-16.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-28.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-39.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-7.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-17.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-29.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-3.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-8.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-18.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-2.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-9.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-19.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-30.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-40.vcf \
#       -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-41.vcf \
#       -o ${output_directory}/fullCohort.vcf \
#       -genotypeMergeOptions UNIQUIFY
#
# module unload GATK/3.8-1-Java-1.8.0_144



# module load ${GATK_module}

# gatk VariantsToTable \
#             -V ${output_directory}/fullCohort.vcf \
#             -F CHROM -F POS -F REF -F ALT -F QUAL -F DP -F GT \
#             -GF AD -GF GT -GF DP \
#             -O ${output_directory}/fullCohort.txt
# ###################################################################################################
### Run HaplotypeCaller on recalibrated samples
# ###################################################################################################
# ###################################################################################################
#
for file in ${output_directory}/do_again/${BASE}*_recalibrated.bam

do

FBASE=$(basename $file _recalibrated.bam)
BASE=${FBASE%_recalibrated.bam}

time gatk HaplotypeCaller \
     -R ${ref_genome} \
     -ERC GVCF \
     -I ${output_directory}/do_again/${BASE}_recalibrated.bam \
     -ploidy 2 \
     -O ${output_directory}/do_again/${BASE}_variants.Recal.g.vcf

done

# time gatk CombineVariants \
#     -R reference.fasta \
#     --variant ${output_directory}/HM-D1-10.vcf \
#     --variant ${output_directory}/HM-D1-11.vcf \
#     -o ${output_directory}/combined10_11.vcf \
#     -genotypeMergeOptions UNIQUIFY

# ###################################################################################################
# ### Aggregate the GVCF files using GenomicsDBImport
# ###################################################################################################
# mkdir ${genomicsdb_workspace_path}
# mkdir ${tmp_DIR}
# #
# time gatk GenomicsDBImport \
#        --genomicsdb-workspace-path ${genomicsdb_workspace_path} \
#        --sample-name-map ${sample_name_map}



time gatk CombineGVCFs \
   -R ${ref_genome} \
   -O /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/cohort.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-A_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-10_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-20_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-31_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-42_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-21_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-32_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-44_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-11_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-22_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-33_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-45_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-12_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-24_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-35_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-46_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-13_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-25_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-36_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-4_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-14_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-26_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-37_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-5_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-15_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-27_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-38_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-6_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-16_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-28_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-39_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-7_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-17_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-29_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-3_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-8_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-18_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-2_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-9_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-19_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-30_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-40_variants.Recal.g.vcf \
   -V /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-41_variants.Recal.g.vcf
#
#        # ###################################################################################################
#        # ### joint genotype vcfs
#        # ###################################################################################################
#
time gatk GenotypeGVCFs \
         -R ${ref_genome} \
         --variant /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/cohort.g.vcf \
         -O ${output_directory}/Full_cohort_D1.vcf





         # ###################################################################################################
         # ### Find coverage and put into 10k chunks
         # ###################################################################################################

         module load ${deeptools_module}


         for file in ${raw_data}/${BASE}*_piped.bam

         do

         FBASE=$(basename $file _piped.bam)
         BASE=${FBASE%_piped.bam}

         bamCoverage -b ${raw_data}/${BASE}_piped.bam -o ${output_directory}/${BASE}.bedgraph -of bedgraph -bs 10000

         done

# time gatk SelectVariants \
#             -R ${ref_genome} \
#             -V ${output_directory}/Full_cohort.vcf \
#             -O ${output_directory}/Full_cohort_noNonVars.vcf \
#             --exclude-non-variants true
#
#           time gatk SelectVariants \
#                         -R ${ref_genome} \
#                         -V ${output_directory}/Full_cohort_noNonVars.vcf \
#                         -O ${output_directory}/Full_cohort_noNonVars_noFils.vcf \
#                         --exclude-filtered true



# time gatk SelectVariants \
#        -R ${ref_genome} \
#        -V ${output_directory}/Full_cohort.vcf \
#        -O ${output_directory}/test.vcf \
#        -sample-expressions '!vc.getGenotype('HM-D1-A').isHet()'

       # time gatk SelectVariants \
       #        -R ${ref_genome} \
       #        -V ${output_directory}/Full_cohort_VF_SV.vcf \
       #        -O ${output_directory}/Full_cohort_VF_SV_noNoCalls.vcf \
       #        --max-nocall-number 45
       #
       #
       #
       #        gatk VariantsToTable \
       #             -V ${output_directory}/Full_cohort_VF_SV_noNoCalls.vcf \
       #             -F CHROM -F POS -F REF -F ALT -F QUAL -F DP \
       #             -GF AD -GF GT -GF DP \
       #             -O ${output_directory}/Full_cohort_VF_SV_noNoCalls.txt


# -G_filter isHet == 1
# -G_filterName
# --filterName One --filterExpression "X < 1" --filterName Two --filterExpression "X > 2"



# gatk VariantFiltration \
# -V ${output_directory}/Full_cohort.vcf \
# -O ${output_directory}/Full_cohort_VF.vcf \
# --genotype-filter-expression "isHet == 1" \
# --genotype-filter-name "isHetFilter"
#
#
# gatk SelectVariants \
# -V ${output_directory}/Full_cohort_VF.vcf \
# --set-filtered-gt-to-nocall \
# -O ${output_directory}/Full_cohort_VF_SV.vcf


### What if I individually genotype all of my GVCFs
# Then in the Ancestor, remove everything that is heterozygous
# Then combine my VCFs based on POS and CHROM
# That should give me a giant file with every site in every line that is NOT variant in the Ancestor
      # in theory...
# env means esclude not variant loci
# ef means exclude filtered loci

            # vc.getGenotype('sample08').isHomVar()
            # vc.getGenotype('sample08').isHomRef()
