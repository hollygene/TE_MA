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
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq"



# #########################################################################################
# #samtools: converts sam files to bam files and sorts them
# #########################################################################################
#
module load ${samtools_module}

# #index reference genome
#
# samtools faidx ${ref_genome}
#
# # #convert sam files to bam files
# for file in ${output_directory}/H0/*_aln.sam
#
# do
#
# FBASE=$(basename $file _aln.sam)
# BASE=${FBASE%_aln.sam}
#
# samtools view -bt ${ref_genome_dir}/*.fai \
# ${output_directory}/H0/${BASE}_aln.sam \
#   > ${output_directory}/H0/${BASE}.bam
#
# done



# ############################
# ### sort the bam files
# ############################

# for file in ${output_directory}/H0/*.bam
#
# do
#
# FBASE=$(basename $file .bam)
# BASE=${FBASE%.bam}
#
# samtools sort -@ 12 -o ${output_directory}/H0/${BASE}.sorted.bam \
#    ${output_directory}/H0/${BASE}.bam
#
# done

# ############################
# ### index the bam files
# ############################

# for file in ${output_directory}/H0/*.sorted.bam
#
# do
#
# FBASE=$(basename $file .sorted.bam)
# BASE=${FBASE%.sorted.bam}
#
# samtools index -@ 12 ${output_directory}/H0/${BASE}.sorted.bam
#
# done
# ###################################################################################################
# ## Picard to mark duplicates
# ###################################################################################################
# #
module load ${picard_module}

for file in ${output_directory}/*.sorted.bam

do

FBASE=$(basename $file .sorted.bam)
BASE=${FBASE%.sorted.bam}

time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar ValidateSamFile \
      I=${output_directory}/${BASE}.sorted.bam \
      MODE=SUMMARY

done

# for file in ${output_directory}/H0/*.sorted.bam
#
# do
#
# FBASE=$(basename $file .sorted.bam)
# BASE=${FBASE%.sorted.bam}
#
# time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=TRUE \
# I=${output_directory}/H0/${BASE}.sorted.bam \
# O=${output_directory}/H0/${BASE}_removedDuplicates.bam \
# M=${output_directory}/H0/${BASE}_removedDupsMetrics.txt
#
# done


###################################################################################################
# Using GATK HaplotypeCaller in GVCF mode
# apply appropriate ploidy for each sample
# will need to do this separtely for haploid and diploid samples
###################################################################################################

# module load ${GATK_module}
#
# ### H0 samples
# for file in ${output_directory}/H0/*_removedDuplicates.bam
#
# do
#
# FBASE=$(basename $file _removedDuplicates.bam)
# BASE=${FBASE%_removedDuplicates.bam}
#
#
# time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -I ${output_directory}/H0/${BASE}_removedDuplicates.bam \
#      -ploidy 1 \
#      -O ${output_directory}/H0/${BASE}_variants.g.vcf
#
# done


###################################################################################################
### Aggregate the GVCF files using GenomicsDBImport
###################################################################################################
# mkdir /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB/tmp
#
# gatk  --java-options "-Xmx4g -Xms4g" \
#        GenomicsDBImport \
#        --genomicsdb-workspace-path /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB \
#        --batch-size 50 \
#        --sample-name-map /home/hcm14449/Github/TE_MA/H0_sample_map.txt \
#        --TMP_DIR:/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/GenDB/tmp \
#        --reader-threads 12
