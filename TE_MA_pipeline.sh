#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N testScriptJuly19
#PBS -l nodes=1:ppn=12
#PBS -l walltime=480:00:00
#PBS -l mem=400gb
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
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out"
# mkdir $output_directory
#location of TRIMMED data to be used in the analysis
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq"


########################################################################################################################
# works
# cd ${data_dir}
#
mkdir ${output_directory}
#
# gunzip /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/*.gz
#
# module load ${fastqc_module}
#
# for file in ${data_dir}/*.fastq
#
# do
#
# FBASE=$(basename $file .fastq)
# BASE=${FBASE%.fastq}
#
# time fastqc -o ${output_directory} ${data_dir}/${BASE}.fastq
#
# done


#######################################################################################
# works
# cd ${data_dir}
#
# mkdir ${trimmed_data}
#
# module load ${trimgalore_module}
#
# #need to remove TruSeq adapters Index 6
#
#
#
# # trim all fastq files
# for file in $data_dir/*.fastq
#
# do
#
# FBASE=$(basename $file .fastq)
# BASE=${FBASE%.fastq}
#
# trim_galore --phred33 -q 20 -o $trimmed_data ${BASE}.fastq
#
# done
#
# module unload ${trimgalore_module}

 #######################################################################################
# works
# thinks the arabidopsis files are still there

module load ${bwa_module}

 #index the ref genome
# bwa index ${ref_genome}

for file in ${raw_data}/*_R1_001.fastq

 do

 FBASE=$(basename $file _R1_001.fastq)
 BASE=${FBASE%_R1_001.fastq}

bwa mem -M -t 12 ${ref_genome} ${raw_data}/${BASE}_R1_001.fastq ${raw_data}/${BASE}_R2_001.fastq > ${output_directory}/${BASE}_aln.sam

 done

#########################################################################################

#samtools
module load ${samtools_module}

#index reference genome

samtools faidx ${ref_genome}

#convert sam files to bam files
for file in ${output_directory}/*_aln.sam

do

FBASE=$(basename $file _aln.sam)
BASE=${FBASE%_aln.sam}

samtools view -bt ${ref_genome_dir}/*.fai \
${output_directory}/${BASE}_aln.sam \
  > ${output_directory}/${BASE}.bam

done

### sort the bam files
for file in ${output_directory}/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}

samtools sort -o ${output_directory}/${BASE}.sorted.bam \
   ${output_directory}/${BASE}.bam

done


#################################################################################################################
#bam to BigWig > for quality control purposes
module load ${bedtools_module}
module load ${python_module}
module load ${samtools_module}
export PATH=${PATH}:${script_location}

## Loop
for file in ${output_directory}/*.sorted.bam

do

FBASE=$(basename $file .sorted.bam)
BASE=${FBASE%.sorted.bam}

python3 ${bamToBigWig} -sort ${ref_genome_dir}/*.fai ${output_directory}/${BASE}.sorted.bam

done

###################################################################################################
## Picard to mark duplicates

module load ${picard_module}

for file in ${output_directory}/*.bam

do

FBASE=$(basename $file .bam)
BASE=${FBASE%.bam}

time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
I=${output_directory}/${BASE}.bam \
O=${output_directory}/${BASE}_markedDuplicates.bam \
M=${output_directory}/${BASE}_markedDupsMetrics.txt

done

###################################################################################################

##Genotype 4 random samples from each experiment
# Using GATK HaplotypeCaller in GVCF mode
# apply appropriate ploidy for each sample

# module load ${GATK_module}
#
# time gatk HaplotypeCaller \
#      -R ${ref_genome} \
#      -I ${H0_bams}/H0-A.bam \
#      -ploidy 1 \
#      -O ${H0_bams}/H0-A_variants.vcf
