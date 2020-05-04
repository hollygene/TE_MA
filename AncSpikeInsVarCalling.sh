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
output_directory="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/Anc_SpikeIns"
mkdir $output_directory
#location of data to be used in the analysis
raw_data="/project/dwhlab/Holly/TE_MA_Paradoxus/Paradoxus_MA/Anc_SpikeIns/Holly_gDNA"

# cd ${raw_data}
cd ${output_directory}

module load ${picard_module}
module load ${bwa_module}
module load ${samtools_module}
module load ${GATK_module}



for file in ${raw_data}/*.fastq.gz

do

FBASE=$(basename $file .fastq.gz)
BASE=${FBASE%_R1_001.fastq}

java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
    FASTQ=${raw_data}/${BASE}R1.fq.gz \
    FASTQ2=${raw_data}/${BASE}R2.fq.gz \
    OUTPUT=${output_directory}/${BASE}_fastqtosam.bam \
    READ_GROUP_NAME=${BASE} \
    SAMPLE_NAME=${BASE} \
    PLATFORM=illumina \
    SEQUENCING_CENTER=GGBC

done
