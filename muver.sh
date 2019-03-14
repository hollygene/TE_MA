#PBS -S /bin/bash
#PBS -N muver
#PBS -q highmem_q
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l mem=200gb
#PBS -l walltime=72:00:00


cd $PBS_O_WORKDIR
muver_module="muver/0.1.0-foss-2016b-Python-2.7.14"
trimgalore_module="Trim_Galore/0.4.5-foss-2016b"
data_dir="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/samples"
ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/S288C_reference_sequence_R64-2-1_20150113.fsa"
ref_genome_dir="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/"
fastq_list="/home/hcm14449/Github/TE_MA/FASTQ_LIST.txt"
control_sample_name="Ancestor"
experiment_directory="/scratch/hcm14449/TE_MA_Paradoxus/Practice/output.3.14.19"
mkdir $experiment_directory
trimmed_data="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/samples/trimmed"
# mkdir $trimmed_data

# module load ${trimgalore_module}
#
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
# cd $ref_genome_dir
# wget https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz

module load ${muver_module}
# java -jar $EBROOTGATK/GenomeAnalysisTK.jar

#index reference genome
# echo ${ref_genome}
muver index-reference ${ref_genome}
#
# create repeat file for reference genome
muver create-repeat-file ${ref_genome}

# run the pipeline
muver run-pipeline ${ref_genome} ${fastq_list} ${control_sample_name} ${experiment_directory}

# muver create-repeat-file /scratch/hcm14449/TE_MA_Paradoxus/Practice/files/samples/Sample2_R1.fastq
#
# muver fit-repeat-indel-rates /scratch/hcm14449/TE_MA_Paradoxus/Practice/output.3/bams/Sample2.bam

# muver run-pipeline /scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome.fa.fai \
# /home/hcm14449/Github/TE_MA/FASTQ_LIST.txt Ancestor /scratch/hcm14449/TE_MA_Paradoxus/Practice/output
