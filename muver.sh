#PBS -S /bin/bash
#PBS -N muver
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=50gb
#PBS -l walltime=480:00:00
module load muver/0.1.0-foss-2016b-Python-2.7.14

cd $PBS_O_WORKDIR
# mkdir "/scratch/hcm14449/TE_MA_Paradoxus/Practice/output"
data_dir="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/samples"
ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome"
fastq_list="/home/hcm14449/Github/TE_MA/FASTQ_LIST.txt"
control_sample_name="Ancestor"
experiment_directory="/scratch/hcm14449/TE_MA_Paradoxus/Practice/output"

# java -jar $EBROOTGATK/GenomeAnalysisTK.jar

#index reference genome
# echo ${ref_genome}
# muver index-reference ${ref_genome}
#
# create repeat file for reference genome
# muver create-repeat-file ${ref_genome}

# run the pipeline
muver run-pipeline ${ref_genome} ${fastq_list} ${control_sample_name} ${experiment_directory}

# muver run-pipeline /scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome.fa.fai \
# /home/hcm14449/Github/TE_MA/FASTQ_LIST.txt Ancestor /scratch/hcm14449/TE_MA_Paradoxus/Practice/output
