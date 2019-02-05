#PBS -S /bin/bash
#PBS -N muver
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=50gb
#PBS -l walltime=480:00:00

cd $PBS_O_WORKDIR

data_dir="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/samples"

ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome.fa"

module load muver/0.1.0-foss-2016b-Python-2.7.14

#index reference genome
muver index-reference ${ref_genome}

# create repeat file for reference genome
# muver create-repeat-file /scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome.fa > /scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome_Repeats.txt
