#PBS -S /bin/bash
#PBS -N j_s_samtools
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00

cd $PBS_O_WORKDIR
