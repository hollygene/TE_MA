#PBS -S /bin/bash
#PBS -q batch
#PBS -N testScriptJuly19
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=100gb
#PBS -M hcm14449@uga.edu
#PBS -m abe

cd /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data
wget http://web.qbcg.uga.edu:8888/fop/82kNLuio/Holly.gz.gz
