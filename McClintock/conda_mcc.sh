
#!/bin/bash
#PBS -q highmem_q
#PBS -N MK_yeast
#PBS -l nodes=1:ppn=12 -l mem=120gb
#PBS -l walltime=30:00:00
#PBS -M hcm14449@uga.edu
#PBS -m abe
#PBS -o /lustre1/hcm14449/mcc/mcc_0823_yeast_bioconda_RM7_all.out
#PBS -e /lustre1/hcm14449/mcc/mcc_0823_yeast_bioconda_RM7_all.out
#PBS -j oe


# load anaconda module
module load Anaconda3/5.0.1


# you should add the part to create new environment with required dependencies and activate the environment
#only need to do create once
conda create -y -n MCCLINTOCK

source activate MCCLINTOCK

#only need to do install once
#repeatmasker might be different for GFF files (diff version)

conda install -y repeatmasker=4.0.7
conda install -y perl-bioperl-run=1.006900
conda install -y bwa=0.7.4
conda install -y bedtools=2.17.0
conda install -y bowtie=1.0.0
conda install -y ucsc-blat=366
conda install -y samtools=0.1.19
conda install -y r=3.5.0
conda install -y fastqc=0.11.2
conda install -y ucsc-twobittofa=366
conda install -y bcftools=1.2
conda install -y exonerate=2.4.0


# export RM
# export PATH=$PATH:/home/sh60271/app/RepeatMasker


# install McC
mcc_dir="/lustre1/hcm14449/mcc/mcc_software_bioconda_test"
#only need to install mcc once
cp -r /home/sh60271/git/mcclintock $mcc_dir
cd $mcc_dir
sh install.sh


# change to your mcc installation dir if you plan to use your installed version
run_dir="/lustre1/hcm14449/TE_MA_Paradoxus/Anc_SpikeIns/DNA_copy/8_23_18"
mkdir -p $run_dir
rm -rf $run_dir/*

#ref genome dir
ref_dir="/lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/cerevisiae"

# use existing data
data_dir="/lustre1/hcm14449/TE_MA_Paradoxus/Anc_SpikeIns/DNA_copy/data"
out_dir=$run_dir/out
mkdir -p $data_dir
mkdir -p $out_dir

# bash $mcc_dir/mcclintock.sh -d -C -r $data_dir/sacCer2.fasta -c \
# $data_dir/sac_cer_TE_seqs.fasta -1 $data_dir/SRR800842_1.fastq -2 \
# $data_dir/SRR800842_2.fastq -p $PBS_NP -b -o $out_dir

#single end
bash $mcc_dir/mcclintock.sh -d -C -r $ref_dir/sacCer2.fasta -c \
$mcc_dir/test/sac_cer_TE_seqs.fasta -1 $data_dir/HM_H0_S16_R1_001.fq \
-p $PBS_NP -b -o $out_dir
