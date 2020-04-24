#PBS -S /bin/bash
#PBS -q batch
#PBS -N MAT_masking
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=12:00:00
#PBS -l mem=80gb
#PBS -M hcm14449@uga.edu
#PBS -m abe

bedtools_module="BEDTools/2.28.0-foss-2018a"
ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/genome.337.fasta"
#directory reference genome is located in
ref_genome_dir="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/"
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas"
to_mask="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/MATlocusToMask.gff3"

module load ${bedtools_module}

bedtools maskfasta -fi ${ref_genome} -bed ${to_mask} -fo ${ref_genome_dir}/337_MATmasked.fasta
