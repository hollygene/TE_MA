#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N genotype
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=72:00:00
#PBS -l mem=200gb
#PBS -M hcm14449@uga.edu
#PBS -m abe

# Variant effect predictor (VEP)
# https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff
# https://wiki.gacrc.uga.edu/wiki/Ensembl-VEP-Sapelo2
# https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html

HTSlib_module="HTSlib/1.8-foss-2018a"
gff_file="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/337.nuclear_genome.est_evidence.gff3"
workDir="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/"
vcf_file="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/H0_noLow_noHigh_redGem_AncCalls_NoHetsVars.vcf"
out_file="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/vep.test.txt"
genome_fasta="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/genome.337.fasta"

module load ${HTSlib_module}
module load gffread/0.9.12-foss-2016b

# grep -v "#" ${gff_file} | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > ${workDir}/337.nuclear_genome.est_evidence.gff3.gz
# tabix -p gff /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/337.nuclear_genome.est_evidence.gff3.gz
# bgzip -c ${genome_fasta} > /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/genome.337.fasta.gz
#
# # need to make a gtf file from gff3 file for use in cache building for VEP
# gffread ${workDir}/337.nuclear_genome.maker.raw.gff3 -T -o ${workDir}/337.nuclear_genome.maker.raw.gtf

# building the cache

# module load VEP/95.0-foss-2018b-Perl-5.28.0

# doing this on my computer now
perl gtf2vep.pl -i /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/337.nuclear_genome.maker.raw.gtf -f /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/genome.337.fasta -d 78 -s s_par_337


./vep -i /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/VCFs/H0/AllVariants/H0_noLow_noHigh_redGem_AncCalls_NoHetsVars.vcf --cache --gtf /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/337.nuclear_genome.maker.raw.gtf --fasta /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/genome.337.fasta




perl variant_effect_predictor.pl -offline -i ${vcf_file} -s s_par_337





singularity exec /usr/local/singularity-images/ensembl-vep.simg vep -i ${vcf_file} -o ${out_file} --gff ${workDir}/337.nuclear_genome.est_evidence.gff3.gz --fasta /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/genome.337.fasta.gz

singularity exec /usr/local/singularity-images/ensembl-vep.simg/ensembl-vep.simg which vep


# singularity exec ./ensembl-vep.simg vep [options]
