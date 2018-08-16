#!/bin/bash
#PBS -q highmem_q
#PBS -N MCL_D20
#PBS -l nodes=1:ppn=4 -l mem=80gb
#PBS -l walltime=24:00:00
#PBS -M hcm14449@uga.edu
#PBS -m abe
#PBS -o /lustre1/hcm14449/TE_MA_Paradoxus/MCL_D20.out
#PBS -e /lustre1/hcm14449/TE_MA_Paradoxus/MCL_D20.out
#PBS -j oe

export PATH=$PATH:/home/sh60271/app/bedtools2/bin
export PATH=/home/sh60271/app/bwa-0.7.4:$PATH
export PATH=$PATH:/home/sh60271/app/RepeatMasker
module load R/3.4.4-foss-2016b-X11-20160819-GACRC
module load BioPerl/1.7.1-foss-2016b-Perl-5.24.1
module load SAMtools/0.1.19-foss-2016b
module load BLAT/3.5-foss-2016b
module load Bowtie/1.2.2-foss-2016b
module load ucsc/359
module load Exonerate/2.4.0-foss-2016b
module load FastQC/0.11.5-Java-1.8.0_144

# module load seqtk

mcc_dir="/lustre1/hcm14449/mcc/mcc8_16_18"
run_dir="/lustre1/hcm14449/TE_MA_Paradoxus/Anc_SpikeIns/DNA_copy"
mkdir -p $run_dir

#data_dir=$run_dir/data
out_dir=$run_dir/out
mkdir -p $out_dir
rm -rf $out_dir/*
#mkdir -p $data_dir
#mkdir -p $out_dir

# unzip fastq read to data folder
#fq1="/lustre1/hcm14449/TE_MA_Paradoxus/Holly_gDNA/HM_H0_S16_R1_001.fastq"

# fq1=$data_dir/$(basename "$fqz1" .gz)
# zcat $fqz1 > $fq1
# fq2=$data_dir/$(basename "$fqz2" .gz)
# zcat $fqz2 > $fq2

#fq1=$data_dir/HM_H0_S16_R1_001.fastq


# Download the reference genome from UCSC (allows easy browsing of results)
#printf "Downloading reference genome...\n\n"

#wget -P $data_dir -nc -q http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz
#tar xvzf $data_dir/chromFa.tar.gz
#rm $data_dir/chromFa.tar.gz
# Combine the chromosomes together
#cat /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/*.1 >  /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/Spar.ref.fa

#ref_dir=/lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa

# Download gff locations of reference TE copies
#wget -P $data_dir -nc -q http://files.figshare.com/287395/File_S2.txt
#awk '{print $3"\treannotate\ttransposable_element\t"$4"\t"$5"\t.\t"$6"\t.\tID="$1;}' $data_dir/File_S2.txt > tmp
#sed '1d;$d' $data_dir/tmp > $data_dir/reference_TE_locations.gff
#rm $data_dir/File_S2.txt
#rm $data_dir/tmp

# TE database
#te_seqs_dir="/home/sh60271/git/transposons/current/D_mel_transposon_sequence_set.fa"
# Remove anything after the # so that forward slashes are not included
#awk -F"#" '{if($0 ~ ">") print $1; else print $0}' $te_seqs_dir > $data_dir/tmp
#mv $data_dir/tmp $data_dir/$(basename "$te_seqs_dir")
#te_seqs_dir=$data_dir/$(basename "$te_seqs_dir")

# Run McC pipeline


#bash $mcc_dir/mcclintock.sh -d -C -o $out_dir -m "retroseq temp ngs_te_mapper te-locate relocate" \
#-r $ref_dir -c $te_seqs_dir -1 $fq1 -p $PBS_NP -b

#sh mcclintock.sh -m "RelocaTE ngs_te_mapper" -r reference.fasta -c te_consensus.fasta -g te_locations.gff -t te_families.tsv \
#-1 sample_1.fastq -2 sample_2.fastq -p 2 -i -b
bash $mcc_dir/mcclintock.sh -m "relocate temp ngs_te_mapper retroseq te-locate" \
-r /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/cerevisiae/sacCer2.fasta -c /home/hcm14449/Github/mcclintock/test/sac_cer_TE_seqs.fasta  \
-t /home/hcm14449/Github/mcclintock/test/sac_cer_te_families.tsv \
-1 /lustre1/hcm14449/TE_MA_Paradoxus/Anc_SpikeIns/DNA_copy/HM_D20_S15_R1_001.fastq -o $out_dir \
-C -d -b



#bash /lustre1/hcm14449/mcc/mcc8_18_18/mcclintock.sh -m "relocate temp ngs_te_mapper retroseq te-locate" \
#-r /lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/cerevisiae/sacCer2.fasta -c /home/hcm14449/Github/mcclintock/test/sac_cer_TE_seqs.fasta  \
#-t /home/hcm14449/Github/mcclintock/test/sac_cer_te_families.tsv \
#-1 /lustre1/hcm14449/TE_MA_Paradoxus/Anc_SpikeIns/DNA_copy/HM_D20_S15_R1_001.fastq -o /lustre1/hcm14449/TE_MA_Paradoxus/Anc_SpikeIns/DNA_copy/out \
#-C -d -b

#-pa 4
