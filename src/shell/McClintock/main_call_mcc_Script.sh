#!/bin/bash
#PBS -q batch
#PBS -N mcc_yeast
#PBS -l nodes=1:ppn=4 -l mem=20gb
#PBS -l walltime=50:00:00
#PBS -M hcm14449@uga.edu
#PBS -m abe
#PBS -o /lustre1/hcm14449/mcc/mcc_yeast_0823.out
#PBS -e /lustre1/hcm14449/mcc/mcc_yeast_0823.out
#PBS -j oe

cd $PBS_O_WORKDIR

module load SAMtools/0.1.19-foss-2016b
module load parallel/20160622-foss-2016b
module load Trim_Galore/0.4.5-foss-2016b

mcc_dir="/lustre1/hcm14449/mcc/mcc_8_22_18"
run_dir="/lustre1/hcm14449/TE_MA_Paradoxus/Anc_SpikeIns/DNA_copy/8_23_18"
mkdir -p $run_dir
rm -rf $run_dir/*

raw_dir="/lustre1/hcm14449/TE_MA_Paradoxus/Anc_SpikeIns/DNA_copy/data"

data_dir=$run_dir/data
out_dir=$run_dir/out
mkdir -p $data_dir
mkdir -p $out_dir

# ############### Download fastq files from EBI ###############
# # SRR4074413
# cd $data_dir
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/003/SRR4074413/SRR4074413_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/003/SRR4074413/SRR4074413_2.fastq.gz .
# # SRR4074412
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/002/SRR4074412/SRR4074412_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/003/SRR4074413/SRR4074412_2.fastq.gz .
# # SRR4074411
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/001/SRR4074411/SRR4074411_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/003/SRR4074413/SRR4074411_2.fastq.gz .
# # SRR4074394
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/004/SRR4074394/SRR4074394_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/004/SRR4074394/SRR4074394_2.fastq.gz .
# # SRR4074385
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/005/SRR4074385/SRR4074385_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/005/SRR4074385/SRR4074385_2.fastq.gz .
# # SRR4074384
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/004/SRR4074384/SRR4074384_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/004/SRR4074384/SRR4074384_2.fastq.gz .
# # SRR4074383
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/003/SRR4074383/SRR4074383_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/003/SRR4074383/SRR4074383_2.fastq.gz .
# # SRR4074358
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/008/SRR4074358/SRR4074358_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/008/SRR4074358/SRR4074358_2.fastq.gz .
# # SRR4074258
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/008/SRR4074258/SRR4074258_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/008/SRR4074258/SRR4074258_2.fastq.gz .
# # SRR4074257
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/007/SRR4074257/SRR4074257_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/007/SRR4074257/SRR4074257_2.fastq.gz .
# # SRR4074256
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/006/SRR4074256/SRR4074256_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/006/SRR4074256/SRR4074256_2.fastq.gz .
# # SRR4074255
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/005/SRR4074255/SRR4074255_1.fastq.gz .
# ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR407/005/SRR4074255/SRR4074255_2.fastq.gz .


############### Read trimming ###############
# trim reads and unzip to data folder
quality=20
stringency=6
length=50

# trim and unzip raw reads
export quality
export stringency
export length
export data_dir

fq1_list=($(ls -d $raw_dir/*_001.fastq))
# pair-end
# parallel trim_galore --paired --fastqc -q $quality --stringency $stringency -o $data_dir {} {=s/_1/_2/=} ::: "${fq1_list[@]}"

# adapt to single end
parallel trim_galore --fastqc -q $quality --stringency $stringency -o $data_dir {} ::: "${fq1_list[@]}"

# # unzip
# fqz_dir=(`ls $data_dir/*val*.fq.gz`)
# export data_dir
# parallel -j $PBS_NP 'gunzip -c {} > $data_dir/$(basename {.} .gz)' ::: "${fqz_dir[@]}"

# ############### Prepare ref and TE library ###############
# prep yeast ref
if [ ! -f $data_dir/sacCer2.fasta ]
then
wget -P $data_dir http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz
tar xvzf $data_dir/chromFa.tar.gz -C $data_dir
rm $data_dir/chromFa.tar.gz
# Combine the chromosomes together
cat $data_dir/chr*fa $data_dir/2micron.fa > $data_dir/sacCer2.fasta
rm $data_dir/chr*fa $data_dir/2micron.fa
fi
yeast_ref=$data_dir/sacCer2.fasta

# yeast TE library
yeast_te_seqs="/home/hcm14449/Github/mcclintock/test/sac_cer_TE_seqs.fasta"

############### Run McC ###############
# set up first sample that other samples depend on
sample_no=1

# Run the pipeline

pay attention to the file name
# pair end mode
for file in `ls -d -1 $data_dir/*_val_1.fq`
do
sample_name=$(basename "$file" _val_1.fq)
fq1=$file
fq2=$(echo ${file}|sed -E "s/_1/_2/g")
log_output=$out_dir/$sample_name.log
if [ $sample_no -eq 1 ]; then
JOBID=`qsub -N ${sample_name} -v mcc_dir="${mcc_dir}",out_dir="${out_dir}",ref_dir="${yeast_ref}",te_library="${yeast_te_seqs}",fq1="${fq1}",fq2="${fq2}" -o "${log_output}" -e "${log_output}" mcc_unit_call_sapelo2_yeast.sh`
else
qsub -W depend=afterok:$JOBID -N ${sample_name} -v mcc_dir="${mcc_dir}",out_dir="${out_dir}",ref_dir="${yeast_ref}",te_library="${yeast_te_seqs}",fq1="${fq1}",fq2="${fq2}" -o "${log_output}" -e "${log_output}" mcc_unit_call_sapelo2_yeast.sh
fi
sample_no=$((sample_no+1))
done
†
# single end mode
# for file in `ls -d -1 $data_dir/*_001_trimmed.fq`
# do
# sample_name=$(basename "$file" _001_trimmed.fq)
# fq1=$file
# log_output=$out_dir/$sample_name.log
# if [ $sample_no -eq 1 ]; then
# JOBID=`qsub -N ${sample_name} -v mcc_dir="${mcc_dir}",out_dir="${out_dir}",ref_dir="${yeast_ref}",te_library="${yeast_te_seqs}",fq1="${fq1}" -o "${log_output}" -e "${log_output}" unit_call_SE.sh`
# else
# qsub -W depend=afterok:$JOBID -N ${sample_name} -v mcc_dir="${mcc_dir}",out_dir="${out_dir}",ref_dir="${yeast_ref}",te_library="${yeast_te_seqs}",fq1="${fq1}" -o "${log_output}" -e "${log_output}" unit_call_SE.sh
# fi
# sample_no=$((sample_no+1))
# done
