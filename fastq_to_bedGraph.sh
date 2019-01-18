#PBS -S /bin/bash
#PBS -N j_s_samtools
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00

cd $PBS_O_WORKDIR

data_dir="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/samples"

#first need to generate alignment (paired-end) using bwa

module load BWA/0.7.15-foss-2016b

# index reference genome using index
# time bwa index /scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome.fa


# bwa mem /scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome.fa \
# Ancestor_R1.fastq Ancestor_R2.fastq > Ancestor.sam

#list files in the directory without extension
cd $data_dir
# this does not work
#  ls -1 | sed 's/\_[a-z*].[a-z]*//g' >
data_dir="/scratch/hcm14449/TE_MA_Paradoxus/Practice/files/samples"

for file in `ls -d -1 $data_dir/*_R1.fastq`
 do
 sample_name=$(basename "$file" _R1.fastq)
 fq1=$file
 fq2=$(echo ${file}|sed -E "s/_R1/_R2/g")

 bwa mem /scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome.fa \
"${fq1}" "${fq2}" > "${fq1}".sam

 done


#
#  if [ $sample_no -eq 1 ]; then
#     JOBID=`qsub -N ${sample_name} -v fq1="${fq1}",fq2="${fq2}" mcc_unit_call_sapelo2_yeast.sh`
#   else
#   qsub -W depend=afterok:$JOBID -N ${sample_name} -v mcc_dir="${mcc_dir}",out_dir="${out_dir}",ref_dir="${yeast_ref}",te_library="${yeast_te_seqs}",fq1="${fq1}",fq2="${fq2}" -o "${log_output}" -e "${log_output}" mcc_unit_call_sapelo2_yeast.sh
#   fi
#   sample_no=$((sample_no+1))
#
# for file in ./*.fastq
#
# do
#
# FBASE=$(basename $file .fastq)
# BASE=${FBASE%.fastq}
#
# trim_galore --phred33 -q 20 -o trimmed ${BASE}.fastq
#
# done
#
# module load SAMtools/1.6-foss-2016b
#
# samtools
