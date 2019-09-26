module load FastQC/0.11.8-Java-1.8.0_144


time fastqc -o /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/QC_Out /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/HM-H0-A_R2_001.fastq


muver_module="muver/0.1.0-foss-2016b-Python-2.7.14-20190318"

module load muver/0.1.0-foss-2016b-Python-2.7.14-20190318
# module load ${GATK_module}
# java -jar $EBROOTGATK/GenomeAnalysisTK.jar

#index reference genome
# echo ${ref_genome}
muver index-reference /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa

# muver index-reference S288C_reference_sequence_R64-2-1_20150113.fa
#
# create repeat file for reference genome
muver create-repeat-file scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa


muver run-pipeline /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa /home/hcm14449/Github/TE_MA/FASTQ_LIST_SHORT_H0 ${control_sample_name_H0} ${experiment_directory_H0}
