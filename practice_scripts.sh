# for practicing on ancestor for H0


bwa_module="BWA/0.7.15-foss-2016b"
#location of samtools module
samtools_module="SAMtools/1.6-foss-2016b"
#location of bedtools module
bedtools_module="BEDTools/2.28.0-foss-2018a"
#location of python module
python_module="Python/3.5.2-foss-2016b"
#location of picard module
picard_module="picard/2.16.0-Java-1.8.0_144"
#location of GATK module
GATK_module="GATK/4.0.3.0-Java-1.8.0_144"
#location of bamtoBigWig script and accessories
script_location="/scratch/hcm14449/TE_MA_Paradoxus/jbscripts"
#location of bam to bigwig script
bamToBigWig="/scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py"
module load SAMtools/1.6-foss-2016b

#index reference genome

samtools faidx /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa

samtools view -bt /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/*.fai \
/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/H0/HM-H0-A_aln.sam \
  > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/H0/HM-H0-A.bam

  samtools sort -o /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/H0/HM-H0-A.sorted.bam \
     /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/H0/HM-H0-A.bam


     module load BEDTools/2.28.0-foss-2018a
     module load Python/3.5.2-foss-2016b
     module load SAMtools/1.6-foss-2016b
     export PATH=${PATH}:/scratch/hcm14449/TE_MA_Paradoxus/jbscripts

     python3 /scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py -sort /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/*.fai /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/H0/HM-H0-A.sorted.bam

     samtools depth \
     /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/H0/HM-H0-A.sorted.bam \
     |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/trimmed/H0/HM-H0-A_depth.txt
###################################################################################
     module load SAMtools/1.6-foss-2016b

     #index reference genome
module load BWA/0.7.15-foss-2016b
mkdir /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out

     # samtools faidx /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa
bwa mem -M /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/HM-D20-A_R1_001.fastq /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/HM-D20-A_R2_001.fastq > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D20-A_aln.sam

     samtools view -bt /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/*.fai \
     /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D20-A_aln.sam \
       > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D20-A.bam

       samtools sort -o /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D20-A.sorted.bam \
          /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D20-A.bam


          module load BEDTools/2.28.0-foss-2018a
          module load Python/3.5.2-foss-2016b
          module load SAMtools/1.6-foss-2016b
          export PATH=${PATH}:/scratch/hcm14449/TE_MA_Paradoxus/jbscripts

          python3 /scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py -sort /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/*.fai /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D20-A.sorted.bam

          samtools depth \
          /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D20-A.sorted.bam \
          |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/HM-D20-A_depth.txt

          ###################################################################################
    module load SAMtools/1.6-foss-2016b

               #index reference genome
          module load BWA/0.7.15-foss-2016b
          # mkdir /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out

               # samtools faidx /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa
          bwa mem -M /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/HM-D1-A_R1_001.fastq /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/IL_Data/GW_run3/00_fastq/HM-D1-A_R2_001.fastq > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D1-A_aln.sam

               samtools view -bt /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/*.fai \
               /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D1-A_aln.sam \
                 > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D1-A.bam

                 samtools sort -o /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D1-A.sorted.bam \
                    /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D1-A.bam


                    module load BEDTools/2.28.0-foss-2018a
                    module load Python/3.5.2-foss-2016b
                    module load SAMtools/1.6-foss-2016b
                    export PATH=${PATH}:/scratch/hcm14449/TE_MA_Paradoxus/jbscripts

                    python3 /scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py -sort /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/*.fai /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D1-A.sorted.bam

                    samtools depth \
                    /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/HM-D1-A.sorted.bam \
                    |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/HM-D1-A_depth.txt



                    for file in /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/*_aln.sam

                    do

                    FBASE=$(basename $file _aln.sam)
                    BASE=${FBASE%_aln.sam}

                    samtools view -bt /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa/*.fai \
                    /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/${BASE}_aln.sam \
                      > /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/${BASE}.bam

                    done
for file in /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/*.sam

do

  FBASE=$(basename $file .sam)
BASE=${FBASE%.sam}

samtools flagstat /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/${BASE}.sam \
> /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D0/${BASE}_stats.txt

done


module load SAMtools/1.6-foss-2016b

samtools flagstat /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-A.bam \
> /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-A.bam_stats.txt

samtools flagstat /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-4.bam \
> /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/H0/HM-H0-4.bam_stats.txt

samtools flagstat /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-A.bam \
> /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-A.bam_stats.txt

samtools flagstat /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-22.bam \
> /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-22.bam_stats.txt

samtools flagstat /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-7.bam \
> /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D1/HM-D1-7.bam_stats.txt

samtools flagstat /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D20/HM-D20-44.bam \
> /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D20/HM-D20-44.bam_stats.txt

samtools flagstat /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D20/HM-D20-11.bam \
> /scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/Out/D20/HM-D20-11.bam_stats.txt
