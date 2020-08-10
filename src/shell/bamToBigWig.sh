

module load BEDTools/2.28.0-foss-2018a

module load Python/3.5.2-foss-2016b

module load SAMtools/1.6-foss-2016b

export PATH=${PATH}:/scratch/hcm14449/TE_MA_Paradoxus/jbscripts

python3 /scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py -sort /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa.fai H0_37_New_S23_R1_001.sorted.bam


## Loop
for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/*.sorted.bam

do

FBASE=$(basename $file .sorted.bam)
BASE=${FBASE%.sorted.bam}

python3 /scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py -sort /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa.fai ${BASE}.sorted.bam

done

### for ancestors
for file in /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_gDNA/*.sorted.bam

do

FBASE=$(basename $file .sorted.bam)
BASE=${FBASE%.sorted.bam}

python3 /scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py -sort /scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/YPS138.genome.fa.fai ${BASE}.sorted.bam

done


### for arabidopsis lines


python3 /scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py -sort /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas.fai /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/Col_200_S24_R1_001.sorted.bam

python3 /scratch/hcm14449/TE_MA_Paradoxus/jbscripts/file_to_bigwig_pe.py -sort /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/TAIR10_chr_all.fas.fai /scratch/hcm14449/TE_MA_Paradoxus/test_Spike_InsJune2019/Holly_spikein_289/Holly_spikein/Arabidopsis/Col_500_S18_R1_001.sorted.bam
