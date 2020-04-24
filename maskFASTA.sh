#PBS -S /bin/bash
#PBS -q batch
#PBS -N MAT_masking
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=12:00:00
#PBS -l mem=80gb
#PBS -M hcm14449@uga.edu
#PBS -m abe

bedtools_module="BEDTools/2.28.0-foss-2018a"
picard_module="picard/2.4.1-Java-1.8.0_144"
bwa_module="BWA/0.7.15-foss-2016b"
ref_genome="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/genome.337.fasta"
#directory reference genome is located in
ref_genome_dir="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/"
raw_data="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas"
to_mask="/scratch/hcm14449/TE_MA_Paradoxus/ref_genome/paradoxus/337Ref/MATlocusToMask.gff3"
D0_ancestor="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas/unmapped_bams/D0-A__markilluminaadapters.bam"
H0_ancestor="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas/unmapped_bams/H0-A__markilluminaadapters.bam"
D1_ancestor="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas/unmapped_bams/D1-A__markilluminaadapters.bam"
D20_ancestor="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas/unmapped_bams/D20-A__markilluminaadapters.bam"
unmapped_bams="/scratch/hcm14449/TE_MA_Paradoxus/Illumina_Data/AllFastas/unmapped_bams"
module load ${bedtools_module}
module load ${picard_module}
module load ${bwa_module}
# bedtools maskfasta -fi ${ref_genome} -bed ${to_mask} -fo ${ref_genome_dir}/337_MATmasked.fasta

###################
####
## Use piped command like in normal pipeline to map to masked reference genome
#create sequence dictionary for masked genome

java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar CreateSequenceDictionary \
      R=${ref_genome_dir}/337_MATmasked.fasta \
      O=${ref_genome_dir}/337_MATmasked.dict


#D0 ancestor

java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
I=${D0_ancestor} \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=${unmapped_bams}/TMP | \
bwa mem -M -t 7 -p ${ref_genome_dir}/337_MATmasked.fasta /dev/stdin| \
java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=${D0_ancestor} \
OUTPUT=${unmapped_bams}/D0_ancestor_MATmasked.bam \
R=${ref_genome_dir}/337_MATmasked.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=${unmapped_bams}/TMP
