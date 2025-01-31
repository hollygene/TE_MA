#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N MAT_masking
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=36:00:00
#PBS -l mem=200gb
#PBS -M hcm14449@uga.edu
#PBS -m abe

bedtools_module="BEDTools/2.28.0-foss-2018a"
picard_module="picard/2.4.1-Java-1.8.0_144"
bwa_module="BWA/0.7.15-foss-2016b"
samtools_module="SAMtools/1.6-foss-2016b"
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

# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar CreateSequenceDictionary \
#       R=${ref_genome_dir}/337_MATmasked.fasta \
#       O=${ref_genome_dir}/337_MATmasked.dict

# #index masked reference
# bwa index ${ref_genome_dir}/337_MATmasked.fasta
#
# #D0 ancestor

### ancestor spike ins

Anc_SpikeIns="/project/dwhlab/Holly/TE_MA_Paradoxus/Paradoxus_MA/Anc_SpikeIns/Holly_gDNA"


# mkdir ${raw_data}/Anc_SpikeIns
cd ${raw_data}/Anc_SpikeIns
module load ${bwa_module}

module load ${samtools_module}


# for file in ${raw_data}/Anc_SpikeIns/*_sorted.bam
#
# do
# FBASE=$(basename $file _sorted.bam)
# BASE=${FBASE%_sorted.bam}

# bwa aln ${ref_genome} ${raw_data}/Anc_SpikeIns/${BASE}.fastq.gz > ${raw_data}/Anc_SpikeIns/${BASE}.sai
# bwa samse ${ref_genome} ${raw_data}/Anc_SpikeIns/${BASE}.sai ${raw_data}/Anc_SpikeIns/${BASE}.fastq.gz > ${raw_data}/Anc_SpikeIns/${BASE}.sam
# samtools view -b ${raw_data}/Anc_SpikeIns/${BASE}.sam > ${raw_data}/Anc_SpikeIns/${BASE}.bam
# samtools sort ${raw_data}/Anc_SpikeIns/${BASE}.bam > ${raw_data}/Anc_SpikeIns/${BASE}_sorted.bam
# samtools index ${raw_data}/Anc_SpikeIns/${BASE}_sorted.bam
#
# done
#
# for file in ${raw_data}/Anc_SpikeIns/*_sorted.bam
#
# do
# FBASE=$(basename $file _sorted.bam)
# BASE=${FBASE%_sorted.bam}
#
# samtools view -b ${BASE}_sorted.bam Spar_III_RaGOO:210706-211230 > ${BASE}.HML.bam
#
# samtools mpileup ${BASE}.HML.bam | awk '{ count++ ; SUM += $4 } END { print "Total: " SUM "\t" "Nucleotides: " count "\t" "Average_coverage: " SUM/count }' >  ${BASE}_HML_Depth.txt
#
# done

########## HMR
#
# for file in ${raw_data}/Anc_SpikeIns/*_sorted.bam
#
# do
# FBASE=$(basename $file _sorted.bam)
# BASE=${FBASE%_sorted.bam}
#
# samtools view -b ${BASE}_sorted.bam Spar_III_RaGOO:304110-304598 > ${BASE}.HMR.bam
#
# samtools mpileup ${BASE}.HMR.bam | awk '{ count++ ; SUM += $4 } END { print "Total: " SUM "\t" "Nucleotides: " count "\t" "Average_coverage: " SUM/count }' >  ${BASE}_HMR_Depth.txt
#
# done


########################################################
# for file in ${unmapped_bams}/*_markilluminaadapters.bam
#
# do
#
# FBASE=$(basename $file _markilluminaadapters.bam)
# BASE=${FBASE%_markilluminaadapters.bam}
#
# samtools sort ${unmapped_bams}/${BASE}_markilluminaadapters.bam > ${unmapped_bams}/${BASE}_sorted.bam
#
# done
########################################################
# for file in ${unmapped_bams}/*_sorted.bam
#
# do
#
# FBASE=$(basename $file _sorted.bam)
# BASE=${FBASE%_sorted.bam}
#
# samtools index ${unmapped_bams}/${BASE}_sorted.bam
#
# done

########################################################
for file in ${unmapped_bams}/*_MATmasked.bam

do

FBASE=$(basename $file _MATmasked.bam)
BASE=${FBASE%_MATmasked.bam}

samtools view -b ${unmapped_bams}/${BASE}_MATmasked.bam Spar_III_RaGOO:210706-211230 > ${unmapped_bams}/${BASE}_HML.bam

done
########################################################
for file in ${unmapped_bams}/*_HML.bam

do

FBASE=$(basename $file _HML.bam)
BASE=${FBASE%_HML.bam}

samtools mpileup ${unmapped_bams}/${BASE}_HML.bam | awk '{ count++ ; SUM += $4 } END { print "Total: " SUM "\t" "Nucleotides: " count "\t" "Average_coverage: " SUM/count }' >  ${unmapped_bams}/${BASE}_HML_Depth.txt

done
########################################################
grep ^ ${unmapped_bams}/*_HML_Depth.txt > ${unmapped_bams}/HML_depth.txt



#######


########################################################
for file in ${unmapped_bams}/*_MATmasked.bam

do

FBASE=$(basename $file _MATmasked.bam)
BASE=${FBASE%_MATmasked.bam}

samtools view -b ${unmapped_bams}/${BASE}_MATmasked.bam Spar_III_RaGOO:304110-304598 > ${unmapped_bams}/${BASE}_HMR.bam

done
########################################################
for file in ${unmapped_bams}/*_HMR.bam

do

FBASE=$(basename $file _HMR.bam)
BASE=${FBASE%_HMR.bam}

samtools mpileup ${unmapped_bams}/${BASE}_HMR.bam | awk '{ count++ ; SUM += $4 } END { print "Total: " SUM "\t" "Nucleotides: " count "\t" "Average_coverage: " SUM/count }' >  ${unmapped_bams}/${BASE}_HMR_Depth.txt

done
########################################################


grep ^ ${unmapped_bams}/*_HMR_Depth.txt > ${unmapped_bams}/HMR_depth.txt




# for file in ${unmapped_bams}/*_markilluminaadapters.bam
#
# do
#
# FBASE=$(basename $file _markilluminaadapters.bam)
# BASE=${FBASE%_markilluminaadapters.bam}
# OUT="${BASE}_pipedNewRef.sh"
# echo "#!/bin/bash" >> ${OUT}
# echo "#PBS -N ${BASE}_pipedNewRef" >> ${OUT}
# echo "#PBS -l walltime=12:00:00" >> ${OUT}
# echo "#PBS -l nodes=1:ppn=1:HIGHMEM" >> ${OUT}
# echo "#PBS -q highmem_q" >> ${OUT}
# echo "#PBS -l mem=150gb" >> ${OUT}
# echo "" >> ${OUT}
# echo "cd ${unmapped_bams}" >> ${OUT}
# echo "module load ${picard_module}" >> ${OUT}
# echo "module load ${bwa_module}" >> ${OUT}
# echo "module load ${samtools_module}" >> ${OUT}
# echo "module load ${GATK_module}" >> ${OUT}
# echo "" >> ${OUT}
# echo "java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
# I=${unmapped_bams}/${BASE}_markilluminaadapters.bam \
# FASTQ=/dev/stdout \
# CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
# TMP_DIR=${unmapped_bams}/TMP | \
# bwa mem -M -t 7 -p ${ref_genome_dir}/337_MATmasked.fasta /dev/stdin| \
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
# ALIGNED_BAM=/dev/stdin \
# UNMAPPED_BAM=${unmapped_bams}/${BASE}_markilluminaadapters.bam \
# OUTPUT=${unmapped_bams}/${BASE}337_MATmasked.bam \
# R=${ref_genome_dir}/337_MATmasked.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true \
# CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
# INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
# PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
# TMP_DIR=${unmapped_bams}/TMP" >> ${OUT}
# qsub ${OUT}
#
# done


#H0 ancestor

# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
# I=${H0_ancestor} \
# FASTQ=/dev/stdout \
# CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
# TMP_DIR=${unmapped_bams}/TMP | \
# bwa mem -M -t 7 -p ${ref_genome_dir}/337_MATmasked.fasta /dev/stdin| \
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
# ALIGNED_BAM=/dev/stdin \
# UNMAPPED_BAM=${H0_ancestor} \
# OUTPUT=${unmapped_bams}/H0_ancestor_MATmasked.bam \
# R=${ref_genome_dir}/337_MATmasked.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true \
# CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
# INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
# PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
# TMP_DIR=${unmapped_bams}/TMP

#D1 ancestor

# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
# I=${D1_ancestor} \
# FASTQ=/dev/stdout \
# CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
# TMP_DIR=${unmapped_bams}/TMP | \
# bwa mem -M -t 7 -p ${ref_genome_dir}/337_MATmasked.fasta /dev/stdin| \
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
# ALIGNED_BAM=/dev/stdin \
# UNMAPPED_BAM=${D1_ancestor} \
# OUTPUT=${unmapped_bams}/D1_ancestor_MATmasked.bam \
# R=${ref_genome_dir}/337_MATmasked.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true \
# CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
# INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
# PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
# TMP_DIR=${unmapped_bams}/TMP
#
# #D20 ancestor
#
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar SamToFastq \
# I=${D20_ancestor} \
# FASTQ=/dev/stdout \
# CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
# TMP_DIR=${unmapped_bams}/TMP | \
# bwa mem -M -t 7 -p ${ref_genome_dir}/337_MATmasked.fasta /dev/stdin| \
# java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144" -jar  \
# /usr/local/apps/eb/picard/2.4.1-Java-1.8.0_144/picard.jar MergeBamAlignment \
# ALIGNED_BAM=/dev/stdin \
# UNMAPPED_BAM=${D20_ancestor} \
# OUTPUT=${unmapped_bams}/D20_ancestor_MATmasked.bam \
# R=${ref_genome_dir}/337_MATmasked.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true \
# CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
# INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
# PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
# TMP_DIR=${unmapped_bams}/TMP
