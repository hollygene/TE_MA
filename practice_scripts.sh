# for practicing on ancestor for H0



module load ${samtools_module}

#index reference genome

samtools faidx ${ref_genome}

#convert sam files to bam files
for file in ${output_directory}/*_aln.sam

do

FBASE=$(basename $file _aln.sam)
BASE=${FBASE%_aln.sam}

samtools view -bt ${ref_genome_dir}/*.fai \
${output_directory}/${BASE}_aln.sam \
  > ${output_directory}/${BASE}.bam

done
