




java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  \
/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar FastqToSam \
    FASTQ=${raw_data}/${BASE}R1.fq.gz \
    FASTQ2=${raw_data}/${BASE}R2.fq.gz \
    OUTPUT=${output_directory}/${BASE}_fastqtosam.bam \
    READ_GROUP_NAME=${BASE} \
    SAMPLE_NAME=${BASE} \
    PLATFORM=illumina \
    SEQUENCING_CENTER=GGBC
