
#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=5 -l mem=40gb
#PBS -l walltime=80:00:00
#PBS -M hcm14449@uga.edu
#PBS -m abe
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




cd $PBS_O_WORKDIR

bash $mcc_dir/mcclintock.sh -d -C -o $out_dir -r $ref_dir -c $te_library -1 $fq1 -p $PBS_NP -b
