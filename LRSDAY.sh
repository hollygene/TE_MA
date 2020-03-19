#!/bin/bash

run_dir="/scratch/jc33471/pilon/337"
data_dir=$run_dir/data
pilon_dir=$run_dir/pilon
mum_dir=$run_dir/mummer
ragoo_dir=$run_dir/ragoo
anno_dir=$run_dir/annotation
script_dir=$run_dir/scripts
LRSDAY_HOME="/home/jc33471/LRSDAY"
ENV_PATH="/home/jc33471/miniconda/envs/LRSDAY2"

mkdir -p $data_dir $pilon_dir $mum_dir $ragoo_dir $anno_dir $script_dir


prefix="337"
# copy LRSDAY script to script_dir
cp $LRSDAY_HOME/scripts/* $script_dir
cp $LRSDAY_HOME/misc/* $script_dir
cp $LRSDAY_HOME/data/* $data_dir

source activate LRSDAY2


cd $anno_dir
# prepare assembly for 337
cp $ragoo_dir/ragoo_output_sr/ragoo.fasta . # change form corr to sr
awk 'BEGIN{RS=">"} /^Spar_/ {print ">"$0}' ragoo.fasta > $anno_dir/genome.337.fasta
# annotation for YPS138 from Yue et al.
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/YPS138.genome.fa.gz
gunzip YPS138.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/YPS138.all_feature.gff.gz
gunzip YPS138.all_feature.gff.gz
# extract centromere
grep centromere $anno_dir/YPS138.all_feature.gff > $anno_dir/YPS138.centromere.gff
bedtools getfasta -fi $anno_dir/YPS138.genome.fa -bed $anno_dir/YPS138.centromere.gff -fo $anno_dir/YPS138.centromere.fa
# rename fasta headers
sed -i -E 's/>chr(.*):.*/>Centromere\1/g' $anno_dir/YPS138.centromere.fa


# LRSDAY script
query=$anno_dir/YPS138.centromere.fa
genome=$anno_dir/genome.337.fasta
exonerate --showvulgar no --showcigar no --showalignment no --showtargetgff yes --bestn 1 $query $genome > $prefix.centromere.exonerate.gff
perl $script_dir/exonerate_gff2gff3.pl  -i $prefix.centromere.exonerate.gff -o $prefix.centromere.gff3.tmp -t $prefix
perl $script_dir/tidy_maker_gff3.pl -r $genome -i  $prefix.centromere.gff3.tmp -o  $prefix.centromere.gff3 -t $prefix

https://sgd-prod-upload.s3.amazonaws.com/S000208228/SGD_features.README
https://sgd-prod-upload.s3.amazonaws.com/S000211940/SGD_features.20160514.tab.gz



## FROM LRSDAY
#######################################
# set project-specific variables
prefix="SK1" # The file name prefix for the processing sample. Default = "SK1" for the testing example.
genome="./../07.Supervised_Final_Assembly/$prefix.assembly.final.fa" # The path of the input genome assembly.
query="$LRSDAY_HOME/data/S288C.centromere.fa" # The S. cerevisiae S288C reference centromere sequences based on Yue et al. (2017) Nature Genetics.
debug="no" # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".




$exonerate_dir/exonerate --showvulgar no --showcigar no --showalignment no --showtargetgff yes --bestn 1 $query $genome >$prefix.centromere.exonerate.gff
perl $LRSDAY_HOME/scripts/exonerate_gff2gff3.pl  -i $prefix.centromere.exonerate.gff -o $prefix.centromere.gff3.tmp -t $prefix
perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl -r $genome -i  $prefix.centromere.gff3.tmp -o  $prefix.centromere.gff3 -t $prefix





cd ./../10.Core_X_Element_Annotation
bash LRSDAY.10.Core_X_Element_Annotation.sh





#LRSDAY pipeline
