awk scripts for parsing .gff3 file(s) for use with karyoploteR

awk '$3 == "expressed_sequence_match" {print $0}' 337.nuclear_genome.est_evidence.txt > 337_esm.txt
awk '{print $1, $4, $5, $9}' 337_esm.txt > 337_esm_red.txt

awk -F"=" '$1=$1' OFS="\t" 337_esm_red.txt > 337_split_esm.txt

awk '{print $1, $2, $3, $6}' 337_split_esm.txt > 337_sites.bed

awk '$1 ~ "##" {print $0}' 337.centromere.gff3 > 337.chrms.bed

awk '{print $2, $3, $4}' 337.chrms.bed > 337_chrm_sizes.bed

awk 'length {$2 ~ /1/}' 337_chrm_sizes.bed > 337_chrm_sizes_red.bed

awk -F "\t" 'NF > 1 {print $0}' 337_chrm_sizes.bed > 337_chrm_sizes_red2.bed

awk 'NF' 337_chrm_sizes_red.bed > 337_chrm_sizes_red2.bed


awk -F" " 'BEGIN { OFS = "\t" } {print $0}' 337_sites.bed > 337_sites_2.bed

sed $'s/ /\t/g' 337_sites.bed > 337_sites_2.bed


for file in /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/BedGraphs/TELocate/D20/raw/*_raw.TY1.bed 

do

FBASE=$(basename $file _raw.TY1.bed)
BASE=${FBASE%_raw.TY1.bed}

awk 'BEGIN{FS=OFS="\t"}{$1=$1"_RaGOO"}1' < ${BASE}_raw.TY1.bed > ${BASE}.bed

done



for file in ${raw_data}/*_R1_001.fastq

do

FBASE=$(basename $file _R1_001.fastq)
BASE=${FBASE%_R1_001.fastq}
