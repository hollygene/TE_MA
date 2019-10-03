## Script to create features file of hard-to-map regions to remove from analysis
#1
# change .gff file into a txt file
mv /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138.all_feature.gff /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138.all_feature.txt
# #extract certain phrases from $3
awk '{if ($3 ~ /centromere/ || $3 ~ /telomere/ || $3 ~ /X_element/ || $3 ~ /long_terminal_repeat/ || $3 ~ /Y_prime_element/ || $3 ~ /LTR/ || $3 ~ /TY/) print $0}' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138.all_feature.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_hard_to_map.txt
#
# # end up with 445 line file of features to remove from analysis
#
# # now want to combine subtelomere file with this features file
# # change .gff file to .txt
mv /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138.subtelomere.gff /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138.subtelomere.txt
#
# # concatenate subtelomere file with hard to map file
cat /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138.subtelomere.txt /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_hard_to_map.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove.txt
#
# # want to remove any lines where $1 == $1 and $2 > $4 in other file and $2 < $5 in other file
# #11
#
# # add header to features file
echo -e "Chromosome\tStrain\tFeature\tStart\tStop\tDot\tStrand\tDot\tInfo" | cat - /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove_header.txt
# #
awk -F'\t' '
FNR == 1  {print; next}   # file1 header
{
  if ($1=="chrI") print "ref|NC_001133|","\t",$0;
  if ($1=="chrII") print "ref|NC_001134|","\t",$0;
  if ($1=="chrIII") print "ref|NC_001135|","\t",$0;
  if ($1=="chrIV") print "ref|NC_001136|","\t",$0;
  if ($1=="chrV") print "ref|NC_001137|","\t",$0;
  if ($1=="chrVI") print "ref|NC_001138|","\t",$0;
  if ($1=="chrVII") print "ref|NC_001139|","\t",$0;
  if ($1=="chrVIII") print "ref|NC_001140|","\t",$0;
  if ($1=="chrIX") print "ref|NC_001141|","\t",$0;
  if ($1=="chrX") print "ref|NC_001142|","\t",$0;
  if ($1=="chrXI") print "ref|NC_001143|","\t",$0;
  if ($1=="chrXII") print "ref|NC_001144|","\t",$0;
  if ($1=="chrXIII") print "ref|NC_001145|","\t",$0;
  if ($1=="chrXIV") print "ref|NC_001146|","\t",$0;
  if ($1=="chrXV") print "ref|NC_001147|","\t",$0;
  if ($1=="chrXVI") print "ref|NC_001148|","\t",$0;
}' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove_header.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove_header_Chrs.txt

# echo -e "Chromosome\tStrain\tFeature\tStart\tStop\tDot\tStrand\tDot\tInfo" | cat - YPS138_seq_to_remove_header_Chrs.txt > YPS138_seq_to_remove_header_Chrs_header.txt

# sort -t $'\t' -k 1,5 YPS138_seq_to_remove_header_Chrs.txt > YPS138_seq_to_remove_header_Chrs_sort.txt

# for mutations file too
awk -F'\t' '
FNR == 1  {print; next}   # file1 header
{
  if ($1=="chrI") print "ref|NC_001133|","\t",$0;
  if ($1=="chrII") print "ref|NC_001134|","\t",$0;
  if ($1=="chrIII") print "ref|NC_001135|","\t",$0;
  if ($1=="chrIV") print "ref|NC_001136|","\t",$0;
  if ($1=="chrV") print "ref|NC_001137|","\t",$0;
  if ($1=="chrVI") print "ref|NC_001138|","\t",$0;
  if ($1=="chrVII") print "ref|NC_001139|","\t",$0;
  if ($1=="chrVIII") print "ref|NC_001140|","\t",$0;
  if ($1=="chrIX") print "ref|NC_001141|","\t",$0;
  if ($1=="chrX") print "ref|NC_001142|","\t",$0;
  if ($1=="chrXI") print "ref|NC_001143|","\t",$0;
  if ($1=="chrXII") print "ref|NC_001144|","\t",$0;
  if ($1=="chrXIII") print "ref|NC_001145|","\t",$0;
  if ($1=="chrXIV") print "ref|NC_001146|","\t",$0;
  if ($1=="chrXV") print "ref|NC_001147|","\t",$0;
  if ($1=="chrXVI") print "ref|NC_001148|","\t",$0;
}' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D20/mutations.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D20/D20_mutations_Chrs.txt

awk -F'\t' '
FNR == 1  {print; next}   # file1 header
{
  if ($1=="chrI") print "ref|NC_001133|","\t",$0;
  if ($1=="chrII") print "ref|NC_001134|","\t",$0;
  if ($1=="chrIII") print "ref|NC_001135|","\t",$0;
  if ($1=="chrIV") print "ref|NC_001136|","\t",$0;
  if ($1=="chrV") print "ref|NC_001137|","\t",$0;
  if ($1=="chrVI") print "ref|NC_001138|","\t",$0;
  if ($1=="chrVII") print "ref|NC_001139|","\t",$0;
  if ($1=="chrVIII") print "ref|NC_001140|","\t",$0;
  if ($1=="chrIX") print "ref|NC_001141|","\t",$0;
  if ($1=="chrX") print "ref|NC_001142|","\t",$0;
  if ($1=="chrXI") print "ref|NC_001143|","\t",$0;
  if ($1=="chrXII") print "ref|NC_001144|","\t",$0;
  if ($1=="chrXIII") print "ref|NC_001145|","\t",$0;
  if ($1=="chrXIV") print "ref|NC_001146|","\t",$0;
  if ($1=="chrXV") print "ref|NC_001147|","\t",$0;
  if ($1=="chrXVI") print "ref|NC_001148|","\t",$0;
}' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/H0/mutations.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/H0/H0_mutations_Chrs.txt

awk -F'\t' '
FNR == 1  {print; next}   # file1 header
{
  if ($1=="chrI") print "ref|NC_001133|","\t",$0;
  if ($1=="chrII") print "ref|NC_001134|","\t",$0;
  if ($1=="chrIII") print "ref|NC_001135|","\t",$0;
  if ($1=="chrIV") print "ref|NC_001136|","\t",$0;
  if ($1=="chrV") print "ref|NC_001137|","\t",$0;
  if ($1=="chrVI") print "ref|NC_001138|","\t",$0;
  if ($1=="chrVII") print "ref|NC_001139|","\t",$0;
  if ($1=="chrVIII") print "ref|NC_001140|","\t",$0;
  if ($1=="chrIX") print "ref|NC_001141|","\t",$0;
  if ($1=="chrX") print "ref|NC_001142|","\t",$0;
  if ($1=="chrXI") print "ref|NC_001143|","\t",$0;
  if ($1=="chrXII") print "ref|NC_001144|","\t",$0;
  if ($1=="chrXIII") print "ref|NC_001145|","\t",$0;
  if ($1=="chrXIV") print "ref|NC_001146|","\t",$0;
  if ($1=="chrXV") print "ref|NC_001147|","\t",$0;
  if ($1=="chrXVI") print "ref|NC_001148|","\t",$0;
}' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D0/mutations.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D0/D0_mutations_Chrs.txt

awk -F'\t' '
FNR == 1  {print; next}   # file1 header
{
  if ($1=="chrI") print "ref|NC_001133|","\t",$0;
  if ($1=="chrII") print "ref|NC_001134|","\t",$0;
  if ($1=="chrIII") print "ref|NC_001135|","\t",$0;
  if ($1=="chrIV") print "ref|NC_001136|","\t",$0;
  if ($1=="chrV") print "ref|NC_001137|","\t",$0;
  if ($1=="chrVI") print "ref|NC_001138|","\t",$0;
  if ($1=="chrVII") print "ref|NC_001139|","\t",$0;
  if ($1=="chrVIII") print "ref|NC_001140|","\t",$0;
  if ($1=="chrIX") print "ref|NC_001141|","\t",$0;
  if ($1=="chrX") print "ref|NC_001142|","\t",$0;
  if ($1=="chrXI") print "ref|NC_001143|","\t",$0;
  if ($1=="chrXII") print "ref|NC_001144|","\t",$0;
  if ($1=="chrXIII") print "ref|NC_001145|","\t",$0;
  if ($1=="chrXIV") print "ref|NC_001146|","\t",$0;
  if ($1=="chrXV") print "ref|NC_001147|","\t",$0;
  if ($1=="chrXVI") print "ref|NC_001148|","\t",$0;
}' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D1/mutations.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D1/D1_mutations_Chrs.txt


# count how many columns there are
# head -1 mutations_Chrs.txt | wc -w
#
#
# #6
# awk -F'\t' '{print $675,$0}' mutations_Chrs.txt > YPS138_seq_to_remove_header_Chrs_2.txt
#
# sed 's/^\<chrI\>/"chrOne"/g' YPS138_seq_to_remove_header.txt > YPS138_seq_to_remove_header_renamedChrs.txt
#
# awk '{sub(/\chrI\>/, "chrOne")}1' YPS138_seq_to_remove_header.txt > YPS138_seq_to_remove_header_renamedChrs.txt

# awk '
#     NR  == 1  {next}          # file2 header
#     FNR == 1  {print; next}   # file1 header
#     FNR == NR {min[FNR]=$2; max[FNR]=$3; next}
#     {
#         for (key in min)
#             if (min[key] <= $2 && $2 <= max[key])
#                 next
#         print
#     }'

awk -F'\t' 'FNR == 1  {print; next}   # file1 header
{print $1,"\t",$5,"\t",$6}' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove_header_Chrs.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove_header_Chrs_short.txt

# sort -t $'\t' -k1 -k2 YPS138_seq_to_remove_header_short.txt > YPS138_seq_to_remove_header_short_sort.txt

# awk '
#     NR  == 1  {next}          # file2 header
#     FNR == 1  {print; next}   # file1 header
#     FNR == NR {min[FNR]=$2; max[FNR]=$3; next}
#     {
#         for (key in min)
#             if (min[key] <= $2 && $2 <= max[key])
#                 next
#         print
#     }' YPS138_seq_to_remove_header_short_sort.txt mutations.txt > mutations_clean_TEST2.txt
#
awk '
        NR  == 1  {next}          # file2 header
        FNR == 1  {print; next}   # file1 header
        FNR == NR {min[FNR]=$2; max[FNR]=$3; next}
        {
            for (key in min)
                if (min[key] <= $3 && $3 <= max[key])
                    next
            print
        }' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove_header_Chrs_short.txt /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D20/D20_mutations_Chrs.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D20/D20_mutations_Chrs_Clean.txt

        awk '
                NR  == 1  {next}          # file2 header
                FNR == 1  {print; next}   # file1 header
                FNR == NR {min[FNR]=$2; max[FNR]=$3; next}
                {
                    for (key in min)
                        if (min[key] <= $3 && $3 <= max[key])
                            next
                    print
                }' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove_header_Chrs_short.txt /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/H0/H0_mutations_Chrs.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/H0/H0_mutations_Chrs_Clean.txt

                awk '
                        NR  == 1  {next}          # file2 header
                        FNR == 1  {print; next}   # file1 header
                        FNR == NR {min[FNR]=$2; max[FNR]=$3; next}
                        {
                            for (key in min)
                                if (min[key] <= $3 && $3 <= max[key])
                                    next
                            print
                        }' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove_header_Chrs_short.txt /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D0/D0_mutations_Chrs.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D0/D0_mutations_Chrs_Clean.txt

                        awk '
                                NR  == 1  {next}          # file2 header
                                FNR == 1  {print; next}   # file1 header
                                FNR == NR {min[FNR]=$2; max[FNR]=$3; next}
                                {
                                    for (key in min)
                                        if (min[key] <= $3 && $3 <= max[key])
                                            next
                                    print
                                }' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove_header_Chrs_short.txt /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D1/D1_mutations_Chrs.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D1/D1_mutations_Chrs_Clean.txt


        # awk -F'\t' '{
        #   if ($1=="chrI") print "ref|NC_001133|","\t",$0;
        #   if ($1=="chrII") print "ref|NC_001134|","\t",$0;
        #   if ($1=="chrIII") print "ref|NC_001135|","\t",$0;
        #   if ($1=="chrIV") print "ref|NC_001136|","\t",$0;
        #   if ($1=="chrV") print "ref|NC_001137|","\t",$0;
        #   if ($1=="chrVI") print "ref|NC_001138|","\t",$0;
        #   if ($1=="chrVII") print "ref|NC_001139|","\t",$0;
        #   if ($1=="chrVIII") print "ref|NC_001140|","\t",$0;
        #   if ($1=="chrIX") print "ref|NC_001141|","\t",$0;
        #   if ($1=="chrX") print "ref|NC_001142|","\t",$0;
        #   if ($1=="chrXI") print "ref|NC_001143|","\t",$0;
        #   if ($1=="chrXII") print "ref|NC_001144|","\t",$0;
        #   if ($1=="chrXIII") print "ref|NC_001145|","\t",$0;
        #   if ($1=="chrXIV") print "ref|NC_001146|","\t",$0;
        #   if ($1=="chrXV") print "ref|NC_001147|","\t",$0;
        #   if ($1=="chrXVI") print "ref|NC_001148|","\t",$0;
        # }' YPS138.subtelomere.txt > YPS138.subtelomere_Chrs.txt
        #
        # awk -F'\t' '{print $1,"\t",$5,"\t",$6}' YPS138.subtelomere_Chrs.txt > YPS138.subtelomere_Chrs_short.txt

        # awk '
        #         FNR == NR {min[FNR]=$2; max[FNR]=$3; next}
        #         {
        #             for (key in min)
        #                 if (min[key] <= $3 && $3 <= max[key])
        #                     next
        #             print
        #         }' YPS138.subtelomere_Chrs_short.txt mutations_Chrs_TEST.txt > mutations_Clean.txt

# typeset sed_cmd=$(sed '1d' /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Ref_Genome/YPS138_seq_to_remove.txt | awk '{print $2","$3"d"}')
# sed "$cmd" /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D20/mutations.txt > /Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/Muver/D20/mutations_Test.txt
#
#
# # using python, test script:
#
# if mutations['Chromosome'] == remove['Chromosome']
#   & mutations[(mutations['Position'] >= remove['Start']) OR mutations['Position'] <= remove['Stop'])]
#
# import pandas as pd
# mutations = pd.read_table("mutations.txt")
# remove = pd.read_table("YPS138_seq_to_remove_header.txt")
#
# mut = open(input('mutations.txt'),'r').readlines()
# mut.readlines()
# rem = open("YPS138_seq_to_remove_header.txt",'rt')
# rem.readlines()
#
# def (test):
#   if mut['Position'] == rem['Position'] and mut['Position'] <= rem['Start'] or mut['Position'] >= rem['Stop']:
#     pass
#   else:
#     print("keep")
