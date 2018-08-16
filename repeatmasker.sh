
reference_folder="/lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/cerevisiae"
consensus_te_seqs="/lustre1/hcm14449/mcc/mcc_8_16_18/test/sac_cer_TE_seqs.fasta"
reference_genome="/lustre1/hcm14449/TE_MA_Paradoxus/ref_genome/cerevisiae/sacCer2.fasta"

RepeatMasker -dir $reference_folder -pa $PBS_NP -lib $consensus_te_seqs -s -gff -nolow -no_is $reference_genome
