






for i in `ls -1 /scratch/jc33471/paradoxusHolly/run0217all/out/Spar | grep -v "reference"`;do

grep "TY1_non-ref" /scratch/jc33471/paradoxusHolly/run0217all/out/Spar/${i}/results/*_telocate_nonredundant.bed > /pathToYourResultFolder/${i}_telocate_red.TY1.bed

done

# or alternatively, using raw output (there might be more predictions, but with potential redundant)


for i in `ls -1 /scratch/jc33471/paradoxusHolly/run0217all/out/Spar | grep -v "reference"`;do

grep "TY1_non-ref" /scratch/jc33471/paradoxusHolly/run0217all/out/Spar/${i}/results/*_telocate_raw.bed > /pathToYourResultFolder/${i}_telocate_raw.TY1.bed

done
