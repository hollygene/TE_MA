### First want to remove lines in the header that are unnecessary
# anything with "Flag" can be removed

# import pandas module
import pandas as pd

# making data frame from text file

data = pd.read_table("mutations.txt")
# to do this with a regex
# data.filter(like='Flag')
data.drop(data.filter(like='Ploidy'),axis=1,inplace=True)
data.drop(data.filter(like='Flag'),axis=1,inplace=True)
data.drop(data.filter(like='Subclonal'),axis=1,inplace=True)
# to save modifications to a new file:
data.to_csv('test.txt', sep='\t', index=False)

# to get # of mutations (all types)
data = pd.read_table("D20_mutations_Chrs_CleanER.txt")
dat = pd.Series(data.values.flatten())
dat.str.count("g.").sum()
dat.str.count("del").sum()
dat.str.count("ins").sum()
dat.str.count(">").sum()

dataH = pd.read_table("H0_mutations_Chrs_Clean.txt")
datH = pd.Series(dataH.values.flatten())
datH.str.count("g.").sum()
datH.str.count("del").sum()
datH.str.count("ins").sum()
datH.str.count(">").sum()

dataD0 = pd.read_table("D0_mutations_Chrs_CleanER.txt")
datD0 = pd.Series(dataD0.values.flatten())
datD0.str.count("g.").sum()
datD0.str.count("del").sum()
datD0.str.count("ins").sum()
datD0.str.count(">").sum()

dataD1 = pd.read_table("D1_mutations_Chrs_CleanER.txt")
datD1 = pd.Series(dataD1.values.flatten())
datD1.str.count("g.").sum()
datD1.str.count("del").sum()
datD1.str.count("ins").sum()
datD1.str.count(">").sum()
