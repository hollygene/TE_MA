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
dat = pd.Series(data.values.flatten())
dat.str.count("g.").sum()
dat.str.count("del").sum()
dat.str.count("ins").sum()
dat.str.count(">").sum()
