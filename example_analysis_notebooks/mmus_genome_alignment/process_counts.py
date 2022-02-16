# Import necessary packages
import pandas as pd
import numpy as np
import glob
from functools import reduce

########################## Raw Read Count ###########################

# Define the path that you are going to take the files from.
path = '/home/kflanaga/scratch' # use your path
all_files = glob.glob(path + "/counts*.txt")

# Create an empty list to put each dataframe into.
li = []

# Fill the list with all of the count dataframes.
for filenames in all_files:
    df = pd.read_csv(filenames,
                 sep='\t',
                header=1)
    li.append(df)

# Remove the uneccessary columns and shorten the sample names.
for i in range(len(li)):
    sample = li[i]
    sample = sample.drop(["Chr","Start","End","Strand","Length"], axis=1)
    sample.columns = sample.columns.str.replace("_Aligned.sortedByCoord.out.bam","")
    li[i] = sample

# Merge all of the dataframes together using an inner join.
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['Geneid'],
                                            how='inner'), li)

# Save the merged dataframe as the raw read count input for RiboDiff
df_merged.to_csv("raw_read_count.txt", sep = '\t', index = False)

########################## Experimental Design ###########################

# Creating a list of all of the Samples
Samples = list(df_merged.columns[1:])

#Using a couple of string searches to record if each sample is a control or not.
Conditions = []
for i in Samples:
    if "Fmr1" in i:
        Conditions.append("treated")
    elif "control" in i:
        Conditions.append("control")

# Using a couple of string searches to record if each sample is made
# of RNA sequences or ribosome profiles
Data_Type = []
for i in Samples:
    if "RNA" in i:
        Data_Type.append("RNA-Seq")
    elif "RPF" in i:
        Data_Type.append("Ribo-Seq")

# Combine all of the lists together into a DataFrame
exp_d = pd.DataFrame({"Samples": Samples,
              "Data_Type": Data_Type,
              "Conditions": Conditions})

# Saving the dataframe as the experimental_design csv
exp_d.to_csv("experimental_design.csv", index = False)
