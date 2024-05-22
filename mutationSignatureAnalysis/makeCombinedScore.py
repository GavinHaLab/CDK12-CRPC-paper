#!/usr/bin/env python3

import os, argparse, pandas as pd

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Input Assignment_Solution file from SigProfiler')
    parser.add_argument('-o', '--output', required=True, help='Output file name')

    args = parser.parse_args()
    return args

args = getArgs()
inputfile = args.input
outputfile = args.output

df = pd.read_csv(inputfile, sep="\t", header=0, index_col=0)
df['combinedScore'] = df['SBS3'] + df['SBS8']
df.drop(['SBS3', 'SBS8'], axis=1, inplace=True)
# now for each row, divide all values by the sum of all values
df = df.div(df.sum(axis=1), axis=0)

# now, subset to only combinedScore column
df = df[['combinedScore']]
# print number of non-zero rows
print(df[df['combinedScore'] > 0].shape[0])
df.to_csv(outputfile, sep="\t", index=True)
