#!/usr/bin/env python3

import os, argparse, pandas as pd

def getArgs():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Input directory containing VCF files')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-n', '--names', required=False, help='File containing sample names to change')
    args=parser.parse_args()
    return args

def changenames(name, names_dict):
    if name in names_dict.keys():
        return names_dict[name]
    else:
        return name

args = getArgs()
inputdir = args.input
outputdir = args.output

names_dict = {}
if args.names:
    with open(args.names, 'r') as f:
        first = True
        for line in f:
            if first:
                first = False
                continue
            line = line.strip().split("\t")
            names_dict[line[0]] = line[1]

if not inputdir.endswith("/"):
    inputdir += "/"

if not outputdir.endswith("/"):
    outputdir += "/"

for file in os.listdir(inputdir):
    path = inputdir + file
    sample = file.split(".vcf")[0]
    outfile = outputdir + sample + ".trimmed.vcf"
    inf = pd.read_csv(path, sep="\t", comment="#", header=0)
    out = pd.DataFrame(columns=['chr', 'pos', 'sample', 'ref', 'alt'])
    out['sample'] = [sample] * len(inf)
    out['chr'] = inf['Chr']
    out['pos'] = inf['Start']
    out['ref'] = inf['Ref']
    out['alt'] = inf['Alt']
    out.to_csv(outfile, sep="\t", index=False, header=False)




