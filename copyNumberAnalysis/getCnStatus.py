#!/usr/bin/env python3

import pandas as pd, os, sys
import argparse

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i',dest="infile",type=str,required=True,help='Input CN matrix file')
    parser.add_argument('-g', '--genes_of_interest',dest="genes",type=str,required=True,help='File with genes of interest')
    parser.add_argument('-p', '--ploidy_file', required=True, help='File with information regaring purity/ploidy per sample (same as used for makeMatrix)')
    parser.add_argument('-gf', '--gene_file', required=True, help='File with gene names and chromosomal locations')
    parser.add_argument('-o',dest="outfile",type=str,required=True,help='Output file')
    
    arg = parser.parse_args()
    
    return arg

args = getArgs()
i = args.infile
o = args.outfile
g = args.genes
p = args.ploidy_file
gf = args.gene_file

sex_gene = pd.read_csv(gf, sep="\t", header=0)
sex_gene = sex_gene[sex_gene['Chr'].isin(['chrX', 'chrY'])]
sex_gene = sex_gene['Gene'].tolist()
sex_gene.append('AREnhancer')

## Read in genes of interest, assign ploidy and purity, add to dataframe with samples as rows and ploidy + genes as columns

## need to change this file now to use median CN across genes instead of ploidy from titan file. 

genes = ['ploidy']
with open(g, 'r') as f:
    for line in f:
        genes.append(line.strip())

ploidy_dict = dict()
samples_orig = {}
with open(p, 'r') as ploidy:
    first = True
    for line in ploidy:
        if first:
            first = False
            continue
        line = line.strip().split('\t')
        sample_old = line[0]
        sample = line[0].split('_cluster')[0]
        samples_orig[sample] = sample_old
        ploidy = float(line[2])
        ploidy_dict[sample] = ploidy

samples = ploidy_dict.keys()
df_out = pd.DataFrame(index=samples, columns=genes)

# for sample in samples:
#     df_out.loc[sample, 'ploidy'] = ploidy_dict[sample]

# sanity check , print out first bit of df out
# print(df_out.head())

## Read in CN matrix, with header and index
cn_mat = pd.read_csv(i, sep='\t', header=0, index_col=0)
# print(cn_mat.index) 
# take median of all gene CNs, set as ploidy
cn_mat['ploidy'] = cn_mat.median(axis=1)
# print(cn_mat['ploidy'])

first_gene = True

print(cn_mat.head())
for gene in genes:
    print(gene)
    # if first_gene:
    #     first_gene = False
    #     continue
    if gene in cn_mat.columns:
        for sample in samples:
            if gene in sex_gene:
                old_name = samples_orig[sample]
                print(f'previous CN {cn_mat.loc[old_name, gene]} for sample {sample}')
                df_out.loc[sample, gene] = cn_mat.loc[old_name, gene] * 2 # multiply 
                print(f'new CN {df_out.loc[sample, gene]}')
                if gene == 'AR' and (sample.startswith('05-144') or sample.startswith('12-011')):
                    # print('here is the result of the multiplication', cn_mat.loc[old_name, gene] * 2)
                    print('working for sex-linked genes', gene, 'with CN', cn_mat.loc[old_name, gene], 'in sample ', sample, 'with ploidy', ploidy_dict[sample])
                    print('here is the result of the multiplication', cn_mat.loc[old_name, gene] * 2)
                continue
            else:
                old_name = samples_orig[sample]
                df_out.loc[sample, gene] = cn_mat.loc[old_name, gene]
                continue
    else:
        print(f'{gene} not in CN matrix')
        continue
        # print(sample, gene)
        # want to see why this is not working
        # old_name = samples_orig[sample]
        # df_out.loc[sample, gene] = cn_mat.loc[old_name, gene]

for sample in samples:
    df_out.loc[sample, 'ploidy'] = cn_mat.loc[samples_orig[sample], 'ploidy']

print(df_out.head())

# divide all rows by sample ploidy values
for sample in samples:
    df_out.loc[sample, genes] = round(df_out.loc[sample, genes].astype(float) / df_out.loc[sample, 'ploidy'], 3)

# round to 2 decimal places 

df_out.to_csv(o, sep='\t', index=True, header=True)