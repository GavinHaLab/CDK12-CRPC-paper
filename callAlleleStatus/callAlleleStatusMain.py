#!/usr/bin/env python3 

__author__ = "Thomas Persse"

import os, argparse, pandas as pd
import re, numpy as np



##############################################
### UPDATED VERSION OF CALLALLELESTATUS.PY ###
###     LAST CHANGES MADE 2023.12.18       ###
##############################################


##############################################
### RUNNING TO-DO LIST FOR THIS SCRIPT:    ###
### 1. Add SV portion to script    DONE    ###
### 2. Add p. or c. to splice variants ADJ ###




def getArgs():
    parser=argparse.ArgumentParser()

    parser.add_argument('-v', '--vcf_dir', required=True, help='Directory containing VCF files')
    parser.add_argument('-c', '--cn_file', required=True, help='File containing CN calls for genes of interest')
    parser.add_argument('-i', '--indel_dir', required=True, help='Directory containing indel files')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-loh', '--loh_file', required=True, help='File containing LOH calls for genes of interest')
    parser.add_argument('-p', '--patients', required=False, help='File containing patient IDs to be included in analysis')
    parser.add_argument('-sv', '--sv_dir', required=False, help='Directory containing SV files')
    parser.add_argument('-ccg', '--cancerCensusGenes', required=False, help='File containing cosmic cancer census genes, to be annotated in output file')
    parser.add_argument('-gf', '--gene_file', required=False, help='File containing genes and coords chr included, currently written to accept ensembl mane format')
    parser.add_argument('-w', '--wordy', required=False, help='Print out more information', action='store_true')
    args=parser.parse_args()
    return args

args = getArgs()
vcf_dir = args.vcf_dir
cn_file = args.cn_file
# matrix_dir = args.matrix_dir
indel_dir = args.indel_dir
output_dir = args.output_dir
loh_file = args.loh_file
if args.patients:
    patients = args.patients
# titan_dir = args.titan_dir
if args.sv_dir:
    sv_dir = args.sv_dir
if args.cancerCensusGenes:
    cancerCensusGenes = args.cancerCensusGenes

if args.gene_file:
    sex_gene = pd.read_csv(args.gene_file, sep="\t", header=0)
    sex_gene = sex_gene[sex_gene['chrom'].isin(['chrX', 'chrY'])]
    sex_gene = sex_gene['name'].tolist()


def clean_index(column):
    # remove _cluster* from index column
    return column.split('_cluster')[0]

if args.wordy:
    print(args.wordy)
    wordy = True
else:
    wordy = False

def find_protein_change_indel(line):
    prots = re.findall(r'p\.\w+\d+\w*(?:\*\d+)?', line) # does this work for snvs?
    if len(prots) == 0:
        return ''
    if len(prots) >= 1: # yet again, changed from > 1 to >= 1... hoping this can rescue some functionality. 
        # print('this is the protein found: ' + prots)
        return prots[0]

def find_protein_change_p53(line):
    if 'NM_000546' in line:
        line = line.split('NM_000546')[1].split(',')[0]
        prot = re.findall(r'p\.\w+\d+\w+', line)
        # print(line)
        if len(prot) == 0:
            return ''
        if len(prot) >= 1:
            return prot[0]
    else:
        prot = re.findall(r'p\.\w+\d+\w+', line)
        if len(prot) == 0:
            return ''
        if len(prot) >= 1:
            # return prot[0]
            return prot[0]

def addCounts(df):
    # using Mutect_Tumor, Varscan_TUmor, and Strelka_Tumor columns, filter out rows where both have < 5 Alt reads and less than 10 total reads
    # first, need to parse out the alt and ref reads from the columns
    # mutect format = GT:AD:AF:DP:F1R2:F2R1:SB, need to parse out AD and DP
    # varscan format = GT:GQ:DP:RD:AD:FREQ:DP4, need to parse out AD and DP
    # strelka format = GT:DP:FDP:SDP:SUBDP:AU:CU:GU:TU, need to parse out AU and TU
    # first, need to check if the columns are empty, if so, return empty dataframe
    tumor_columns = ['Mutect_Tumor', 'Varscan_Tumor', 'Strelka_Tumor']
    df[tumor_columns] = df[tumor_columns].fillna(0)
    
    if df['Mutect_Tumor'].notna().all(): # if all values are not NaN
        df['Mutect_AD'] = df['Mutect_Tumor'].apply(lambda x: x.split(':')[1].split(',')[1] if x != 0 else 0)
        df['Mutect_DP'] = df['Mutect_Tumor'].apply(lambda x: x.split(':')[3] if x != 0 else 0)
    else:
        df['Mutect_AD'] = 0
        df['Mutect_DP'] = 0

    if df['Varscan_Tumor'].notna().all(): #
        df['Varscan_AD'] = df['Varscan_Tumor'].apply(lambda x: x.split(':')[4] if x != 0 else 0)
        df['Varscan_DP'] = df['Varscan_Tumor'].apply(lambda x: x.split(':')[2] if x != 0 else 0)
    else:
        df['Varscan_AD'] = 0
        df['Varscan_DP'] = 0

    if df['Strelka_Tumor'].notna().all():
        df['Strelka_AD'] = df['Strelka_Tumor'].apply(lambda x: x.split(':')[5].split(',')[1] if x != 0 else 0)
        df['Strelka_DP'] = df['Strelka_Tumor'].apply(lambda x: x.split(':')[1] if x != 0 else 0)
    else:
        df['Strelka_AD'] = 0
        df['Strelka_DP'] = 0


    # convert to numeric
    df['Mutect_AD'] = pd.to_numeric(df['Mutect_AD'])
    df['Mutect_DP'] = pd.to_numeric(df['Mutect_DP'])
    df['Varscan_AD'] = pd.to_numeric(df['Varscan_AD'])
    df['Varscan_DP'] = pd.to_numeric(df['Varscan_DP'])
    df['Strelka_AD'] = pd.to_numeric(df['Strelka_AD'])
    df['Strelka_DP'] = pd.to_numeric(df['Strelka_DP'])

    df['Max_AD'] = df[['Mutect_AD', 'Varscan_AD', 'Strelka_AD']].max(axis=1)
    df['Max_DP'] = df[['Mutect_DP', 'Varscan_DP', 'Strelka_DP']].max(axis=1)
    # make all gnomAD NaNs 0
    df['gnomAD_genome_ALL'] = df['gnomAD_genome_ALL'].fillna(0)
    df = df[df['gnomAD_genome_ALL'] < 0.1]
    return df
    # for sanity, print out mean, median, and quantiles of Max_AD and Max_DP
    # print(df['Max_AD'].mean())
    # print(df['Max_AD'].median())
    # print(df['Max_AD'].quantile([0.25, 0.75]))
    # print(df['Max_DP'].mean())
    # print(df['Max_DP'].median())
    # print(df['Max_DP'].quantile([0.25, 0.75]))
    # return df

def addCountsIndel(df):
    ### need to add this function to the indel portion of the script, filtering out rows where AD < 5 and DP < 10
    # using svaba calls for now, so input format is GT:AD:DP:GQ:PL:SR:CR:LR:LO
    # we need AD and DP, will add as column
    df['AD'] = df['Otherinfo14'].str.split(':').str[1]
    df['DP'] = df['Otherinfo14'].str.split(':').str[2]
    df['AD'] = pd.to_numeric(df['AD'])
    df['DP'] = pd.to_numeric(df['DP'])

    # for sanity, print out mean, median, and quantiles of AD and DP
    # print(df['AD'].mean())
    # print(df['AD'].median())
    # print(df['AD'].quantile([0.25, 0.75]))
    # print(df['DP'].mean())
    # print(df['DP'].median())
    # print(df['DP'].quantile([0.25, 0.75]))

    return df


def formatIndex(line):
    # will apply to index column, split on _ and return first two elements joined by : 
    return ':'.join(line.split('_')[:2])

def callLOH(df_col):
    # look across column, report 'LOH' if value is 1, empty string otherwise
    if df_col == 1:
        return 'LOH'
    else:
        return ''

def callCNV(df_col):
    # look across column, annotate with gain, amp, loss, del. leave out neutral, since we only want evidence
    if df_col >= 2.5:
        return 'AMP'
    elif df_col >= 1.5 and df_col < 2.5:
        return 'GAIN'
    elif df_col <= 0.5 and df_col > 0:
        return 'DEL'
    elif df_col == 0:
        return 'HomDel'
    else:
        return ''

def getPatientSamples(pat_file, samples_list):
    ## patient ID is subset of sample ID, so we can use this to subset samples
    pat_samps = []
    with open(patients, 'r') as f:
        for line in f:
            line = line.strip()
            for sample in samples_list:
                if line in sample:
                    pat_samps.append(sample)
    return pat_samps


##  get indel filepaths:
indel_locs = {}

# for file in os.listdir(indel_dir):
#     # print(dir)
#     # for file in os.listdir(d):
#     # print('yes')
#     name = file
#     if ".DS_Store" in name:
#         continue
#     file = file + '/' + file + '.svaba.somatic.indel.hg38_multianno.txt'
#     # print(file)
#     # if file.endswith('svaba.somatic.indel.hg38_multianno.txt'):
#         # print('yes')
#         # print(dir)
#     sample = file.split('.')[0]
#     indel_locs[name] = indel_dir + file
    # print(indel_locs[name])

## read in CN matrix, which is also where we get our subset of genes
cn = pd.read_csv(cn_file, sep="\t", header=0, index_col=0)
samples = cn.index
sample_list = list(samples)
genes = cn.columns
gene_list = list(genes)
if 'ploidy' in genes:
    gene_list = gene_list[1:]

# if patients file is provided, subset samples to patients of interest
if args.patients:
    samples = getPatientSamples(patients, samples)
    sample_list = list(samples)
    # subset cn matrix to samples of interest
    cn = cn.loc[samples, :]
    # subset indel_locs to samples of interest
    indel_locs = {k: v for k, v in indel_locs.items() if k in samples}


# create empty dictionary to store evidence per sample, per gene
evidence = {}
for sample in sample_list:
    evidence[sample] = {}
    for gene in gene_list:
        evidence[sample][gene] = {'CNV': [], 'LOH': []}

# want to parse cancer census genes, if provided, make two lists: TSGs and OGs

# ic = ic[ic['Role in Cancer'].str.contains('oncogene|fusion', case=False)]
# print(ic['Role in Cancer'].unique())


# XXX: work on this
if args.cancerCensusGenes:
    tsgs = []
    ogs = []
    ic = pd.read_csv(args.cancerCensusGenes, sep="\t", header=0, index_col=0)
    print(ic)
    ic = ic.fillna('NA')
    oncogenes = ic[ic['Role in Cancer'].str.contains('oncogene', case=False)]
    tsgs = ic[ic['Role in Cancer'].str.contains('TSG', case=False)]
    tsgs = tsgs.index.tolist()
    oncogenes = oncogenes.index.tolist()
    print(len(oncogenes))
    print(len(tsgs))
    onco = set(oncogenes)
    tsg = set(tsgs)
    # determine difference between oncogenes and tsgs
    print(onco.difference(tsg))
    print(len(onco.difference(tsg)))
    print(tsg.difference(onco))
    print(len(tsg.difference(onco)))
    if 'TP53' in oncogenes and 'TP53' in tsgs:
        # remove TP53 from oncogenes
        oncogenes.remove('TP53')
    print('not differentiating tsg and oncogenes for now')
    # currently omitting oncogene/tsg differentiation, due to mixed reports for 
for gene in gene_list:
    if gene in oncogenes and gene in tsgs:
        print('Gene ' + gene + ' is in both oncogenes and TSGs')
    elif gene in oncogenes and gene not in tsgs:
        print('Gene ' + gene + ' is in oncogenes')
    elif gene in tsgs and gene not in oncogenes:
        print('Gene ' + gene + ' is in TSGs')
    else:
        print('Gene ' + gene + ' is in neither oncogenes nor TSGs')


##### LOH PORTION #####

print('\n Calling LOH \n')

# read in loh matrix, subset to genes of interest
loh = pd.read_csv(loh_file, sep="\t", header=0, index_col=0)
loh.index = loh.index.str.split('_cluster').str[0]
loh = loh[gene_list]
loh = loh.applymap(callLOH) # applymap applies function to each element in dataframe

# add loh to evidence dictionary
print(gene_list)
        

for sample in sample_list:
    for gene in gene_list:
        # print('final one', type(loh))
        loh_val = loh.loc[sample, gene]
        if loh_val == '':
            continue
        else:
            evidence[sample][gene]['LOH']=(loh.loc[sample, gene])

print('\n Calling CNVs \n')

##### CN PORTION #####

# read in cn matrix, annotate with cn_gain, cn_loss, cn_neutral, cn_amp, cn_del
cn = pd.read_csv(cn_file, sep="\t", header=0, index_col=0)
# already subset to genes of interest.
# apply callCNV function to each column
cn = cn.applymap(callCNV)

# add cn to evidence dictionary
for sample in sample_list:
    for gene in gene_list:
        if cn.loc[sample, gene] == '':
            continue
        else:
            evidence[sample][gene]['CNV']=(cn.loc[sample, gene])


##### COMBINE CNV, LOH HERE #####

# adju for gain, amp, loss, del, homdel, loh

for sample in evidence:
    for gene in evidence[sample]:
        cnv = evidence[sample][gene]['CNV']
        loh = evidence[sample][gene]['LOH']
        if cnv == [] and loh == []: # if both are empty, set to empty list
            evidence[sample][gene] = []
        elif cnv == [] and loh != []: # if cnv is empty, but loh is not, set to loh
            evidence[sample][gene] = ['LOH']
        elif cnv != [] and loh == []: # if loh is empty, but cnv is not, set to cnv
            evidence[sample][gene] = [cnv]
        elif cnv != [] and loh != []: # if both are not empty, set to both
            print(cnv, loh)
            if cnv == 'HomDel' and loh == 'LOH':
                evidence[sample][gene] = ['HomDel']
            else:
                evidence[sample][gene] = [f'{cnv}_LOH']
            # if cnv == 'DEL' and loh == 'LOH': # if both are del, set to homdel
            #     evidence[sample][gene] = ['DEL_LOH']
            # elif cnv == 'AMP' and loh == 'LOH':
            #     evidence[sample][gene] = [cnv, loh]
        else:
            print('something went wrong')
            # print(evidence[sample][gene])

##### SNV/INDEL PORTION #####

print('\n Calling Indels/SNVs \n')

# read in snv files, subset to genes of interest
# this portion isn't working as needed... need to figure out how to subset to genes of interest
snv_paths = [os.path.join(vcf_dir, f) for f in os.listdir(vcf_dir) if f.endswith('twoormore.txt')]
# print(snv_paths)
for snv_path in snv_paths:
    snv = pd.read_csv(snv_path, sep="\t", header=0, index_col=0)
    # fill spaces with underscores
    snv = snv.replace(' ', '_', regex=True)
    sample_name = snv_path.split('/')[-1].split('_twoormore')[0]
    if sample_name.startswith('~$'):
        continue

    if sample_name == 'FH0826_E_C1':
        print('\n\n\nTP53 here')
        print(snv)
        print(snv[snv['Gene.refGene'] == 'TP53'])
        print('\n\n\n')
    try:
        snv = addCounts(snv)
    except:
        print(f'broken for this sample: {sample_name}, unsure why')
        print(snv.columns)
        exit()
    # now, filter out rows where Max_AD < 5 and Max_DP < 10
    if sample_name == 'FH0826_E_C1':
        print('\n\n\nTP53 here')
        print(snv)
        # print(snv['Gene.refGene']== 'TP53')
        # print only rows where gene is TP53
        print(snv[snv['Gene.refGene'] == 'TP53'])
        print('\n\n\n')
    snv = snv[(snv['Max_AD'] >= 5) & (snv['Max_DP'] >= 10)]

    # then, remove columns other than gene, Func.refGene, ExonicFunc.refGene, AAChange.refGene, CLNSIG, Max_AD, Max_DP
    snv = snv[['Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'CLNSIG', 'Max_AD', 'Max_DP']]
    # snv = snv[['Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'CLNSIG']]
    # print(sample_name)
    # subset to rows with gene col values in gene_list
    snv = snv[snv['Gene.refGene'].isin(gene_list)] # subset to calls that are not synonymous
    snv = snv[snv['ExonicFunc.refGene'] != 'synonymous_SNV']
    snv = snv[snv['CLNSIG'] != 'Benign']
    snv = snv[snv['CLNSIG'] != 'Likely_benign']
    if sample_name == 'FH0826_E_C1':
        print('\n\n\nTP53 here')
        print(snv)
        print(snv[snv['Gene.refGene'] == 'TP53'])
        print('\n\n\n')
    # wnat to see what the columns are here
    # 
    # add snv to evidence dictionary
    for gene in gene_list:
        subset = snv[snv['Gene.refGene'] == gene]
        # print(subset['AAChange.refGene'].tolist())
        exonic = subset[subset['Func.refGene'] == 'exonic']
        splicing = subset[subset['Func.refGene'] == 'splicing']

        try: # have had some issues with missing columns, so this is a temporary fix
            if gene == 'TP53':
                changes = exonic['AAChange.refGene'].apply(find_protein_change_p53)
            else:
                changes = exonic['AAChange.refGene'].apply(find_protein_change_indel)
        except:
            print('No AAChange.refGene column found for sample ' + sample_name + ', and gene ' + gene)
            continue

        changes = changes[changes != ''] # remove empty strings
        # add varID to changes while applying formatIndex function
        if wordy:
            changes = changes + ';' + exonic.index.map(formatIndex)
        else:
            changes = changes
        # print(splicing)
        # splicing = splicing + ';' + splicing.index.map(formatIndex)
        splicing = splicing['Func.refGene'] + ';' + splicing.index.map(formatIndex)
        changes = changes.tolist()
        if len(splicing) > 0:
            # print(splicing['AAChange.refGene'])
            # print(len(splicing))
            for splice in splicing.tolist():
                changes.append(splice)
            # changes.append('splice')
        if changes == []:
            print('no snvs found for sample ' + sample_name + ' and gene ' + gene)
            continue
        # changes = ":".join(changes)
        for change in changes:
            if sample_name in evidence:
                evidence[sample_name][gene].append(change)
                print(evidence[sample_name][gene])
            else:
                print('Sample ' + sample_name + ' not found in evidence dictionary')

print('\n Calling Indels \n')


### INDEL PORTION ###
### SWAPPING FOR NEW METHOD USING COMBINED CALLS, SO WILL REWRITE THIS SECTION ###
# read in indel files, subset to genes of interest
# for sample in indel_locs:
#     ipath = indel_locs[sample]
#     # print(ipath)
#     try:
#         indel = pd.read_csv(ipath, sep="\t", header=0, index_col=0)
#     except:
#         print('No indel file found for sample ' + sample)
#         continue
#     sample_name = sample
#     # subset to rows with gene col values in gene_list
#     indel = indel[indel['Gene.refGene'].isin(gene_list)]
#     indel = addCountsIndel(indel)

#     # now, filter out rows where AD < 5 and DP < 10
#     indel = indel[(indel['AD'] >= 5) & (indel['DP'] >= 10)]
#     # ensure no benign calls go through
#     indel = indel[indel['CLNSIG'] != 'Benign']
#     # loop through gene list in subset df, apply find_protein_change_indel function to AAChange.refGene column, add to evidence dictionary
#     for gene in gene_list:
#         subset = indel[indel['Gene.refGene'] == gene]
#         # subset = subset[(subset['Func.refGene'] == 'exonic') | (subset['Func.refGene'] == 'splicing')]
#         exonic = subset[subset['Func.refGene'] == 'exonic']
#         splicing = subset[subset['Func.refGene'] == 'splicing']
#         # print(len(splicing))
#         if gene == 'TP53':
#             changes = exonic['AAChange.refGene'].apply(find_protein_change_p53)
#         else:
#             changes = exonic['AAChange.refGene'].apply(find_protein_change_indel)
#         # parse splicing column, add 'splice' to changes list
#         changes = changes[changes != '']
#         if wordy:
#             changes = changes + ';' + exonic.index.map(formatIndex)
#         else:
#             changes = changes

#         splicing = splicing['Func.refGene'] + ';' + splicing.index.map(formatIndex)
#         changes = changes.tolist()
#         if len(splicing) > 0:
#             # print(splicing['AAChange.refGene'])
#             # print(len(splicing))
#             for splice in splicing.tolist():
#                 changes.append(splice)
#         if changes == []:
#             continue
#         for change in changes:
#             if sample_name in evidence:
#                 evidence[sample_name][gene].append(change)
#                 print(change)
#             else:
#                 print('Sample ' + sample_name + ' not found in evidence dictionary')

### END OF OLD INDEL PORTION ###
# indel_dir = [os.path.join(indel_dir, f) for f in os.listdir(indel_dir) if f.endswith('union_indels.txt')]
# print(indel_dir)
for file in os.listdir(indel_dir):
    print(file)
    if file.endswith('.txt'):
        sample_name= file.split('/')[-1].split('.svaba')[0]
        print('Currently analyzing ' + sample_name)
        file = os.path.join(indel_dir, file)
        # now, read indel file, subset to genes of interest, and add to evidence dictionary
        # note that here, we have already parsed values, so for aachange we use p.value column, check if splicing using Func.refGene column and adding chr:start to the call
        indel = pd.read_csv(file, sep="\t", header=0, index_col=0)
        indel = indel[indel['Gene.refGene'].isin(gene_list)]
        indel = indel[indel['CLNSIG'] != 'Benign']
        indel = indel[indel['CLNSIG'] != 'Likely_benign']
        # indel = indel[indel['CLNSIG'] != 'NA']

        for gene in gene_list:
            subset = indel[indel['Gene.refGene'] == gene]
            exonic = subset[subset['Func.refGene'] == 'exonic']
            splicing = subset[subset['Func.refGene'] == 'splicing']
            if gene == 'TP53':
                changes = exonic['AAChange.refGene'].apply(find_protein_change_p53)
            else:
                changes = exonic['AAChange.refGene'].apply(find_protein_change_indel)
            # splicing = splicing['Func.refGene'] + ';' + splicing['Chr'] + ':' + splicing['Start']
            # Chr is in index, so need to use index instead of calling column when reporting splicing site
            splice = splicing['Func.refGene'] + ';' + splicing.index + ':' + splicing['Func.refGene']
            changes = changes.tolist()
            if len(splicing) > 0:
                for splice in splicing.tolist():
                    changes.append(splice)
            if changes == []:
                continue
            for change in changes:
                if sample_name in evidence:
                    evidence[sample_name][gene].append(change)
                else:
                    print('Sample ' + sample_name + ' not found in evidence dictionary')




if args.sv_dir:
    print('\n Calling SVs \n')
    # sv dir should be structured as svDir/sample.svabaTitan.genes.sv.txt
    sv_files = [os.path.join(sv_dir, f) for f in os.listdir(sv_dir) if f.endswith('genes.sv.txt')]
    for sv_file in sv_files:
        # will have columns 'gene1' and 'gene2', so we need to loop through each gene of interest, subset to it, and add to evidence dictionary
        # will probably want to make id column from the following columns: gene1:gene2; type, (chromosome_1:start_1,chromosome_2:start_2)

        sample = sv_file.split('/')[-1].split('.')[0]
        print('Currently analyzing ' + sample)
        sv = pd.read_csv(sv_file, sep="\t", header=0, index_col=0)
        for gene in gene_list:
            # subset = sv[sv['gene1'] == gene]
            # subset2 = sv[sv['gene2'] == gene]
            # subset3
            # need to do subsetting so that gene1 and gene2 are not the same gene for 1 and 2, then create new subset3 that is for gene1 and gene2 being the same gene
            subset = sv[(sv['gene1'] == gene) & (sv['gene2'] != gene)]
            # replace gene2 NAs with 'Intergenic'
            subset['gene2'] = subset['gene2'].fillna('Intergenic')
            subset2 = sv[(sv['gene2'] == gene) & (sv['gene1'] != gene)]
            subset2['gene1'] = subset2['gene1'].fillna('Intergenic')
            subset3 = sv[(sv['gene1'] == gene) & (sv['gene2'] == gene)]
            # want to replace nas in gene1 and gene2 with 'None'
            # subset = subset.append(subset2) can't append, need to concatenate
            # subset = pd.concat([subset, subset2])
            subset = pd.concat([subset, subset2, subset3])
            subset = subset.fillna('None')
            susbset = subset.drop_duplicates()
            # print(subset)
            if subset.empty:
                continue
            else:
                if wordy:
                    subset['id'] = subset['gene1'] + ':' + subset['gene2'] + ' ' + subset['type'] + ',(' + subset['chromosome_1'] + ':' + subset['start_1'].astype(str) + ',' + subset['chromosome_2'] + ':' + subset['start_2'].astype(str) + ')' # this is inducing nan values, need to fix
                else:
                    subset['id'] = subset['gene1'] + ':' + subset['gene2'] + ' ' + subset['type']
                # debug introduced nan values, so print subset where id = nan
                
                # print(subset[subset['id'].isna()])
                # print(subset['id'])
                if sample in evidence:
                    evidence[sample][gene].extend(subset['id'].tolist())
                    # print(evidence[sample][gene])
                else:
                    print('Sample ' + sample + ' not found in evidence dictionary')


# create new dataframe with samples as rows and genes as columns

# first, create new columns (2 per gene), with status and 0 or 1 if gene is mutated. Colnames will be geneName and geneName_altered
## 2024.01.17: change to have three columns per gene: status (BAL/MAL/Intact), allele A (mutation 1), and allele B (mutation 2)
# gene_list_new = [gene + '_altered' for gene in gene_list]
gene_list_alleleA = [gene + '_alleleA' for gene in gene_list]
gene_list_alleleB = [gene + '_alleleB' for gene in gene_list if gene not in sex_gene]
gene_status = [gene + '_status' for gene in gene_list]
gene_list_new = gene_status + gene_list_alleleA + gene_list_alleleB
# gene_list = gene_list_new + gene_list
# order by gene
gene_list_new.sort()

print('writing to file, please wait...')

final = pd.DataFrame(columns=gene_list_new, index=samples)
# populate dataframe with evidence
# need to add section that calls homdels as BAL, not MAL
for sample in evidence:
    for gene in evidence[sample]:
        # this is where to add TSG/OG differentiation, alternate reporting...
        # print(evidence[sample][gene])
        if len(evidence[sample][gene]) == 0:
            final.loc[sample, gene + '_status'] = 'Intact'
            final.loc[sample, gene + '_alleleA'] = 'Intact'
            if gene not in sex_gene:
                final.loc[sample, gene + '_alleleB'] = 'Intact'
            # do we need to list alleles as none/NA here? check emails
            # final.loc[sample, gene + '_altered'] = 0
        else:
            try:
                if gene not in oncogenes:
                    if len(evidence[sample][gene]) == 1:
                        if 'HomDel' in evidence[sample][gene] or 'HomDel_LOH' in evidence[sample][gene]:
                            final.loc[sample, gene + '_status'] = 'BAL'
                            final.loc[sample, gene + '_alleleA'] = 'HomDel'
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleB'] = 'HomDel'
                        elif evidence[sample][gene][0] == 'GAIN' or evidence[sample][gene][0] == 'AMP':
                            final.loc[sample, gene + '_status'] = evidence[sample][gene][0]
                            final.loc[sample, gene + '_alleleA'] = evidence[sample][gene][0]
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleB'] = 'Intact'
                            
                        else:
                            final.loc[sample, gene + '_status'] = 'MAL'
                            final.loc[sample, gene + '_alleleA'] = evidence[sample][gene][0]
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleB'] = 'Intact'
                    else:
                        if 'HomDel' in evidence[sample][gene] or 'HomDel_LOH' in evidence[sample][gene]:
                            # remove loh from list
                            final.loc[sample, gene + '_status'] = 'BAL'
                            final.loc[sample, gene + '_alleleA'] = 'HomDel'
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleB'] = 'HomDel'
                            # evidence[sample][gene] = [x for x in evidence[sample][gene] if x != 'LOH']
                            # final.loc[sample, gene] = f'BAL; {",".join(evidence[sample][gene])}'
                        elif 'GAIN' in evidence[sample][gene][0] or 'AMP' in evidence[sample][gene][0]:
                            final.loc[sample, gene + '_status'] = evidence[sample][gene][0]
                            final.loc[sample, gene + '_alleleA'] = evidence[sample][gene][0]
                            if len(evidence[sample][gene]) == 2:
                                if gene not in sex_gene:
                                    final.loc[sample, gene + '_alleleB'] = 'Intact'
                            else:
                                if gene not in sex_gene:
                                    final.loc[sample, gene + '_alleleB'] = ",".join(evidence[sample][gene][1:])
                                else:
                                    final.loc[sample, gene + '_alleleA'] = ",".join(evidence[sample][gene][0:])

                        else:
                            final.loc[sample, gene + '_status'] = 'BAL'
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleA'] = evidence[sample][gene][0]
                                final.loc[sample, gene + '_alleleB'] = ",".join(evidence[sample][gene][1:])
                            else:
                                final.loc[sample, gene + '_alleleA'] = ",".join(evidence[sample][gene][0:])
                else:
                    # need to check for gain/amp, but not set status to MAL/BAL
                    if len(evidence[sample][gene]) == 1:
                        if 'HomDel' in evidence[sample][gene] or 'HomDel_LOH' in evidence[sample][gene]:
                            final.loc[sample, gene + '_status'] = 'NA'
                            final.loc[sample, gene + '_alleleA'] = 'HomDel'
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleB'] = 'HomDel'
                            print('the unlikely happened')
                            print(evidence[sample][gene], gene, sample)
                        elif evidence[sample][gene][0] == 'GAIN' or evidence[sample][gene][0] == 'AMP':
                            final.loc[sample, gene + '_status'] = evidence[sample][gene][0]
                            final.loc[sample, gene + '_alleleA'] = evidence[sample][gene][0]
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleB'] = 'Intact'
                        else:
                            final.loc[sample, gene + '_status'] = 'Intact'
                            final.loc[sample, gene + '_alleleA'] = evidence[sample][gene][0]
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleB'] = 'Intact'
                    else:
                        if 'HomDel' in evidence[sample][gene] or 'HomDel_LOH' in evidence[sample][gene]:
                            final.loc[sample, gene + '_status'] = 'NA'
                            final.loc[sample, gene + '_alleleA'] = 'HomDel'
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleB'] = 'HomDel'
                            print('the unlikely happened')
                        elif 'GAIN' in evidence[sample][gene][0] or 'AMP' in evidence[sample][gene][0]:
                            final.loc[sample, gene + '_status'] = evidence[sample][gene][0]
                            final.loc[sample, gene + '_alleleA'] = evidence[sample][gene][0]
                            if len(evidence[sample][gene]) == 2:
                                if gene not in sex_gene:
                                    final.loc[sample, gene + '_alleleB'] = evidence[sample][gene][1]
                                else:
                                    final.loc[sample, gene + '_alleleA'] = ",".join(evidence[sample][gene][0:])
                            else:
                                if gene not in sex_gene:
                                    final.loc[sample, gene + '_alleleB'] = ",".join(evidence[sample][gene][1:])
                                else:
                                    final.loc[sample, gene + '_alleleA'] = ",".join(evidence[sample][gene][0:])
                        else:
                            final.loc[sample, gene + '_status'] = 'Intact'
                            final.loc[sample, gene + '_alleleA'] = evidence[sample][gene][0]
                            if gene not in sex_gene:
                                final.loc[sample, gene + '_alleleB'] = ",".join(evidence[sample][gene][1:])
                            else:
                                final.loc[sample, gene + '_alleleA'] = ",".join(evidence[sample][gene][0:])
            except:
                print('Error with sample ' + sample + ' and gene ' + gene)
                print(evidence[sample][gene])
# write dataframe to file
if wordy:
    final.to_csv(output_dir + '/AlleleStatusVerbose.txt', sep="\t", index=True, header=True)
else:
    final.to_csv(output_dir + '/AlleleStatusARfix.txt', sep="\t", index=True, header=True)
# final.to_csv(output_dir + '/haffnerGeneStatus.txt', sep="\t", index=True, header=True)
# python3 /fh/fast/ha_g/projects/ProstateTAN/analysis_WES/CDK12_project/callAlleleStatusMain.py -v /fh/fast/ha_g/projects/ProstateTAN/analysis_WES/SNV_calling/combined_TwoOrMore -c /fh/fast/ha_g/projects/ProstateTAN/analysis_WES/haffnerProject/testCN.txt -i /fh/fast/ha_g/projects/ProstateTAN/analysis_WES/CN-SV_calling/SVABA/ -o ../haffnerProject/geneStatus/ -loh /fh/fast/ha_g/projects/ProstateTAN/analysis_WES/CN-SV_calling/TitanCNA/GeneCnStatus/TAN_WES_geneLOH.txt -p /fh/fast/ha_g/projects/ProstateTAN/analysis_WES/CDK12_project/haffnerPatients.txt