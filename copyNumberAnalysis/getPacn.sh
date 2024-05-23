#!/usr/bin/bash

python3 getCnStatus.py -i /path/to/makeMatrix/outputCN.txt -g genesOfInterest.txt -p purityPloidy.txt -gf /path/to/geneFile.txt -o /path/to/outputFile.txt

# note: gene file ( -gf) used was the mane transcript file from ensembl
