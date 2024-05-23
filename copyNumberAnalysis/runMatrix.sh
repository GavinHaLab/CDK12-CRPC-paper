#!/usr/bin/bash

ml R/3.6.2-foss-2019b-fh1

Rscript makeMatrixFromTITAN-ICHOR_segBased_hg38.R /path/to/titan/run/dir/ /path/to/canonical/transcripts/file.txt /path/to/purity/ploidy/file.txt common Gene 1000 nameOfRun