#!/usr/bin/bash

# see script for environment dependencies

/usr/bin/time -v python3 callAlleleStatusMain.py -v /path/to/annovar/snv/dir -c /path/to/ploidy/adjusted/cn/calls.txt -i /path/to/indel/calls/dir/ -o /path/to/output/directory/ -loh /path/to/loh/calls/from/makeMatrix/endswithLOH.txt -gf /path/to/transcript/file.txt -ccg /path/to/list/of/cosmic/census/genes.tsv -sv /path/to/annotated/sv/calls/dir/