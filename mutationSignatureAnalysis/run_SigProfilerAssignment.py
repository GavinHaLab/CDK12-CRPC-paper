#!/usr/bin/python
__author__="Pushpa Itagi, PhD"
__purpose__="To generate plots from existing vcfs to fit to the COSMIC SBS/DBS/INDEL signatures. Please select the right functions and output and input directory variables accordingly"

## EDITS MADE BY THOMAS PERSSE

print("loading modules")
import sys
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
import pandas as pd

import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
import sys
import os
print(sys.path)

build = "GRCh38"
from SigProfilerMatrixGenerator import install as genInstall

#Function for the SBS signatures
def sigProfilerAssignment(inputdir,project_name):
    #set directories and paths to signatures and samples
    print(" In sig profiler sigProfilerAssignment  Function")
    '''
    dir_inp     = '../data/'
    samples_file     = dir_inp+outdir+"/output/SBS/Allsamples.SBS96.all"
    output_dir      =  dir_inp+outdir+"input_SigProfiler/"
    '''
    output_dir = inputdir
    #Analysis of SP Assignment
    Analyze.cosmic_fit( inputdir,
                    output_dir,
                    genome_build=build,
                    cosmic_version=3.4,
                    verbose=True,
                    context_type="96",
                    input_type="vcf",
                    collapse_to_SBS96=True,
                    make_plots=True,
                    exclude_signature_subgroups=None,
                    sample_reconstruction_plots=True,
                    exome=True, ## disable if running onWGS
                    export_probabilities=True)

if __name__=="__main__":
    input_dir = "/path/to/input/dir/with/vcfs"
    sbs_dir=input_dir;sbs_project_name="TAN-WES-no40"
    sigProfilerAssignment(sbs_dir,sbs_project_name) # for sbs-singlet base substitutions
