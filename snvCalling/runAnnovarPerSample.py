#!/usr/bin/python
__author__="Pushpa Itagi"
__purpose__=" Run annovar on these files"
_comments_="18May2023 script to process the mutect2(or similar callers) vcf output and annotate them by ANNOVAR *******"

#Load modules
import os
import sys
import argparse


#please load the module annovar if using this script directly by typing : ml annovar on the shell
#******* Change this if needed *******
#Refers to the file extensions that are created after running mutation callers
file_pattern1 = ".vcf.gz" #file pattern to search for
###file_pattern2 = "filtered" #skip these filtered vcf files as they have additional mutations which are not pass
input_path_to_humandb = "/path/to/annovar/database/HG38" #databases for annovar files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process files.')
    parser.add_argument('--input_vcf_file_path',help='please enter a path to the vcf file',required=True)
    args = parser.parse_args()
    input_vcf_file_path = args.input_vcf_file_path
    out_file = input_vcf_file_path.split(file_pattern1)[0]
    #******* Change this if needed *******
    command1 = "table_annovar.pl "+input_vcf_file_path+" "+input_path_to_humandb+" -buildver hg38 -out "+out_file+" -remove -protocol refGene,cytoBand,esp6500siv2_all,avsnp144,avsnp150,ALL.sites.2015_08,gnomad_genome,gnomad_exome,exac03,clinvar_20190305,dbnsfp33a,gnomad312_genome,cosmic70 -operation g,r,f,f,f,f,f,f,f,f,f,f,f  -nastring . -polish -vcfinput"
    ##print(command1)
    print(" Annotated file for input file -",input_vcf_file_path," is here: ",out_file)
    print("\n")
    os.system(command1) #execute the command
