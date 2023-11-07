#!/usr/bin/env python3

#packages
import argparse
import timeit
import time
import os

#other code
from parse_padloc_and_defensefinder_output import parse_padloc_and_defensefinder_output
from run_dependencies import run_dependencies

#times function
start = timeit.default_timer()

header = "Predict immune systems from nucleic acid or protein fasta file."

#TODO: check runable with both single file or entire dir
#TODO: check installation guide
#TODO: per system output
#TODO: remove meta
#TODO: make padloc forcibly overwrite
#TODO: make runable from path when single file is used
#TODO: Make file possible for single file input
#TODO: trial run split_big_files_into
#TODO: Does defense finder now run with spaces in directory?
#TODO: delete prot_dir thing
#TODO: remove type

##################################
#Set commandline function options#
##################################

parser = argparse.ArgumentParser(description=header)
parser.add_argument("-i", "--input_db",help="File or directory to be analysed in fasta format.")
parser.add_argument("-o", "--output_dir",help="Location to save output. Default: results",required = False) 
parser.add_argument("-c", "--high_confidence",help="When high confidence output (i.e. found by both PADLOC and DefenseFinder is required)",action="store_true", required = False)
parser.add_argument("-p", "--protein",help="File format of input file is not nucleic acids, but protein. In this case, gff needs to be used.", action="store_true", required = False)
parser.add_argument("-gff",help="GFF file of proteins.", required = False)
parser.add_argument("-cpu",help="How many cores to use for DefenseFinder and PADLOC. Default: max.", required = False) 
parser.add_argument("--no_padloc", help="Don't run PADLOC (e.g. because you don't want the output or if already run and stored in right location).",action="store_true", required = False)
parser.add_argument("--no_defensefinder", help="Don't run defensefinder (e.g. because you don't want the output or if already run and stored in right location).",action="store_true", required = False)
parser.add_argument("-v", "--verbose",help="Verbose output.",action="store_true", required = False)
parser.add_argument("-q", "--quiet",help="Quiet output.",action="store_true", required = False)
parser.add_argument("-m", "--meta",help="When running on a metagenome. Uses metagenome version of prodigal.",action="store_true", required = False)
parser.add_argument("-f","--force", help="Forcibly overwrites prodigal, defensefinder and PADLOC files that already exist.",action="store_true", required = False)
parser.add_argument("--split_big_files_into", help="Give a size to splits big files (e.g. very big metagenomes) into for DefenseFinder, as current DefenseFinder version can hang on them. Recommended size: 8000000. If nothing is provided, files will not be split.", required = False)

args = parser.parse_args()

#############
#parse input#
#############
input_db = args.input_db
if os.path.isfile(input_db):
    input_dir = os.path.dirname(input_db)
    input_db = [os.path.basename(input_db)]
elif os.path.isdir(input_db):
    input_dir = input_db[:]
    input_db = os.listdir(input_db)
else:
    raise Exception(input_db+" is not a file or directory.")

if args.output_dir:
    output_dir = args.output_dir
else:
    output_dir = "results"
try:
    os.mkdir(output_dir)
except:
    pass

if not args.protein:
    prot_dir = os.path.join(output_dir,"proteins")
    try:
        os.mkdir(prot_dir)
    except:
        pass
    gff_dir = os.path.join(output_dir,"gff")
    try:
        os.mkdir(gff_dir)
    except:
        pass
else:
    gff_dir = args.gff
    prot_dir = input_dir

if not args.cpu:
    cpu = str(os.cpu_count())
else:
    cpu = args.cpu

if args.meta:
    from Bio import SeqIO
    meta = True
    split_file_dir = os.path.join(output_dir,"split_files")
    try:
        os.mkdir(split_file_dir)
    except:
        pass
else:
    meta = False

############################
#make remaining result dirs#
############################
padloc_output_dir = os.path.join(output_dir,"padloc_output")
df_output_dir = os.path.join(output_dir,"defense_finder_output")
combined_output_dir = os.path.join(output_dir,"combined_output")

df_systems_dir = os.path.join(df_output_dir,"systems")
df_genes_dir = os.path.join(df_output_dir,"genes")
df_hmmer_dir= os.path.join(df_output_dir,"hmmer")

try:
    os.mkdir(combined_output_dir)
    if not args.no_padloc:
        os.mkdir(padloc_output_dir)
    if not args.no_defensefinder:
        os.mkdir(df_output_dir)
except:
    pass
try:
    if not args.no_defensefinder:
        os.mkdir(df_systems_dir)
        os.mkdir(df_genes_dir)
        os.mkdir(df_hmmer_dir)
except:
    pass

if args.split_big_files_into:
    split_big_file_number = int(args.split_big_files_into)
else:
    split_big_file_number = 0
    
################
#run programmes#
################
accessions = []
for file in input_db:
    accessions.append(os.path.splitext(file)[0])
    run_dependencies (file, input_dir, prot_dir, gff_dir, df_output_dir,
                      df_systems_dir, df_genes_dir, df_hmmer_dir,
                      padloc_output_dir, cpu, args.protein, args.no_defensefinder,
                      args.no_padloc, meta, time, split_big_file_number, args.force, args.verbose)

##############
#parse output#
##############
parse_padloc_and_defensefinder_output(accessions,output_dir,padloc_output_dir,df_output_dir,combined_output_dir,high_confidence=args.high_confidence,output_file="output_summary.csv",ref_file="immune_system_list_reference.csv")


stop = timeit.default_timer()
print('Analysis finished. Processing time: %.2fs, %.2fm' % (stop - start, (stop - start) / 60))

