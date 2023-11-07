#find_prokaryotic_immune_systems
Automated script that runs DefenseFinder and PADLOC and merges their output into one database

#When using this software, please reference
- Olijslager, L.H., Weijers, D., & Swarts, D.C., Abundance of prokaryotic immune systems is correlated with host optimal growth temperature., unpublished.

- Tesson, F., Hervé, A., Mordret, E., Touchon, M., d’Humières, C., Cury, J., & Bernheim, A. (2022). Systematic and quantitative view of the antiviral arsenal of prokaryotes. Nature communications, 13(1), 2561.

- Abby, S. S., Néron, B., Ménager, H., Touchon, M., & Rocha, E. P. (2014). MacSyFinder: a program to mine genomes for molecular systems with an application to CRISPR-Cas systems. PloS one, 9(10), e110726.

- Payne, L. J., Todeschini, T. C., Wu, Y., Perry, B. J., Ronson, C. W., Fineran, P. C., ... & Jackson, S. A. (2021). Identification and classification of antiviral defence systems in bacteria and archaea with PADLOC reveals new system types. Nucleic Acids Research, 49(19), 10868-10878.

#installing the software using conda

```sh
conda create -n prok_def -c conda-forge -c bioconda -c padlocbio padloc
conda activate prok_def
padloc --db-update

pip3 install mdmparis-defense-finder
defense-finder update
```

If this command fails, try:
```sh
macsydata install -U -u --org mdmparis defense-finder-models
```
And, finally:
```sh
conda deactivate
```

#########
#options#
#########
```sh
$find_immune_systems.py [-h] [-i INPUT_DB] [-o OUTPUT_DIR] [-c] [-p] [-gff GFF] [-cpu CPU] [--no_padloc]
                              [--no_defensefinder] [-v] [-q] [-m] [-f] [--split_big_files_into SPLIT_BIG_FILES_INTO]

Predict immune systems from nucleic acid or protein fasta file.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DB, --input_db INPUT_DB
                        File or directory to be analysed in fasta format.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Location to save output. Default: results
  -c, --high_confidence
                        When high confidence output (i.e. found by both PADLOC and DefenseFinder is required)
  -p, --protein         File format of input file is not nucleic acids, but protein. In this case, gff needs to be
                        used.
  -gff GFF              GFF file of proteins.
  -cpu CPU              How many cores to use for DefenseFinder and PADLOC. Default: max.
  --no_padloc           Don't run PADLOC (e.g. because you don't want the output or if already run and stored in right
                        location).
  --no_defensefinder    Don't run defensefinder (e.g. because you don't want the output or if already run and stored
                        in right location).
  -v, --verbose         Verbose output.
  -q, --quiet           Quiet output.
  -m, --meta            When running on a metagenome. Uses metagenome version of prodigal.
  -f, --force           Forcibly overwrites prodigal, defensefinder and PADLOC files that already exist.
  --split_big_files_into SPLIT_BIG_FILES_INTO
                        Give a size to splits big files (e.g. very big metagenomes) into for DefenseFinder, as current
                        DefenseFinder version can hang on them. Recommended size: 8000000. If nothing is provided,
                        files will not be split.
```

###############
#example usage#
###############
#start a session using this command
```sh
conda activate prok_def
```

#Run on genome files saved in a directory
```sh
find_immune_systems -i genomes_directory
```

#Run on specfific genome file
```sh
find_immune_systems -i genomes_directory/genome.fna
```

#Use on protein sequence, bypassing translation by prodigal
```sh
find_immune_systems -p -i proteins_dir -gff gff_dir
```

#Recalculate combined summary file, without rerunning prodigal, DefenseFinder, and PADLOC
```sh
find_immune_systems -p -i proteins_dir -gff gff_dir --no_padloc --no_defensefinder
```

#end a session using these commands
```sh
conda deactivate
```

########
#output#
########
The programme creates:

- combined_output: directory with, for each sequence file provided, a csv file noting, for each immune system found, the system family and subtype identified by DefenseFinder and PADLOC, the operon (noting the gene identifiers of all proteins included in the system), how the genes are called by DefenseFinder and how they're called by PADLOC.
- defense_finder_output: directory with the output provided by DefenseFinder
- gff: directory with gff files as created by Prodigal
- padloc_output: directory with the output provided by PADLOC
- proteins: directory with protein files as created by Prodigal
- out_summary.csv: comma seperated file containing a summary of the combined output, noting down how often immune system families were found in what sequence file.

#################
#Troubleshooting#
#################
- Macsyfinder (DefenseFinder dependable) can have trouble running from on a directory in spaces, even if the spaces aren't included in the given path. Make sure none of the directories you're using include spaces.
- At the time of release the Macsyfinder (DefenseFinder dependable) version DefenseFinder is working with can have trouble with some big files, eternally hanging on the command. If you're running e.g. big metagenome files, pass a number into --split_big_files_into, which cleaves the files into chunks of this size. In our experience, 8000000 works well.
- If the code gives "defence finder ref error: [system_name]" or "PADLOC ref error: [system_name]", this is an indication that a found system is not in the reference file (immune_system_list_reference.csv). This can be either because DefenseFinder or PADLOC, respectively, got updated since the making of this code, or the system is relatively uncommon. Simply add the system to the reference file and rerun the code to fix this issue. 
- If the code gives "WARNING: same system identified twice. Might need to update reference file", this indicates that the same operon is identified by as a different immune system family by DefenseFinder and PADLOC. The first possible cause for this is either the system families are synonymous, but because of an update or in the case of a less common system, this isn't indicated in the reference file (immune_system_list_reference.csv) yet. Simply update the reference file and rerun the code. The second explanation is that one of the two is wrong. In this case, disregard the incorrect identification before proceeding with your analysis.
- If you are running miniconda and your receive errors concerning r-rmarkdown, try running ```sh conda install -c r r-rmarkdown```
