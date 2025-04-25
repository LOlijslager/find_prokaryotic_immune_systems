# find_prokaryotic_immune_systems
Automated script that runs DefenseFinder and PADLOC and merges their output into one database

# When using this software, please reference
- Olijslager, L. H., Weijers, D., & Swarts, D. C. (2024). Distribution of specific prokaryotic immune systems correlates with host optimal growth temperature. NAR Genomics and Bioinformatics, 6(3), lqae105.

- Tesson, F., Hervé, A., Mordret, E., Touchon, M., d’Humières, C., Cury, J., & Bernheim, A. (2022). Systematic and quantitative view of the antiviral arsenal of prokaryotes. Nature communications, 13(1), 2561.

- Abby, S. S., Néron, B., Ménager, H., Touchon, M., & Rocha, E. P. (2014). MacSyFinder: a program to mine genomes for molecular systems with an application to CRISPR-Cas systems. PloS one, 9(10), e110726.

- Payne, L. J., Todeschini, T. C., Wu, Y., Perry, B. J., Ronson, C. W., Fineran, P. C., ... & Jackson, S. A. (2021). Identification and classification of antiviral defence systems in bacteria and archaea with PADLOC reveals new system types. Nucleic Acids Research, 49(19), 10868-10878.

If also running the custom systems, include:
van den Berg, D. F., Costa, A. R., Esser, J. Q., Stanciu, I., Geissler, J. Q., Zoumaro-Djayoon, A. D., ... & Brouns, S. J. (2024). Bacterial homologs of innate eukaryotic antiviral defenses with anti-phage activity highlight shared evolutionary roots of viral defenses. Cell Host & Microbe, 32(8), 1427-1443.

# installing the software using conda
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

# options
```sh
$find_immune_systems.py [-h] [-i INPUT_DB] [-o OUTPUT_DIR] [-c] [-p] [-gff GFF] [--custom_system_dir CUSTOM_SYSTEM_DIR]     
                              [--conflict_winner CONFLICT_WINNER] [-cpu CPU] [--no_padloc] [--no_defensefinder]                   
                              [--reference_file REFERENCE_FILE] [--skip_DS_candidates] [--cut_no_island CUT_NO_ISLAND] [-v] [-q]  
                              [-m] [-f] [--split_big_files_into SPLIT_BIG_FILES_INTO] 

Predict immune systems from nucleic acid or protein fasta file. Version 0.3.0.

options:
  -h, --help            show this help message and exit
  -i, --input_db INPUT_DB
                        File or directory to be analysed in fasta format.
  -o, --output_dir OUTPUT_DIR
                        Location to save output. Default: results
  -c, --high_confidence
                        When high confidence output (i.e. found by both PADLOC and DefenseFinder is required)
  -p, --protein         File format of input file is not nucleic acids, but protein. In this case, gff needs to be used.
  -gff GFF              GFF file of proteins.
  --custom_system_dir CUSTOM_SYSTEM_DIR
                        Directory with custom HMM-profiles to search for. Also requires systems.csv with the system compositions
                        and bitscores, formatted as in the systems.csv example file. Warning: these systems will only be checked
                        for perfect overlap.
  --conflict_winner CONFLICT_WINNER
                        In case of conflict between identified systems, should DefenseFinder (DF) or PADLOC (PL) overrule?
                        Default: DF.
  -cpu CPU              How many cores to use for DefenseFinder and PADLOC. Default: max.
  --no_padloc           Don\'t run PADLOC or don\'t include the results if already run).
  --no_defensefinder    Don\'t run defensefinder or don\'t include the results if already run.
  --reference_file REFERENCE_FILE
                        Provide filepath for the systems reference file. If not, reference file is searched for in the current
                        working directory and the directory this script is located in.
  --skip_DS_candidates  Skips defence system candidates identified by PADLOC (i.e. the HEC and PDC systems).
  --cut_no_island CUT_NO_ISLAND
                        Removes matches that are not found within 20 proteins of another immune system. Provide system names
                        seperated by \',\' (e.g. SoFIC,PD-T4-7). Provide \'All\' to do this for all systems (Keep in mind this will
                        drastically decrease your output (e.g. CRISPR-Cas systems usually aren\'t)).
  -v, --verbose         Verbose output.
  -q, --quiet           Quiet output.
  -m, --meta            When running on a metagenome. Uses metagenome version of prodigal.
  -f, --force           Forcibly overwrites prodigal, defensefinder and PADLOC files that already exist.
  --split_big_files_into SPLIT_BIG_FILES_INTO
                        Give a size to splits big files (e.g. very big metagenomes) into for DefenseFinder, as current
                        DefenseFinder version can hang on them. Recommended size: 8000000. If nothing is provided, files will not
                        be split.
```
# example usage

start a session using this command
```sh 
conda activate prok_def
```

#Example 1: Run on genome files saved in a directory
```sh 
./find_immune_systems.py -i example_usage/example_genomes
```

#Example 2: Run on specfific genome file
```sh 
./find_immune_systems.py -i example_usage/example_genomes/GCF_000264495.1_ASM26449v1_genomic.fna
```

#Example 3: Use on protein sequence, bypassing translation by prodigal
```sh 
./find_immune_systems.py -p -i example_usage/example_proteomes/ -gff example_usage/example_gff/
```

#Example 4: Our recommended settings, which include some added systems not currently in PADLOC or DefenceFinder, and which predicts only proven defence system proteins and which reduces predictions that appear to be over-predicted in our hands by requiring them to be found in defence islands.
```sh 
./find_immune_systems.py -i example_usage/example_genomes --custom_system_dir custom_systems --skip_DS_candidates --cut_no_island SoFIC,PD-T4-7
```

#end a session using this command
```
sh conda deactivate
```

# output
The programme creates:

- combined_output: directory with, for each sequence file provided, a csv file noting, for each immune system found, the system family and subtype identified by DefenseFinder and PADLOC, the operon (noting the gene identifiers of all proteins included in the system), how the genes are called by DefenseFinder and how they're called by PADLOC.
- defense_finder_output: directory with the output provided by DefenseFinder
- gff: directory with gff files as created by Prodigal
- padloc_output: directory with the output provided by PADLOC
- proteins: directory with protein files as created by Prodigal
- out_summary.csv: comma seperated file containing a summary of the combined output, noting down how often immune system families were found in what sequence file.

# Troubleshooting
- Macsyfinder (DefenseFinder dependable) can have trouble running from on a directory in spaces, even if the spaces aren't included in the given path. Make sure none of the directories you're using include spaces.
- At the time of release the Macsyfinder (DefenseFinder dependable) version DefenseFinder is working with can have trouble with some big files, eternally hanging on the command. If you're running e.g. big metagenome files, pass a number into --split_big_files_into, which cleaves the files into chunks of this size. In our experience, 8000000 works well.
- If the code gives "defence finder ref error: [system_name]" or "PADLOC ref error: [system_name]", this is an indication that a found system is not in the reference file (immune_system_list_reference.csv). This can be either because DefenseFinder or PADLOC, respectively, got updated since the making of this code, or the system is relatively uncommon. Simply add the system to the reference file and rerun the code to fix this issue. 
- If the code gives "WARNING: same system identified twice. Might need to update reference file", this indicates that the same operon is identified by as a different immune system family by DefenseFinder and PADLOC. The first possible cause for this is either the system families are synonymous, but because of an update or in the case of a less common system, this isn't indicated in the reference file (immune_system_list_reference.csv) yet. Simply update the reference file and rerun the code. The second explanation is that one of the two is wrong. In this case, disregard the incorrect identification before proceeding with your analysis.
- If you are running miniconda and your receive errors concerning r-rmarkdown, try running ```sh conda install -c r r-rmarkdown```
