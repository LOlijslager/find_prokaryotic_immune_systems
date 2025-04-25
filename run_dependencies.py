import os

from run_prodigal import run_prodigal
from run_defense_finder import run_defense_finder
from run_padloc import run_padloc
from run_hmmer import run_hmmer

def run_dependencies (file, input_dir, prot_dir, gff_dir, df_output_dir, df_systems_dir, df_genes_dir,df_hmmer_dir,
                      padloc_output_dir, HMMER_output_dir, cpu, custom_system_dir, protein,
                      no_defensefinder, no_padloc,meta, time, split_big_file_number, force, verbose):
    #################
    #predict protein#
    #################
    if not protein:
        genome_file = file
        genome_dir = input_dir
        prot_file,gff_file = run_prodigal(genome_file,genome_dir,prot_dir,gff_dir,meta,force,verbose)
        if verbose:
            print (file, " proteins being predicted using prodigal, time: ", time.strftime("%H:%M:%S", time.localtime()))
        if not os.path.isfile(os.path.join(prot_dir,prot_file)):
            raise Exception(genome_file+" raised an exception in prodigal. Is this file a valid nucleic acid file?")
    else:
        prot_file = file
        base_name = os.path.splitext(prot_file)[0]
        gff_file = base_name+".gff"
    if not os.stat(os.path.join(prot_dir,prot_file)).st_size == 0:
        ####################
        #run defense finder#
        ####################
        if not no_defensefinder:
            if verbose:
                print (file, " immune systems being predicted using defensefinder, time: ", time.strftime("%H:%M:%S", time.localtime()))
            run_defense_finder(prot_file,prot_dir,df_output_dir,df_systems_dir,df_genes_dir,df_hmmer_dir,split_big_file_number,cpu,force,verbose)

        ############
        #run PADLOC#
        ############
        if not no_padloc:
            if verbose:
                print (file, " immune systems being predicted using PADLOC time: ", time.strftime("%H:%M:%S", time.localtime()))
            run_padloc(prot_file,prot_dir,gff_file,gff_dir,padloc_output_dir,cpu,force,verbose)

        ##########################################################
        #run HMMER on any additional systems provided by the user#
        ##########################################################
        if custom_system_dir:
            if verbose:
                print (file, " additional immune systems being predicted using HMMER time: ", time.strftime("%H:%M:%S", time.localtime()))
            run_hmmer(custom_system_dir,custom_system_dir,prot_file,prot_dir,HMMER_output_dir,cpu,force,verbose)
            
