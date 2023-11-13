import os
import subprocess

def run_prodigal(genome_file,genome_dir,prot_dir,gff_dir,meta,force,verbose):
    base_name = os.path.splitext(genome_file)[0]
    prot_file = base_name+".faa"
    prot_file_path = os.path.join(prot_dir,base_name+".faa")
    gff_file = base_name+".gff"
    gff_file_path = os.path.join(gff_dir,base_name+".gff")
    genome_file = os.path.join(genome_dir,genome_file)

    if force or (not os.path.isfile(prot_file_path)):
        
        if meta:
            m = " -p meta"
        else:
            m = ""
        cmd = 'prodigal -q -f gff -a "'+prot_file_path+'" -i "'+genome_file+'" -o "'+gff_file_path+'"'+m
        if verbose:
            print (cmd)
        subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL)
    return prot_file, gff_file
