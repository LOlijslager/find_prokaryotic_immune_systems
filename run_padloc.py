import os
import subprocess

def run_padloc(prot_file,prot_dir,gff_file,gff_dir,output_dir,cpu,force,verbose):
    base_name = os.path.splitext(prot_file)[0]
    if force or (not os.path.isfile(os.path.join(output_dir,base_name+".domtblout"))):

        gff_file_path = os.path.join(gff_dir,gff_file)
        prot_file_path = os.path.join(prot_dir,prot_file)

        cmd = "padloc --faa "+prot_file_path+" --gff "+gff_file_path+" --outdir "+output_dir+" --fix-prodigal --cpu "+cpu
        if verbose:
            print (cmd)
            subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL)
        else:
            subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
