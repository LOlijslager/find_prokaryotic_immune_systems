import os
import subprocess


def run_hmmer(custom_system_dir,HMMER_systems_dir,prot_file,prot_dir,output_dir,cpu,force,verbose):
    prot_file_path = os.path.join(prot_dir,prot_file)
    base_name = os.path.splitext(prot_file)[0]

    #run HMMER
    for file in os.listdir(HMMER_systems_dir):
        if file.endswith (".hmm"):
            hmm_profile = file
            hmm_file_path = os.path.join(HMMER_systems_dir,hmm_profile)
            hmm_name = os.path.splitext(hmm_profile)[0]
            output_file_path = os.path.join(output_dir,base_name+"_"+hmm_name+".out")
            
            if force or (not os.path.isfile(output_file_path)): #Only run if required
                cmd = "hmmsearch --noali --domtblout "+output_file_path+" --cpu "+cpu+" "+hmm_file_path+" "+prot_file_path
                if verbose:
                    print (hmm_profile, cmd)
                subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL)
