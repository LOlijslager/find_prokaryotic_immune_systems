import os
import subprocess
import shutil
    
def run_defense_finder(prot_file,prot_dir,output_dir,systems_dir,genes_dir,hmmer_dir,split_big_file_number,cpu,force,verbose):
    prot_file_path = os.path.join(prot_dir,prot_file)
    base_name = os.path.splitext(prot_file)[0]

    final_systems_file = os.path.join(systems_dir,base_name+"_defense_finder_systems.tsv")
    final_genes_file = os.path.join(genes_dir,base_name+"_defense_finder_genes.tsv")
    final_hmmer_file = os.path.join(hmmer_dir,base_name+"_defense_finder_hmmer.tsv")

    if force or (not os.path.isfile(final_systems_file)):
        systems_file = os.path.join(output_dir,base_name+"_defense_finder_systems.tsv")
        genes_file = os.path.join(output_dir,base_name+"_defense_finder_genes.tsv")
        hmmer_file = os.path.join(output_dir,base_name+"_defense_finder_hmmer.tsv")
            
        if (split_big_file_number!=0) and (os.path.getsize(os.path.join(prot_dir,prot_file)) > split_big_file_number): #current version of defensefinder uses a version of macsyfinder which cannot handle big files
            split_file_dir = os.path.join(output_dir, "split_genome_files")
            split_files = []
            if not os.path.isdir(split_file_dir):
                os.mkdir(split_file_dir)
            with open(os.path.join(prot_dir,prot_file)) as myfile:
                lines = myfile.readlines()
                line_number = len(lines)/(os.path.getsize(os.path.join(prot_dir,prot_file))/8000000)+1
                n = 0
                split_file_number = 1
                split_file_path = os.path.join(split_file_dir,base_name+"_part_"+str(split_file_number)+".faa")
                split_file = open(split_file_path, "w")
                split_files.append(split_file_path)
                for line in lines:
                    if n < line_number:
                        n += 1
                    else:
                        if line.startswith(">"):
                            n = 0
                            split_file.close()
                            split_file_number+=1
                            split_file_path = os.path.join(split_file_dir,base_name+"_part_"+str(split_file_number)+".faa")
                            split_file = open(split_file_path, "w")
                            split_files.append(split_file_path)
                    split_file.write(line)
                split_file.close()
            start_1 = True
            for split_file_path in split_files:
                cmd = 'defense-finder run "'+split_file_path+'" -o "'+output_dir+'" --workers '+cpu
                if verbose:
                    print (cmd)
                subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL)
                
                with open(os.path.join(systems_dir,"systems_"+base_name+".tsv"),"a+") as o_f:
                    with open(systems_file) as i_f:
                        start_2 = True
                        for line in i_f:
                            if start_1 == True:
                                o_f.write(line)
                            elif start_2 == False:
                                o_f.write(line)
                            start_2 = False
                os.remove(systems_file)
                            
                with open(os.path.join(genes_dir,"genes_"+base_name+".tsv"),"a+") as o_f:
                    with open(genes_file) as i_f:
                        start_2 = True
                        for line in i_f:
                            if start_1 == True:
                                o_f.write(line)
                            elif start_2 == False:
                                o_f.write(line)
                            start_2 = False
                os.remove(genes_file)

                with open(os.path.join(hmmer_dir,"hmmer_"+base_name+".tsv"),"a+") as o_f:
                    with open(hmmer_file) as i_f:
                        start_2 = True
                        for line in i_f:
                            if start_1 == True:
                                o_f.write(line)
                            elif start_2 == False:
                                o_f.write(line)
                            start_2 = False
                            
                start_1 = False
                os.remove(hmmer_file)

                os.remove(split_file_path)
                os.remove(split_file_path+".idx")
            os.rmdir(split_file_dir)
        else:
            cmd = 'defense-finder run "'+prot_file_path+'" -o "'+output_dir+'" --workers '+cpu
            if verbose:
                print (cmd)
            subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL)
        
            shutil.move(systems_file, final_systems_file)
            shutil.move(genes_file, final_genes_file)
            shutil.move(hmmer_file, final_hmmer_file)
            os.remove(prot_file_path+".idx")
        
