import os
import csv

def parse_df(df_dir,ref_dict):
    #parse defence finder
    df_db_dict = {}
    for df_file in os.listdir(os.path.join(df_dir,"systems")):
        if df_file.endswith(".tsv"):
            accession = df_file.replace("_defense_finder_systems.tsv","")
            df_db_dict[accession] = []
            df = open(os.path.join(df_dir,"systems",df_file), mode = "r", encoding='utf-8-sig')
            df_list = csv.reader(df,quotechar='"',delimiter="\t",quoting=csv.QUOTE_ALL, skipinitialspace=True)

            start = True
            for line_list in df_list:
                if start:
                    for i in range(len(line_list)):
                        if line_list[i] == "type":
                            system_index = i
                        elif line_list[i] == "subtype":
                            subtype_index = i
                        elif line_list[i] == "protein_in_syst":
                            genes_index = i
                        elif line_list[i] == "name_of_profiles_in_sys":
                            gene_names_index = i

                    start = False
                else:
                    #get relevant values
                    system = line_list[system_index]
                    if system in ref_dict:
                        if system != "skip":
                            system = ref_dict[system]
                    else:
                        print ("defence finder ref error:", system)
                    subtype = line_list[subtype_index]
                    gene_list = ";".join(line_list[genes_index].split(","))
                    gene_names = line_list[gene_names_index]

                    df_db_dict[accession] = df_db_dict[accession]+[[system,subtype,gene_list,gene_names]]
    return df_db_dict
