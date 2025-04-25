import os
import csv

def parse_padloc(padloc_dir,ref_dict,skip_DS_candidates):
    ###parse PADLOC
    padloc_db_dict = {}
    if os.path.exists(padloc_dir):
        for padloc_file in os.listdir(padloc_dir):
            if padloc_file.endswith(".csv"):
                accession = padloc_file.replace("_padloc.csv", "")
                padloc_db_dict[accession] = []
                padloc_file = open(os.path.join(padloc_dir,padloc_file), mode = "r", encoding='utf-8-sig')
                padloc_list = csv.reader(padloc_file,quotechar='"',delimiter=",",quoting=csv.QUOTE_ALL, skipinitialspace=True)

                start = True
                system_numbers = []
                systems_dict = {}
                for line_list in padloc_list:
                    if start:
                        for i in range(len(line_list)):
                            if line_list[i] == "system":
                                system_index = i
                            if line_list[i] == "system.number":
                                system_number_index = i
                            if line_list[i] == "target.name":
                                gene_index = i
                            if line_list[i] == "protein.name":
                                protein_name_index = i
                        start = False
                    else:
                        #get relevant values
                        system = line_list[system_index].split("_")[0]
                        subtype = line_list[system_index]
                        gene = line_list[gene_index]
                        protein_name = line_list[protein_name_index]

                        #check if this instance of the system was already added
                        system_number = line_list[system_number_index]
                        if subtype != "DRT_other": #these systems also exist as retron XIII and are identified solely as such by the online padloc version
                            if system_number not in systems_dict:
                                if system in ref_dict:
                                    system = ref_dict[system]
                                elif system.startswith("HEC") or system.startswith("PDC"):
                                    if skip_DS_candidates:
                                        system = "skip"
                                else:
                                    print ("padloc ref error:", system)
                                if system != "skip": #Like DMS with look like these are unfinished potential partial systems (or just plain duplicates), not necessary to investigate on full genomes like these
                                    systems_dict[system_number] = [system, subtype, [gene],[protein_name]]
                            else:
                                systems_dict[system_number][2] = systems_dict[system_number][2] + [gene]
                                systems_dict[system_number][3] = systems_dict[system_number][3] + [protein_name]
                #remove systems identified more than once that are identified as the same master family (happens for e.g. nhi and AbiO which are counted as the same system by DefenseFinder)
                temp = []
                res = dict()
                for key, val in systems_dict.items():
                    if val not in temp:
                        temp.append(val)
                        res[key] = val
                systems_dict = res.copy()
                #make final list of systems
                for key in systems_dict:
                    system = systems_dict[key][0]
                    subtype = systems_dict[key][1]
                    gene_list = ";".join(systems_dict[key][2])
                    protein_name_list = ";".join(systems_dict[key][3])
                    padloc_db_dict[accession]= padloc_db_dict[accession]+[[system,subtype,gene_list,protein_name_list]]
                padloc_file.close()
    return padloc_db_dict
