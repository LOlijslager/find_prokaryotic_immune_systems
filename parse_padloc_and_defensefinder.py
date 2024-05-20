import os
import csv

##################################
#Counting how many hits DF&PADLOC#
##################################

def parse_padloc_and_defensefinder_output(accessions,output_dir,padloc_dir,df_dir,combined_output_dir,
                                          conflict_winner,quiet,high_confidence,
                                          output_file="output_summary.csv",ref_file="immune_system_list_reference.csv",
                                          warning_file = "warnings.txt"):
    
    warning_file = os.path.join(output_dir,"warnings.txt")
    ref_dict = load_reference_file(ref_file)
    padloc_db_dict = parse_padloc(padloc_dir,ref_dict)
    df_db_dict = parse_df(df_dir,ref_dict)

    db_dictionary,defence_system_list=combine_df_padloc_output(accessions,padloc_db_dict,df_db_dict,combined_output_dir,high_confidence,warning_file,conflict_winner,quiet)
    write_summary_file(output_dir,db_dictionary,defence_system_list,output_file)

def load_reference_file(ref_file):
    #Load file with all the name differences between PADLOC and Defence Finder
    with open(ref_file) as ref_list_file:
        ref_list = csv.reader(ref_list_file)
        start = True
        shared_systems = [] #list with systems both softwares look for
        ref_dict={} #dictionary with all reference values to equalise the systems
        for line in ref_list:
            if not start:
                if "" not in line:
                    shared_systems.append(line[0])
                if line[1] != "":
                    if "/" not in line[1]: #multiple names
                        ref_dict[line[1]] = line[0]
                    else:
                        names = line[1].split("/")
                        for name in names:
                            ref_dict[name] = line[0]
                if line[2] != "":
                    if "/" not in line[2]: #multiple names
                        ref_dict[line[2]] = line[0]
                    else:
                        names = line[2].split("/")
                        for name in names:
                            ref_dict[name] = line[0]
            start = False
    return ref_dict

def parse_padloc(padloc_dir,ref_dict):
    ###parse PADLOC
    padloc_db_dict = {}
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

def parse_df(df_dir,ref_dict):
    #parse defence finder
    df_db_dict = {}
    for df_file in os.listdir(os.path.join(df_dir,"systems")):
        if df_file.endswith(".tsv"):
            accession = df_file.replace("systems_","").replace(".tsv","")
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

def combine_df_padloc_output(accessions,padloc_db_dict,df_db_dict,combined_output_dir,
                             high_confidence,warning_file,conflict_winner,quiet):
    #How many systems are found in total
    defence_system_list = [] #all system families for summary file
    db_dictionary = {} #all systems for summary file
    
    for accession in accessions:
        combined_dict = {}
        db_dictionary[accession] = []

        PL_dict = {}
        if accession in padloc_db_dict.keys():
            for system in padloc_db_dict[accession]:
                PL_dict[";".join(system)] = system

        DF_dict = {}
        if accession in df_db_dict.keys():
            for system in df_db_dict[accession]:
                DF_dict[";".join(system)] = system

        PL_dict = remove_internal_conflict(PL_dict, DF_dict, "PL")
        DF_dict = remove_internal_conflict(DF_dict, PL_dict, "DF")

        if accession in padloc_db_dict.keys():
            for system in PL_dict.keys():
                system_db = PL_dict[system]
                combined_dict = test_system_for_overlap(system_db, combined_dict, "PL", DF_dict, conflict_winner, accession, quiet, warning_file)

        if accession in df_db_dict.keys():
            for system in DF_dict.keys():
                system_db = DF_dict[system]
                combined_dict = test_system_for_overlap(system_db, combined_dict, "DF", PL_dict, conflict_winner, accession, quiet, warning_file)

        with open(os.path.join(combined_output_dir,accession+".csv"), "w") as output_file:
            output_file.write("accession,system family DefenseFinder,system family PADLOC,subtype DefenseFinder,subtype PADLOC,genes,genes DefenseFinder,genes PADLOC\n")
            
            if high_confidence: #in case high confidence (i.e. found by both PADLOC and DF) is desired
                for genes in combined_dict.keys():
                    if ("DF" not in genes.keys()) or ("PL" not in genes.keys()):
                        combined.pop(genes)
                        
            #combine data
            for genes in combined_dict.keys():
                info = combined_dict[genes]

                #make summary file for accession
                system_PL = "N.A."
                system_DF = "N.A."
                subtype_PL = "N.A."
                subtype_DF = "N.A."
                genes_DF = "N.A."
                genes_PL = "N.A."

                if "DF" in info.keys():
                    system_DF = info["DF"][0]
                    subtype_DF = info["DF"][1]
                    genes_DF = info["DF"][3].replace(",",";")
                    if system_DF not in defence_system_list:
                        defence_system_list.append(system_DF)

                if "PL" in info.keys():
                    system_PL = info["PL"][0]
                    subtype_PL = info["PL"][1]
                    genes_PL = "N.A."
                    genes_PL = info["PL"][3].replace(",",";")
                    if system_PL not in defence_system_list:
                        defence_system_list.append(system_PL)

                output_file.write("%s,%s,%s,%s,%s,%s,%s,%s\n"%(accession,system_DF,system_PL,subtype_DF,subtype_PL,genes,genes_DF,genes_PL))

                #make db for summary file
                if system_PL == system_DF:
                    db_dictionary[accession].append(system_PL)
                else:
                    if system_PL != "N.A.":
                        db_dictionary[accession].append(system_PL)
                    if system_DF != "N.A.":
                        db_dictionary[accession].append(system_DF)
    defence_system_list.sort()
    return db_dictionary,defence_system_list

def remove_internal_conflict(db_from_self, db_from_other, testing):
    new_db_from_self = dict()
    #test internal conflict (systems with differen families assigned but the exact same genes within the DF or PL database)
    removed_list = []
    for system1 in db_from_self.keys():
        remove = False
        internal_conflict = False
        for system2 in db_from_self.keys():
            genes_1 = db_from_self[system1][2]
            genes_2 = db_from_self[system2][2]
            genes_set_1 = set(genes_1.split(";"))
            genes_set_2 = set(genes_2.split(";"))
            if genes_set_1 == genes_set_2:
                if system1 != system2:
                    internal_conflict = True
                    db = db_from_self[system1]
                    quiet = "not_main"
                    accession = "not_main"
                    warning_file = "not_main"
                    found_by_both_b, part_of_a_whole, part_of_a_whole_longest, overlap_not_part_of_a_whole, different_family_same_system, overlap_not_part_of_a_whole, different_family, no_conflict = find_overlap(db, genes_set_1, db_from_other, accession, testing, quiet, warning_file)
                    if genes_set_1 in removed_list: #previously merged
                        remove = True
                    elif db[0] == db_from_self[system2][0]: #same system family, differen subsystem: merge
                        new_db_from_self[system1] = db[0], db[1]+"/"+db_from_self[system2][1], db[2], db[3]+"/"+db_from_self[system2][3]
                        removed_list.append(genes_set_1)
                    elif found_by_both_b:
                        new_db_from_self[system1] = db_from_self[system1] #keep unchanged
                    elif different_family_same_system:
                        remove = True #do not test
                    else: #unique find, merge
                        new_db_from_self[system1] = db[0]+"/"+db_from_self[system2][0], db[1]+"/"+db_from_self[system2][1], db[2], db[3]+"/"+db_from_self[system2][3]
                        removed_list.append(genes_set_1)
        if internal_conflict == False:
            new_db_from_self[system1] = db_from_self[system1]
    return new_db_from_self

def find_overlap(db, genes_set, db_from_other, accession, testing, quiet, warning_file):
    no_conflict = False
    found_by_both_b = False
    part_of_a_whole_longest = False
    part_of_a_whole = False
    overlap_not_part_of_a_whole = False
    different_family_same_system = False
    different_family = False
    internal_conflict = False
    overlap = False

    #test for conflict between PL and DF
    conflicting_systems = []
    for system in db_from_other.keys():
        genes_from_other = db_from_other[system][2]
        genes_from_other_set = set(genes_from_other.split(";"))
        combined_genes_set = genes_set.union(genes_from_other_set)
        if len (combined_genes_set) < (len(genes_set)+len(genes_from_other_set)):
            overlap = True
            conflicting_systems.append(db_from_other[system]) #define which system(s) conflict

    if overlap: #Does the system have overlap with a system from the other?
        for conflicting_system in conflicting_systems:
            genes_from_other_set = set(conflicting_system[2].split(";"))
            combined_genes_set = genes_set.union(genes_from_other_set)
            if db[0] == conflicting_system[0]: #Do these systems share the same system family?
                if genes_set == genes_from_other_set: #Do they have exactly the same genes?
                    found_by_both_b = True
                else:
                    if (combined_genes_set == genes_set) or (combined_genes_set == genes_from_other_set): #is one contained in the other?
                        part_of_a_whole = True
                        if combined_genes_set == genes_set: #keep longest genes
                            part_of_a_whole_longest = True
                    else:
                        overlap_not_part_of_a_whole = True

            else:
                if genes_set == genes_from_other_set: #do they have exactly the same genes?
                    different_family_same_system = True
                    if quiet != "not_main":
                        if not quiet:
                            if testing == "DF":
                                print ("Warning: same system identified twice. Might need to update reference file. Pertaining: sequence %s, DefenseFinder %s, PADLOC %s, genes %s"%(accession, db[0], conflicting_system[0], db[2]))
                        with open (warning_file, "a+") as w_f:
                            w_f.write("Warning: same system identified twice. Might need to update reference file. Pertaining: sequence %s, DefenseFinder %s, PADLOC %s, genes %s \n"%(accession, db[0], conflicting_system[0], db[2]))
                else:
                    overlap_not_part_of_a_whole = True
                    different_family = True
    else:
        no_conflict = True

    return found_by_both_b, part_of_a_whole, part_of_a_whole_longest, overlap_not_part_of_a_whole, different_family_same_system, overlap_not_part_of_a_whole, different_family, no_conflict


def test_system_for_overlap(db, combined_dict, testing, db_from_other, conflict_winner, accession, quiet, warning_file):
    genes = db[2]
    genes_set = set(genes.split(";"))
    conflicting_systems = []

    found_by_both_b, part_of_a_whole, part_of_a_whole_longest, overlap_not_part_of_a_whole, different_family_same_system, overlap_not_part_of_a_whole, different_family, no_conflict = find_overlap(db, genes_set, db_from_other, accession, testing, quiet, warning_file)

    tester_wins = False
    conflict_winner_wins = False
    if found_by_both_b:
        tester_wins = True
    elif no_conflict:
        tester_wins = True
    elif part_of_a_whole:
        if part_of_a_whole_longest:
            tester_wins = True
    elif different_family_same_system or overlap_not_part_of_a_whole:
        conflict_winner_wins = True
    
    if conflict_winner_wins:
        if testing == "PL" and conflict_winner == "PL":
            if genes in combined_dict.keys():
                combined_dict[genes]["PL"] = db
            else:
                combined_dict[genes] = {"PL":db}
        elif testing == "DF" and conflict_winner == "DF":
            if genes in combined_dict.keys():
                combined_dict[genes]["DF"] = db
            else:
                combined_dict[genes] = {"DF":db}
    elif tester_wins:
        if testing == "PL":
            if genes in combined_dict.keys():
                combined_dict[genes]["PL"] = db
            else:
                combined_dict[genes] = {"PL":db}
        else:
            if genes in combined_dict.keys():
                combined_dict[genes]["DF"] = db
            else:
                combined_dict[genes] = {"DF":db}
        
    return combined_dict

def write_summary_file(output_dir,db_dictionary,defence_system_list,output_file):
    #export relevant data to csv
    db_file = os.path.join(output_dir,output_file)
    with open(db_file, "w") as output_file:
        defence_systems = ",".join(defence_system_list)
        header = "accession,%s\n"%(defence_systems)
        output_file.write(header)
        for accession in db_dictionary:
            acc_defence_systems = ""
            if accession in db_dictionary:
                for defence_system in defence_system_list:
                    if defence_system in db_dictionary[accession]:
                        acc_defence_systems += str(db_dictionary[accession].count(defence_system))+","
                    else:
                        acc_defence_systems += ","
            output_file.write('%s,%s\n'%(accession,acc_defence_systems))

if __name__ == "__main__":
    output_dir = "results"
    padloc_dir = os.path.join(output_dir,"padloc_output")
    df_dir = os.path.join(output_dir,"defense_finder_output")
    combined_output_dir = os.path.join(output_dir,"combined_output")
    high_confidence = False
    parse_padloc_and_defensefinder_output(output_dir,padloc_dir,df_dir,combined_output_dir,high_confidence)
