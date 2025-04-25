import os
import csv

from parse_defensefinder import parse_df
from parse_padloc import parse_padloc
from parse_hmmer import parse_custom_systems

#################################################
#Counting how many hits DF&PADLOC&custom_systems#
#################################################

def parse_output(accessions,output_dir,padloc_dir,df_dir,hmm_file_dir,combined_output_dir,
                 custom_system_dir,conflict_winner,quiet,skip_DS_candidates,
                 no_defensefinder,no_padloc,cut_no_island,high_confidence,ref_file,
                 output_file="output_summary.csv",warning_file = "warnings.txt"):
    
    warning_file = os.path.join(output_dir,"warnings.txt")
    
    if os.path.isfile("immune_system_list_reference.csv"):
        ref_file = "immune_system_list_reference.csv"
    elif os.path.isfile(os.path.join(os.path.dirname(os.path.realpath(__file__)),"immune_system_list_reference.csv")):
        ref_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),"immune_system_list_reference.csv")
    else:
        raise Exception("Cannot find immune system reference file.")
        
    ref_dict = load_reference_file(ref_file)
    if no_padloc:
        padloc_db_dict = {}
    else:
        padloc_db_dict = parse_padloc(padloc_dir,ref_dict,skip_DS_candidates)
    if no_defensefinder:
        df_db_dict = {}
    else:
        df_db_dict = parse_df(df_dir,ref_dict)
        
    if custom_system_dir != None:
        hmmer_db_dict = parse_custom_systems(accessions,custom_system_dir,hmm_file_dir)
    else:
        hmmer_db_dict = {}
    
    db_dictionary,defence_system_list=combine_output(accessions,padloc_db_dict,df_db_dict,hmmer_db_dict,combined_output_dir,
                                                     high_confidence,cut_no_island,warning_file,conflict_winner,quiet)
    
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


def combine_output(accessions,padloc_db_dict,df_db_dict,hmmer_db_dict,combined_output_dir,high_confidence,
                   cut_no_island,warning_file,conflict_winner,quiet):
    #How many systems are found in total
    defence_system_list = [] #all system families for summary file
    db_dictionary = {} #all systems for summary file
    
    for accession in accessions:
        combined_dict = {}
        db_dictionary[accession] = []

        #rewrite dictionaries
        PL_dict = {}
        if accession in padloc_db_dict.keys():
            for system in padloc_db_dict[accession]:
                PL_dict[";".join(system)] = system

        DF_dict = {}
        if accession in df_db_dict.keys():
            for system in df_db_dict[accession]:
                DF_dict[";".join(system)] = system

        #remove internal conflict
        PL_dict = remove_internal_conflict(PL_dict, DF_dict, "PL")
        DF_dict = remove_internal_conflict(DF_dict, PL_dict, "DF")

        #remove conflict between padloc and defencefinder and all three individual db to combined dict
        if accession in padloc_db_dict.keys():
            for system in PL_dict.keys():
                system_db = PL_dict[system]
                combined_dict = test_system_for_overlap(system_db, combined_dict, "PL", DF_dict, conflict_winner, accession, quiet, warning_file)

        if accession in df_db_dict.keys():
            for system in DF_dict.keys():
                system_db = DF_dict[system]
                combined_dict = test_system_for_overlap(system_db, combined_dict, "DF", PL_dict, conflict_winner, accession, quiet, warning_file)

        if accession in hmmer_db_dict.keys():
            for system_db in hmmer_db_dict[accession]:
                genes = system_db[2]
                if genes in combined_dict.keys():
                    if not quiet:
                        print ("Warning: customly searched for system overlaps perfect with system in DefenseFinder or PADLOC. custom system kept. Pertaining: sequence %s, conflicting_system %s, genes %s"%(accession, combined_dict[genes], genes))
                    with open (warning_file, "a+") as w_f:
                        w_f.write("Warning: customly searched for system overlaps perfect with system in DefenseFinder or PADLOC. custom system kept. Pertaining: sequence %s, conflicting_system %s, genes %s"%(accession, combined_dict[genes], genes))
                    combined_dict[genes]["AS"] = system_db
                else:
                    combined_dict[genes] = {"AS":system_db}

        #write output file
        with open(os.path.join(combined_output_dir,accession+".csv"), "w") as output_file:
            output_file.write("accession,system family DefenseFinder,system family PADLOC,system family custom,subtype DefenseFinder,subtype PADLOC,subtype custom,present in defence island?, genes,genes DefenseFinder,genes PADLOC, genes custom\n")           
            
            if high_confidence: #in case high confidence (i.e. found by both PADLOC and DF or custom system) is desired
                genes_list = list(combined_dict.keys())
                for genes in genes_list:
                    if "AS" in combined_dict[genes]:
                        pass
                    elif ("DF" not in combined_dict[genes]) or ("PL" not in combined_dict[genes]):
                        combined_dict.pop(genes)

            # Check whether proteins systems are present in a defence island
            defence_systems = list(combined_dict.keys())
            defence_island_presence = check_proximity(defence_systems)

            #add defence island information to the information dictionary
            for system_index in range(len(defence_systems)):
                defence_system = defence_systems[system_index]
                present_in_defence_island = defence_island_presence[system_index]
                combined_dict[defence_system]["present_in_defence_island"] = present_in_defence_island

                #removes proteins that are not found in a defence island (i.e. within 20 proteins of another system
                skip = False
                if cut_no_island == True:
                    if not present_in_defence_island:
                        combined_dict.pop(defence_system)
                elif type(cut_no_island) == list:
                    if not present_in_defence_island:
                        if "PL" in combined_dict[defence_system]:
                            if combined_dict[defence_system]["PL"][0] in cut_no_island:
        ##                            print (DS_dict)
    ##                            print (combined_dict[defence_system]["PL"][0])  
                                combined_dict.pop(defence_system)
                                skip = True
                        if not skip:
                            if "DF" in combined_dict[defence_system]:
                                if combined_dict[defence_system]["DF"][0] in cut_no_island:
                                    print (combined_dict[defence_system]["DF"][0])  
                                    combined_dict.pop(defence_system)
                    
##                    print (system_fam, cut_no_island)
                        
            #combine data
            for genes in combined_dict.keys():
                info = combined_dict[genes]

                #make summary file for accession
                system_PL = "N.A."
                system_DF = "N.A."
                system_AS = "N.A."  
                subtype_PL = "N.A."
                subtype_DF = "N.A."
                subtype_AS = "N.A."
                genes_DF = "N.A."
                genes_PL = "N.A."
                genes_AS = "N.A."
                defence_island = str(info["present_in_defence_island"])

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
                        
                if "AS" in info.keys():
                    system_AS = info["AS"][0]
                    subtype_AS = info["AS"][1]
                    genes_AS = "N.A."
                    genes_AS = info["AS"][3].replace(",",";")
                    if system_AS not in defence_system_list:
                        defence_system_list.append(system_AS)

                output_file.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(accession,system_DF,system_PL,system_AS,subtype_DF,subtype_PL,subtype_AS,defence_island,genes,genes_DF,genes_PL,genes_AS))

                #make db for summary file
                if system_AS != "N.A.":
                    db_dictionary[accession].append(system_AS)
                elif system_PL == system_DF:
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
        elif testing == "DF":
            if genes in combined_dict.keys():
                combined_dict[genes]["DF"] = db
            else:
                combined_dict[genes] = {"DF":db}
    return combined_dict

# Function to parse protein systems and return a list of tuples (contig, position)
def parse_protein_systems(protein_system):
    proteins = protein_system.split(';')
    parsed_proteins = []
    for protein in proteins:
        contig = '_'.join(protein.split('_')[:-1])
        position = int(protein.split('_')[-1])
        parsed_proteins.append((contig, position))
    return parsed_proteins

# Function to check if another system is within 20 proteins on the same contig
def check_proximity(protein_systems):
    parsed_systems = [parse_protein_systems(system) for system in protein_systems]
    
    results = []
    
    for i, system in enumerate(parsed_systems):
        contigs_positions = {contig: [] for contig, _ in system}
        for contig, position in system:
            contigs_positions[contig].append(position)
        
        is_within_20 = False
        
        for j, other_system in enumerate(parsed_systems):
            if i != j:
                for contig, position in other_system:
                    if contig in contigs_positions:
                        if any(abs(position - p) <= 20 for p in contigs_positions[contig]):
                            is_within_20 = True
                            break
                if is_within_20:
                    break
        
        results.append(is_within_20)
    
    return results


def write_summary_file(output_dir,db_dictionary,defence_system_list,output_file):
    #export relevant data to csv
    db_file = os.path.join(output_dir,output_file)
    with open(db_file, "w") as output_file:
        defence_systems = ",".join(defence_system_list)
        header = "accession,%s,\n"%(defence_systems)
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
