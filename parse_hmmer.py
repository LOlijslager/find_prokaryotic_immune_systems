import csv
import os
from collections import defaultdict
from Bio import SearchIO
from Bio import SeqIO

def parse_custom_systems(accessions,custom_system_dir,hmm_file_dir):
    #final_systems_file = os.path.join(HMMER_systems_dir,base_name+"_HMMER_systems.tsv")
    hmmer_db_dict = {}

    systems_dict = parse_custom_systems_metadata(custom_system_dir)

    for accession in accessions:
        for system in systems_dict.keys():
            metadatas = check_system_presence(systems_dict,system,accession,hmm_file_dir)
            if metadatas != []:
                for metadata in metadatas:
                    if accession in hmmer_db_dict.keys():
                        hmmer_db_dict[accession].append(metadata)
                    else:
                        hmmer_db_dict[accession] = [metadata]
    return hmmer_db_dict


def check_system_presence(systems_dict, system, accession, hmm_file_dir):
    gene_list = []
    gene_names = []
    gene_positions = []
    metadatas = []

    for i in range(len(systems_dict[system]["HMM-profile"])):
        hmm_profile = systems_dict[system]["HMM-profile"][i]
        gene_name = hmm_profile.replace(".hmm", "")
        subtype = systems_dict[system]["Subtype"][i]
        bitscore = int(systems_dict[system]["Bitscore"][i])
        required = systems_dict[system]["Required"][i]

        hmm_name = os.path.splitext(hmm_profile)[0]
        output_file_path = os.path.join(hmm_file_dir, accession + "_" + hmm_name + ".out")
        gene_present, gene = parse_hmmer(output_file_path, bitscore)

        if gene_present:
            gene_list.append(gene)
            gene_names.append(gene_name)
            gene_position = int(gene.split('_')[-1])  # Extract the position number
            gene_positions.append((gene_position, required, gene, gene_name, subtype))

    # Sort genes by position
    gene_positions.sort()

    # Group genes into clusters within 20 positions of each other
    current_cluster = []
    for i, (position, required, gene, gene_name, subtype) in enumerate(gene_positions):
        if not current_cluster or position - current_cluster[-1][0] <= 20:
            current_cluster.append((position, required, gene, gene_name, subtype))
        else:
            # Process the current cluster
            if all_required_genes_present(current_cluster, systems_dict[system]["HMM-profile"], systems_dict[system]["Required"]):
                cluster_genes = [g[2] for g in current_cluster]
                cluster_gene_names = [g[3] for g in current_cluster]
                cluster_subtype = current_cluster[0][4]
                metadata = [system, cluster_subtype, ";".join(cluster_genes), ",".join(cluster_gene_names)]
                metadatas.append(metadata)
            current_cluster = [(position, required, gene, gene_name, subtype)]
            
    # Process the last cluster
    if current_cluster and all_required_genes_present(current_cluster, systems_dict[system]["HMM-profile"], systems_dict[system]["Required"]):
        cluster_genes = [g[2] for g in current_cluster]
        cluster_gene_names = [g[3] for g in current_cluster]
        cluster_subtype = current_cluster[0][4]
        metadata = [system, cluster_subtype, ";".join(cluster_genes), ",".join(cluster_gene_names)]
        metadatas.append(metadata)
    return metadatas

def all_required_genes_present(cluster, hmm_profiles, required_flags):
    required_gene_names = []
    for hmm_profile, required in zip(hmm_profiles, required_flags):
        if required == "FALSE":
            required = False
        else:
            required = True
        if required:
            gene_name = hmm_profile.replace(".hmm", "")
            required_gene_names.append(gene_name)


    cluster_gene_names = []
    for _, _, _, gene_name, _ in cluster:
        cluster_gene_names.append(gene_name)

    all_genes_present = True
    for gene in required_gene_names:
        if gene not in cluster_gene_names:
            all_genes_present = False

    return all_genes_present


def parse_hmmer(output_file_path,bitscore):
    gene_present = False
    gene = ""
    gene_name = ""
    for qresult in SearchIO.parse(output_file_path, "hmmsearch3-domtab"):
        match = ""
        method = ""
        for hit in qresult:
            if hit.bitscore > bitscore:
                gene_present = True
                gene = hit.id
##                for HSP in hit:
##                    for HSPFragment in HSP:
####                        print (proteome_name,hmm_ dir,hit.id,hit.bitscore,hit.evalue,HSP.acc_avg)
##                        for seq_record in SeqIO.parse(proteome_file_path, "fasta"):
##                            if seq_record.id == hit.id:
##                                with open("HMMER/protein_hits/"+proteome_name+"_"+hmm_dir+"_"+str(n)+".fasta", "w") as output_handle:
##                                    SeqIO.write(seq_record, output_handle, "fasta")
##                                n += 1
    return gene_present,gene

def parse_custom_systems_metadata(custom_system_dir):
    #Make a dictionary out of the custom immune systems
    custom_system_file_path = os.path.join(custom_system_dir, "systems.csv")
    systems_dict = defaultdict(lambda: defaultdict(list))

    with open(custom_system_file_path, mode='r', newline='', encoding='utf-8-sig') as file:
        reader = csv.DictReader(file)
        for row in reader:
            system = row.pop(next(iter(row)))  # Get the first column value as the key
            for key, value in row.items():
                systems_dict[system][key].append(value)
    systems_dict = {k: dict(v) for k, v in systems_dict.items()}

    return systems_dict


