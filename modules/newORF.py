import os
from Bio import SeqIO
from collections import defaultdict

def gc1_ORFs(outputFolder,nORF):
    
    orfs_dict = defaultdict(list)

    # read files
    for fasta in os.listdir(outputFolder):
        if fasta.endswith("_ORFgc1.fasta"):
            input_fasta = os.path.join(outputFolder, fasta)


    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.description
        contig = header.split(".")[0]  # get the first part of description
        orfs_dict[contig].append(record)

    
    maiores_orfs = []

    # process each contig
    for contig, orfs in orfs_dict.items():
        # sorting
        orfs.sort(key=lambda x: len(x.seq), reverse=True)
        # select biggest orfs
        maiores_orfs.extend(orfs[:int(nORF)])

    # write new fasta w/ biggest orfs
    with open(input_fasta, "w") as output_handle:
        SeqIO.write(maiores_orfs, output_handle, "fasta")


def gc5_ORFs(outputFolder,nORF):
    
    orfs_dict = defaultdict(list)

    
    for fasta in os.listdir(outputFolder):
        if fasta.endswith("_ORFgc5.fasta"):
            input_fasta = os.path.join(outputFolder, fasta)

    
    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.description
        contig = header.split(".")[0]  
        orfs_dict[contig].append(record)

    
    maiores_orfs = []

    
    for contig, orfs in orfs_dict.items():
        
        orfs.sort(key=lambda x: len(x.seq), reverse=True)
        # Selecionar as maiores ORFs
        maiores_orfs.extend(orfs[:int(nORF)])

    
    with open(input_fasta, "w") as output_handle:
        SeqIO.write(maiores_orfs, output_handle, "fasta")



def gc11_ORFs(outputFolder,nORF):
    
    orfs_dict = defaultdict(list)

    
    for fasta in os.listdir(outputFolder):
        if fasta.endswith("_ORFgc11.fasta"):
            input_fasta = os.path.join(outputFolder, fasta)

    
    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.description
        contig = header.split(".")[0] 
        orfs_dict[contig].append(record)

    
    maiores_orfs = []

    
    for contig, orfs in orfs_dict.items():
        
        orfs.sort(key=lambda x: len(x.seq), reverse=True)
        
        maiores_orfs.extend(orfs[:int(nORF)])

    
    with open(input_fasta, "w") as output_handle:
        SeqIO.write(maiores_orfs, output_handle, "fasta")

    



