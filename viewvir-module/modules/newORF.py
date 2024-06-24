import os
from Bio import SeqIO
from collections import defaultdict

def gc1_ORFs(outputFolder):
    # Dicionário para armazenar as ORFs por contig
    orfs_dict = defaultdict(list)

    # Ler arquivos
    for fasta in os.listdir(outputFolder):
        if fasta.endswith("_ORFgc1.fasta"):
            input_fasta = os.path.join(outputFolder, fasta)

    # Ler o arquivo FASTA
    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.description
        contig = header.split(".")[0]  # Pega o contig (primeira parte da descrição)
        orfs_dict[contig].append(record)

    # Lista para armazenar os registros selecionados
    maiores_orfs = []

    # Processar cada contig
    for contig, orfs in orfs_dict.items():
        # Ordenar as ORFs por tamanho (comprimento da sequência)
        orfs.sort(key=lambda x: len(x.seq), reverse=True)
        # Selecionar as duas maiores ORFs
        maiores_orfs.extend(orfs[:2])

    # Escrever as ORFs selecionadas em um novo arquivo FASTA
    with open(input_fasta, "w") as output_handle:
        SeqIO.write(maiores_orfs, output_handle, "fasta")


def gc5_ORFs(outputFolder):
    # Dicionário para armazenar as ORFs por contig
    orfs_dict = defaultdict(list)

    # Ler arquivos
    for fasta in os.listdir(outputFolder):
        if fasta.endswith("_ORFgc5.fasta"):
            input_fasta = os.path.join(outputFolder, fasta)

    # Ler o arquivo FASTA
    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.description
        contig = header.split(".")[0]  # Pega o contig (primeira parte da descrição)
        orfs_dict[contig].append(record)

    # Lista para armazenar os registros selecionados
    maiores_orfs = []

    # Processar cada contig
    for contig, orfs in orfs_dict.items():
        # Ordenar as ORFs por tamanho (comprimento da sequência)
        orfs.sort(key=lambda x: len(x.seq), reverse=True)
        # Selecionar as duas maiores ORFs
        maiores_orfs.extend(orfs[:2])

    # Escrever as ORFs selecionadas em um novo arquivo FASTA
    with open(input_fasta, "w") as output_handle:
        SeqIO.write(maiores_orfs, output_handle, "fasta")



def gc11_ORFs(outputFolder):
    # Dicionário para armazenar as ORFs por contig
    orfs_dict = defaultdict(list)

    # Ler arquivos
    for fasta in os.listdir(outputFolder):
        if fasta.endswith("_ORFgc11.fasta"):
            input_fasta = os.path.join(outputFolder, fasta)

    # Ler o arquivo FASTA
    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.description
        contig = header.split(".")[0]  # Pega o contig (primeira parte da descrição)
        orfs_dict[contig].append(record)

    # Lista para armazenar os registros selecionados
    maiores_orfs = []

    # Processar cada contig
    for contig, orfs in orfs_dict.items():
        # Ordenar as ORFs por tamanho (comprimento da sequência)
        orfs.sort(key=lambda x: len(x.seq), reverse=True)
        # Selecionar as duas maiores ORFs
        maiores_orfs.extend(orfs[:2])

    # Escrever as ORFs selecionadas em um novo arquivo FASTA
    with open(input_fasta, "w") as output_handle:
        SeqIO.write(maiores_orfs, output_handle, "fasta")

    



