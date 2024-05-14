import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from intervaltree import IntervalTree

def process_bed_files(bed_folder):
    dicionario = {}
    lista_nomes = ["None"]

    for bed_file in os.listdir(bed_folder):
        if bed_file.endswith(".bed"):
            file_path = os.path.join(bed_folder, bed_file)
            fasta = pd.read_table(file_path, header=None, sep="\t")

            for i in range(len(fasta)):
                tamanho = int(fasta[3][i].split(";")[2].replace("ORF_len=", ""))
                contig = str(fasta[0][i])
                inicio = int(fasta[1][i])
                final = int(fasta[2][i])

                if contig not in lista_nomes:
                    if tamanho > 150:
                        dicionario[contig] = IntervalTree()
                        dicionario[contig][inicio:final+1] = (inicio, final)
                        lista_nomes.append(contig)
                else:
                    if tamanho > 150:
                        dicionario[contig][inicio:final+1] = (inicio, final)

    return dicionario

def process_fasta_files(fasta_folder, dicionario, output_folder):
    for fasta_file in os.listdir(fasta_folder):
        if fasta_file.endswith(".fasta"):
            file_path = os.path.join(fasta_folder, fasta_file)
            output_file = os.path.join(output_folder, fasta_file)
            with open(output_file, "w") as fasta_arch:
                for record in SeqIO.parse(file_path, "fasta"):
                    contig = str(record.id)
                    if contig in dicionario:
                        numerador = 0
                        for interval in sorted(dicionario[contig]):
                            inicio, final = interval.begin, interval.end - 1
                            nucleotide_seq = record.seq[inicio:final+1] if final + 1 == len(record.seq) else record.seq[inicio:final+3]
                            if len(nucleotide_seq) % 3 != 0:
                                nucleotide_seq = nucleotide_seq[:-(len(nucleotide_seq) % 3)]  # Truncar para múltiplo de 3
                            protein_seq = Seq(nucleotide_seq).translate(to_stop=True)
                            fasta_arch.write(f">{contig}.{numerador}\n{protein_seq}\n")
                            numerador += 1

def main(bed_folder, fasta_folder, output_folder):
    dicionario = process_bed_files(bed_folder)
    process_fasta_files(fasta_folder, dicionario, output_folder)

if __name__ == "__main__":
    bed_folder = "ViewVir-results/ORFs/"  # Substitua pelo caminho real para a pasta BED
    fasta_folder = "ViewVir-results/"  # Substitua pelo caminho real para a pasta FASTA
    output_folder = "ViewVir-results/ORFs/"  # Substitua pelo caminho real para a pasta de saída
    os.makedirs(output_folder, exist_ok=True)
    main(bed_folder, fasta_folder, output_folder)
