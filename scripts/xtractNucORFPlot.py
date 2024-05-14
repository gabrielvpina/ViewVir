import pandas as pd
from Bio import SeqIO, Seq
from intervaltree import IntervalTree
import os
from dna_features_viewer import GraphicFeature, GraphicRecord

def process_bed_files(bed_folder):
    dicionario = {}
    lista_nomes = []

    for bed_file in os.listdir(bed_folder):
        if bed_file.endswith(".bed"):
            file_path = os.path.join(bed_folder, bed_file)
            fasta = pd.read_table(file_path, header=None, sep="\t")

            for i in range(len(fasta)):
                tamanho = int(fasta.iloc[i, 3].split(";")[2].replace("ORF_len=", ""))
                contig = str(fasta.iloc[i, 0])
                inicio = int(fasta.iloc[i, 1])
                final = int(fasta.iloc[i, 2])

                if contig not in lista_nomes:
                    if tamanho > 150:
                        dicionario[contig] = IntervalTree()
                        dicionario[contig][inicio:final+1] = (inicio, final)
                        lista_nomes.append(contig)
                else:
                    if tamanho > 150:
                        contida = False
                        for interval in dicionario[contig]:
                            # Verifica se a nova sequência está contida em uma já existente
                            if inicio >= interval.begin and final <= interval.end:
                                contida = True
                                break
                        if not contida:
                            # Remove intervalos menores que estão sobrepostos com a nova sequência
                            dicionario[contig].chop(inicio, final+1)
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
                            fasta_arch.write(f">{contig}.{numerador}\n{nucleotide_seq}\n")
                            numerador += 1


def main(bed_folder, fasta_folder, output_folder):
    dicionario = process_bed_files(bed_folder)
    process_fasta_files(fasta_folder, dicionario, output_folder)

    for fasta_file in os.listdir(fasta_folder):
        if fasta_file.endswith(".fasta"):
            file_path = os.path.join(fasta_folder, fasta_file)
            for record in SeqIO.parse(file_path, "fasta"):
                contig = str(record.id)
                if contig in dicionario:
                    length = len(record.seq)
                    features = []

                    for interval in sorted(dicionario[contig]):
                        inicio, final = interval.begin, interval.end - 1
                        features.append(GraphicFeature(start=inicio, end=final, strand=+1, color="#ffd700"))

                    record_graphic = GraphicRecord(sequence_length=length, features=features)
                    ax, _ = record_graphic.plot(figure_width=5)
                    ax.figure.savefig(os.path.join(output_folder, f'sequence_{contig}.pdf'), bbox_inches='tight')

if __name__ == "__main__":
    bed_folder = "ViewVir-results/ORFs/"  # Substitua pelo caminho real para a pasta BED
    fasta_folder = "ViewVir-results/"  # Substitua pelo caminho real para a pasta FASTA
    output_folder = "ViewVir-results/Gene-View/"  # Substitua pelo caminho real para a pasta de saída
    os.makedirs(output_folder, exist_ok=True)
    main(bed_folder, fasta_folder, output_folder)
