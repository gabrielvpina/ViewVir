import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from intervaltree import IntervalTree


def bedgc1(outputFolder):
    # Inicialização dos dicionários e listas
    dicioGC1 = {}
    listGC1 = []

    # Iterar sobre os arquivos BED
    bed_folder = outputFolder
    for bedfile in os.listdir(bed_folder):
        if bedfile.endswith("_ORFgc1.bed"):
            outputFolder_path = os.path.join(bed_folder, bedfile)
            fasta_bed = pd.read_table(outputFolder_path, header=None, sep="\t")

            # Iterar sobre as linhas do arquivo BED
            for i in range(len(fasta_bed)):
                frame = int(fasta_bed.iloc[i, 3].split(";")[3].replace("ORF_frame=", ""))
                tamanho = int(fasta_bed.iloc[i, 3].split(";")[2].replace("ORF_len=", ""))
                contig = str(fasta_bed.iloc[i, 0])
                inicio = int(fasta_bed.iloc[i, 1])
                final = int(fasta_bed.iloc[i, 2])

                # Verificação e armazenamento do maior intervalo para cada contig
                if contig not in listGC1:
                    if tamanho >= 300:
                        dicioGC1[contig] = IntervalTree()
                        dicioGC1[contig][inicio:final] = (inicio, final, frame)
                        listGC1.append(contig)
                        print(f"Added new contig {contig} with interval {inicio}-{final}")
                else:
                    if tamanho >= 300:
                        current_intervals = list(dicioGC1[contig])
                        if not current_intervals:
                            dicioGC1[contig][inicio:final] = (inicio, final, frame)
                            print(f"Added interval for {contig}: {inicio}-{final}")
                        else:
                            max_interval = max(current_intervals, key=lambda x: x.end - x.begin)
                            if (final - inicio) > (max_interval.end - max_interval.begin):
                                dicioGC1[contig].remove(max_interval)
                                dicioGC1[contig][inicio:final] = (inicio, final, frame)
                                print(f"Replaced interval for {contig}: {inicio}-{final} replacing {max_interval}")

    # Criação do arquivo de saída para as ORFs
    output_fasta = os.path.join(outputFolder, "ORFgc1.fasta")
    with open(output_fasta, "w") as fasta_arch:
        # Processamento dos arquivos FASTA
        fasta_files = [os.path.join(outputFolder, fasta) for fasta in os.listdir(outputFolder) if fasta.endswith(".fasta")]
        for fasta_file in fasta_files:
            print(f"Processing {fasta_file}")
            # Iterar sobre os registros FASTA no arquivo
            for record in SeqIO.parse(fasta_file, "fasta"):
                contig = str(record.id)
                # Verificar se o contig está presente no dicionário
                if contig in dicioGC1:
                    numerador = 0
                    # Iterar sobre os intervalos para o contig específico
                    for interval in sorted(dicioGC1[contig]):
                        inicio, final = interval.begin, interval.end
                        # Extrair a sequência do registro FASTA com base nos intervalos
                        nucleotide_seq = record.seq[inicio:final + 1] #record.seq[inicio:final + 1] if final + 1 == len(record.seq) else record.seq[inicio:final + 3]
                        if len(nucleotide_seq) % 3 != 0:
                            nucleotide_seq = nucleotide_seq[:-(len(nucleotide_seq) % 3)]  # Truncar para múltiplo de 3
                            protein_seq = Seq(nucleotide_seq).translate()
                        # Escrever a sequência no arquivo de saída
                        fasta_arch.write(f">{contig}.{numerador}\n{nucleotide_seq}\n")
                        numerador += 1
                        print(f"Written ORF for {contig}.{numerador}")
