import os
import pandas as pd
from Bio import SeqIO


# Listar arquivos .bed
def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]

#arquivo = list_full_paths("ViewVir-results/ORFs/")

arquivo = input(" -> Arquivo BED do orfipy:\n")
arquivo = arquivo.replace("\\", "/")
arquivo = arquivo.replace("\"", "")
fasta = pd.read_table(arquivo, header=None, sep="\t").sort_values(by=0)



def loop_interno(contig):
    global dicionario
    global inicio
    global final

    atualizador = 0
    sobrepor = 0

    for j in range(len(dicionario[contig])):
        tamanho_grava = dicionario[contig][j]
        tamanho_grava_ini = tamanho_grava[0]
        tamanho_grava_fim = tamanho_grava[1]

        if len(set(list(range(inicio, final + 1))) & set(list(range(tamanho_grava_ini, tamanho_grava_fim)))) < 1:
            for k,q in dicionario[contig]:
                if len(set(list(range(inicio, final + 1))) & set(list(range(k, q)))) > 0:
                    sobrepor += 1
            # adiciona sequências novas que não são sobrepostas
            if list([inicio, final]) not in dicionario[contig]:
                if sobrepor == 0:
                    dicionario[contig] += [list([inicio, final])]
                    atualizador = 1
                    break
        else:
            if len(set(list(range(inicio, final + 1)))) > len(set(list(range(tamanho_grava_ini, tamanho_grava_fim)))):
                if list([inicio, final]) not in dicionario[contig]:
                    dicionario[contig][j] = list([inicio, final])
                else:
                    if [inicio, final] != [tamanho_grava_ini, tamanho_grava_fim]:
                        atualizador = 1
                        del dicionario[contig][j]
                        break
            else:
                break
    if atualizador == 1:
        loop_interno(contig)


dicionario = {}
lista_nomes = ["None"]

for i in range(len(fasta)):
    tamanho = str(fasta[3][i]).split(";")[2].replace("ORF_len=","")
    contig = str(fasta[0][i])
    inicio = int(fasta[1][i])
    final = int(fasta[2][i])

    # caso não esteja no dicionário
    if contig not in str(lista_nomes):
        if int(tamanho) > 150:
            dicionario[contig] = [list([inicio,final])]
            lista_nomes.append(contig)
    #caso esteja no dicionário
    else:
        if int(tamanho) > 150:
            loop_interno(contig)
for i in dicionario:
    print(str(i) + str(dicionario[i]))
arquivo_1 = input(" -> Arquivo com fasta:\n")
arquivo_1 = arquivo_1.replace("\\", "/")
arquivo_1 = arquivo_1.replace("\"", "")
fasta_arch = open("orfs_atualizada_virus_ferment_str11.fasta","w")
for i in SeqIO.parse(arquivo_1, "fasta"):
    numerador = 0
    for j in range(len(dicionario[str(i.id)])):
        if len(i.seq) == int(dicionario[str(i.id)][j][1]+1):
            fasta_arch.write(">"+str(i.id)+"&numerador"+str(numerador)+"\n"+str(i.seq)[int(dicionario[str(i.id)][j][0]):int(dicionario[str(i.id)][j][1])]+"\n")
            numerador += 1
        else:
            fasta_arch.write(">"+str(i.id)+"&numerador"+str(numerador)+"\n"+str(i.seq)[int(dicionario[str(i.id)][j][0]):int(dicionario[str(i.id)][j][1])]+"\n")
            numerador += 1