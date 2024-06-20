import os
import pandas as pd

ncbiSpecie = pd.read_csv("data/NCBI_virSpecies.csv", names=["Species", "Genome.composition"])
ncbiNames = pd.read_csv("data/NCBI_virName.csv", names=["Species", "Genome.composition"])
allVirus = pd.concat([ncbiNames, ncbiSpecie])

def process_diamondTbl(file_path,vvfolder):

    # Pasta para os resultados
    viewvirFolder = "ViewVir-results"
    if not os.path.exists(viewvirFolder):
        os.makedirs(viewvirFolder)

    # Processando a Tabela
    # Leia o arquivo
    file = pd.read_csv(file_path, sep="\t")

    # operação de junção
    file = pd.merge(file, allVirus, on="Species", how="left")

    # ordem das colunas
    outrasCols = list(set(file.columns) - {"Genome.composition"})
    file = file[["QuerySeq", "SubjectSeq", "QseqLength", "SseqLength", "Pident", "Evalue", "QCover",
                      "SubjTitle", "Species", "Genome.composition", "FullQueryLength"]]

    # Filtro de sequencias
    file = file[file["QseqLength"] >= 500]

    # Remover duplicatas
    file = file.drop_duplicates()
    # Obter o nome base do arquivo de entrada
    inputBasename = os.path.basename(file_path)
    # Caminho do arquivo de saída
    output_file = os.path.join(viewvirFolder, inputBasename.replace(".tsv", "_processed.tsv"))
    
    # Escreva os dados em um arquivo de saída dentro do diretório criado
    file.to_csv(output_file, sep="\t", index=False)


    # Filtro para vírus de RNA
    tableRNA = file[file["Genome.composition"].isin(["ssRNA(-)", "ssRNA(+/-)", "ssRNA(+)", "ssRNA-RT", "RNA", "unknown", "NA"])]

    # Dados em um arquivo de saída separado para vírus não DNA
    nonDNA_table = os.path.join(viewvirFolder, inputBasename.replace(".tsv", "_nonDNA.tsv"))
    tableRNA.to_csv(nonDNA_table, sep="\t", index=False)


    #### Generate fasta file
    #vvfolder = "ViewVir-results"
    inputBasename = os.path.basename(vvfolder)
    files = os.listdir(vvfolder)

    output_pre = os.path.join(vvfolder, inputBasename.replace("_nonDNA.tsv", "_pre.tsv"))
    #tableRNA = pd.read_csv(vvfolder + nonDNA)
    arquivo_PRE = tableRNA[["QuerySeq", "FullQueryLength"]]
    arquivo_PRE.to_csv(output_pre, sep="\t", index=False)

    # Processar _pre.tsv para gerar RNA-virus.fasta
    samp = os.path.basename(output_pre).replace("_pre.tsv", "")
    with open(output_pre, "r") as infile, open(f"{vvfolder}/{samp}_nonDNA.fasta", "w") as outfile:
        lines = infile.readlines()
        for line in lines[1:]:  # Ignorar a primeira linha
            line = line.strip().replace('\t', '\n')
            outfile.write(f">{line}\n")
    
    os.remove(output_pre)