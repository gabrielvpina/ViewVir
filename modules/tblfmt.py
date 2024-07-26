import os
import pandas as pd
from Bio import SeqIO

ncbiSpecie = pd.read_csv("data/NCBI-species_name.csv", names=["Species","Family","Genome.composition"])
ncbiNames = pd.read_csv("data/NCBI-organism_name.csv", names=["Species","Family","Genome.composition"])
allVirus = pd.concat([ncbiNames, ncbiSpecie])

def viralFilter(vvfolder):

    # Pasta para os resultados
    viewvirFolder = vvfolder
    if not os.path.exists(viewvirFolder):
        os.makedirs(viewvirFolder)

    for file_path in os.listdir(vvfolder): 
        if file_path.endswith("_proc.tsv"):
            # Processando a Tabela
            # Leia o arquivo
            file = pd.read_csv(os.path.join(vvfolder,file_path), sep="\t",on_bad_lines='skip')

            # operação de junção
            file = pd.merge(file, allVirus, on="Species", how="left")

            # ordem das colunas
            outrasCols = list(set(file.columns) - {"Genome.composition"})
            file = file[["QuerySeq", "SubjectSeq", "QseqLength", "SseqLength", "Pident", "Evalue", "QCover",
                      "SubjTitle", "Species", "Family", "Genome.composition", "FullQueryLength"]]

            # Filtro de sequencias
            file = file[file["QseqLength"] >= 500]

            # Remover duplicatas
            file = file.drop_duplicates()
            # Obter o nome base do arquivo de entrada
            inputBasename = os.path.basename(file_path)
            # Caminho do arquivo de saída
            output_file = os.path.join(viewvirFolder, inputBasename.replace("_proc.tsv", "_processed.tsv"))
    
            # Escreva os dados em um arquivo de saída dentro do diretório criado
            file.to_csv(output_file, sep="\t", index=False)


            # Filtro para vírus de RNA
            tableRNA = file[file["Genome.composition"].isin(["ssRNA(-)", "ssRNA(+/-)", "ssRNA(+)", "ssRNA-RT", "RNA", "unknown", "NA"])]

            # Dados em um arquivo de saída separado para vírus não DNA
            nonDNA_table = os.path.join(viewvirFolder, inputBasename.replace("_proc.tsv", "_nonDNA.tsv"))
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
            os.remove(os.path.join(vvfolder,file_path))



def renameFasta(outputFolder):

    for fasta in os.listdir(outputFolder):
        if fasta.endswith(".fasta"):
            input_fasta = os.path.join(outputFolder, fasta)
    # Lista para armazenar os registros modificados
    registros_modificados = []

    # Iterar sobre os registros no arquivo FASTA de entrada
    for i, record in enumerate(SeqIO.parse(input_fasta, "fasta"), start=1):
        # Modificar o cabeçalho
        novo_id = f"contig_{i:02d}"
        record.id = novo_id
        record.description = novo_id  # Atualizar a descrição também

        # Adicionar o registro modificado à lista
        registros_modificados.append(record)

    fastafile = os.path.join
    # Escrever os registros modificados no novo arquivo FASTA
    with open(input_fasta, "w") as output_handle:
        SeqIO.write(registros_modificados, output_handle, "fasta")

        
