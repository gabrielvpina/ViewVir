import os
import pandas as pd
import plotly.express as px
from plotly.offline import plot
import subprocess
from dna_features_viewer import GraphicFeature, GraphicRecord
from intervaltree import IntervalTree
from Bio import SeqIO, Seq



##########################################################################################

# Leitura dos dados de vírus
ictvMSL = pd.read_csv("data/ICTV_MSL.csv", names=["Species", "Genome.composition"])
ictvVMR_species = pd.read_csv("data/ICTV_VMR-species.csv")
ictvVMR_virNames = pd.read_csv("data/ICTV_VMR-virNames.csv")
ncbiSpecie = pd.read_csv("data/NCBI_virSpecies.csv", names=["Species", "Genome.composition"])
ncbiNames = pd.read_csv("data/NCBI_virName.csv", names=["Species", "Genome.composition"])

# Merge de todas as tabelas
allVirus = pd.concat([ncbiNames, ncbiSpecie])

# Listar arquivos .tsv
def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]

arquivos = list_full_paths("diamond-processed/")

# Crie a pasta para os resultados
viewvirFolder = "ViewVir-results"
if not os.path.exists(viewvirFolder):
    os.makedirs(viewvirFolder)

############################################################################################

# Main loop

for arquivo_path in arquivos:
    # Processando a Tabela
    # Leia o arquivo
    arquivo = pd.read_csv(arquivo_path, sep="\t")

    # operação de junção
    arquivo = pd.merge(arquivo, allVirus, on="Species", how="left")

    # ordem das colunas
    outrasCols = list(set(arquivo.columns) - {"Genome.composition"})
    arquivo = arquivo[["QuerySeq", "SubjectSeq", "QseqLength", "SseqLength", "Pident", "Evalue", "QCover",
                      "SubjTitle", "Species", "Genome.composition", "FullQueryLength"]]

    # Filtro de sequencias
    arquivo = arquivo[arquivo["QseqLength"] >= 500]

    # Remover duplicatas
    arquivo = arquivo.drop_duplicates()

    # nome do arquivo de entrada
    inputBasename = os.path.basename(arquivo_path)

    # Subdiretório com o nome do arquivo dentro de "ViewVir-results"
    outputFolder = os.path.join(viewvirFolder, inputBasename.replace(".fasta_merged.fasta.tsv", ""))
    os.makedirs(outputFolder, exist_ok=True)

    # Caminho completo do arquivo de saída dentro do diretório criado
    dmndTable_outPath = os.path.join(outputFolder, inputBasename.replace(".fasta_merged.fasta.tsv", "_diamond.tsv"))

    scatterPlt_outPath = os.path.join(outputFolder, inputBasename.replace(".fasta_merged.fasta.tsv", "_ScatterPlt.html"))

    # Escreva os dados em um arquivo de saída dentro do diretório criado
    arquivo.to_csv(dmndTable_outPath, sep="\t", index=False)

    # Filtro para vírus de RNA
    tableRNA = arquivo[arquivo["Genome.composition"].isin(["ssRNA(-)", "ssRNA(+/-)", "ssRNA(+)", "ssRNA-RT", "RNA", "unknown", "NA"])]
    tableRNA = tableRNA.drop_duplicates()

    # Dados em um arquivo de saída separado para vírus não DNA
    nonDNA_table = os.path.join(outputFolder, inputBasename.replace(".fasta_merged.fasta.tsv", "_nonDNA.tsv"))
    tableRNA.to_csv(nonDNA_table, sep="\t", index=False)



    ##############################################################################################
    ############## GENERATING FASTA ##############################################################


    output_pre = os.path.join(outputFolder, inputBasename.replace(".fasta_merged.fasta.tsv", "_pre_fasta.tsv"))
    arquivo_PRE = tableRNA[["QuerySeq", "FullQueryLength"]]
    arquivo_PRE.to_csv(output_pre, sep="\t", index=False)

    # Processar _pre_fasta.tsv para gerar RNA-virus.fasta
    samp = os.path.basename(output_pre).replace("_pre_fasta.tsv", "")
    with open(output_pre, "r") as infile, open(f"{outputFolder}/{samp}_nonDNA.fasta", "w") as outfile:
        lines = infile.readlines()
        for line in lines[1:]:  # Ignorar a primeira linha
            line = line.strip().replace('\t', '\n')
            outfile.write(f">{line}\n")
    
    os.remove(output_pre)

    ###############################################################################################
    ############################## SCATTER PLOT    ################################################
     
    
    

    # Gráfico Interativo
    arquivo["MatchSequence"] = arquivo["Species"] + " -> " + arquivo["SubjTitle"]
    arquivo["Genome.composition"] = arquivo["Genome.composition"].fillna("NA")

    fig = px.scatter(arquivo, x="QseqLength", y="MatchSequence", size="QCover", color="Genome.composition",
                     hover_data=["Pident", "Evalue"])
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(yaxis={'categoryorder': 'total ascending'}, title='ViewVir interactive scatter plot')
    plot(fig, filename=scatterPlt_outPath)





    ##################################################################################################
    ################################## ORFIPY  #######################################################
    
    lista_arquivos = [os.path.join(outputFolder, file)
                      for file in os.listdir(outputFolder)
                      if file.endswith(".fasta")]

    for arquivo in lista_arquivos:
        entrada = arquivo
        nome = os.path.basename(arquivo).replace("_nonDNA.fasta", "")
        saida = outputFolder

        # Executa orfipy
        orfipy_gc1 = [
            "orfipy", "--partial-3", "--partial-5", "--outdir", saida, entrada, "--bed", f"{nome}_ORFsGC1.bed"
        ]
        subprocess.run(orfipy_gc1, check=True)

        # Executa orfipy
        orfipy_gc5 = [
            "orfipy", "--table", "5", "--partial-3", "--partial-5", "--outdir", saida, entrada, "--bed",  f"{nome}_ORFsGC5.bed"
        ]
        subprocess.run(orfipy_gc5, check=True)

    # Remove arquivos .log
    log_files = [os.path.join(outputFolder, file) for file in os.listdir(outputFolder) if file.endswith(".log")]
    for log_file in log_files:
        os.remove(log_file)


    #############################################################################################
    ######################### SUBFILE ##############################################

    bed_files = [os.path.join(outputFolder, bed) for bed in os.listdir(outputFolder) if bed.endswith("ORFsGC1.bed")]
    for bed_file in bed_files:
        # Pasta de saída na mesma pasta que o arquivo FASTA
        output_subfolder = os.path.join(outputFolder, os.path.splitext(os.path.basename(bed_file))[0])  
        os.makedirs(output_subfolder, exist_ok=True)  # Criar a pasta de saída, se ainda não existir


     #############################################################################################
     ######################### BED GC1 FILE PROCESSING ##############################################

     # Processamento dos arquivos BED
        dicioGC1 = {}
        listGC1 = []

        # Iterar sobre os arquivos BED
        bed_folder = os.path.dirname(bed_file)
        for bed_file in os.listdir(bed_folder):
          if bed_file.endswith("ORFsGC1.bed"):
            bed_file_path = os.path.join(bed_folder, bed_file)
            fasta_bed = pd.read_table(bed_file_path, header=None, sep="\t")

            # Iterar sobre as linhas do arquivo BED
            for i in range(len(fasta_bed)):
                tamanho = int(fasta_bed.iloc[i, 3].split(";")[2].replace("ORF_len=", ""))
                contig = str(fasta_bed.iloc[i, 0])
                inicio = int(fasta_bed.iloc[i, 1])
                final = int(fasta_bed.iloc[i, 2])

                if contig not in listGC1:
                    if tamanho > 150:
                        dicioGC1[contig] = IntervalTree()
                        dicioGC1[contig][inicio:final+1] = (inicio, final)
                        listGC1.append(contig)
                else:
                    if tamanho > 150:
                        contida = False
                        for interval in dicioGC1[contig]:
                            # Verifica se a nova sequência está contida em uma já existente
                            if inicio >= interval.begin and final <= interval.end:
                                contida = True
                                break
                        if not contida:
                            # Remove intervalos menores que estão sobrepostos com a nova sequência
                            dicioGC1[contig].chop(inicio, final+1)
                            dicioGC1[contig][inicio:final+1] = (inicio, final)


    fasta_files = [os.path.join(outputFolder, fasta) for fasta in os.listdir(outputFolder) if fasta.endswith(".fasta")]
    for fasta_file in fasta_files:

    # Processamento dos arquivos FASTA
     output_subfolder_fasta = os.path.join(output_subfolder, "processed.fasta")  # Arquivo de saída de FASTA processado
     with open(output_subfolder_fasta, "w") as fasta_arch:
        # Iterar sobre os registros FASTA no arquivo
        for record in SeqIO.parse(fasta_file, "fasta"):
            contig = str(record.id)
            # Verificar se o contig está presente no dicionário
            if contig in dicioGC1:
                numerador = 0
                # Iterar sobre os intervalos para o contig específico
                for interval in sorted(dicioGC1[contig]):
                    inicio, final = interval.begin, interval.end - 1
                    # Extrair a sequência do registro FASTA com base nos intervalos
                    nucleotide_seq = record.seq[inicio:final+1] if final + 1 == len(record.seq) else record.seq[inicio:final+3]
                    if len(nucleotide_seq) % 3 != 0:
                        nucleotide_seq = nucleotide_seq[:-(len(nucleotide_seq) % 3)]  # Truncar para múltiplo de 3
                    # Escrever a sequência no arquivo de saída
                    fasta_arch.write(f">{contig}.{numerador}\n{nucleotide_seq}\n")
                    numerador += 1

    # Plotagem dos gráficos
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig = str(record.id)
        if contig in dicioGC1:
            length = len(record.seq)
            features = []

            for interval in sorted(dicioGC1[contig]):
                inicio, final = interval.begin, interval.end - 1
                features.append(GraphicFeature(start=inicio, end=final, strand=+1, color="#ffd700"))

            record_graphic = GraphicRecord(sequence_length=length, features=features)
            ax, _ = record_graphic.plot(figure_width=5)
            plot_output_path = os.path.join(output_subfolder, f'sequence_{contig}.pdf')  # Caminho de saída específico
            ax.figure.savefig(plot_output_path, bbox_inches='tight')





 ###########################################################################################################################


    #############################################################################################
    ######################### BED GC5 FILE PROCESSING ###########################################
    

    bedgc5_files = [os.path.join(outputFolder, bed) for bed in os.listdir(outputFolder) if bed.endswith("ORFsGC5.bed")]
    for bedgc5_file in bedgc5_files:
        out_subfoldergc5 = os.path.join(outputFolder, os.path.splitext(os.path.basename(bedgc5_file))[0])  
        os.makedirs(out_subfoldergc5, exist_ok=True) 
    # Processamento dos arquivos BED
        dicioGC5 = {}
        listGC5 = []

    # Iterar sobre os arquivos BED
    bedgc5_folder = os.path.dirname(bedgc5_file)
    for bedgc5_file in os.listdir(bedgc5_folder):
        if bedgc5_file.endswith("ORFsGC5.bed"):
            bedgc5_file_path = os.path.join(bedgc5_folder, bedgc5_file)
            fasta_bedgc5 = pd.read_table(bedgc5_file_path, header=None, sep="\t")

            # Iterar sobre as linhas do arquivo BED
            for i in range(len(fasta_bedgc5)):
                tamanho = int(fasta_bedgc5.iloc[i, 3].split(";")[2].replace("ORF_len=", ""))
                contig = str(fasta_bedgc5.iloc[i, 0])
                inicio = int(fasta_bedgc5.iloc[i, 1])
                final = int(fasta_bedgc5.iloc[i, 2])

                if contig not in listGC5:
                    if tamanho > 150:
                        dicioGC5[contig] = IntervalTree()
                        dicioGC5[contig][inicio:final+1] = (inicio, final)
                        listGC5.append(contig)
                else:
                    if tamanho > 150:
                        contida = False
                        for interval in dicioGC5[contig]:
                            # Verifica se a nova sequência está contida em uma já existente
                            if inicio >= interval.begin and final <= interval.end:
                                contida = True
                                break
                        if not contida:
                            # Remove intervalos menores que estão sobrepostos com a nova sequência
                            dicioGC5[contig].chop(inicio, final+1)
                            dicioGC5[contig][inicio:final+1] = (inicio, final)


    fastagc5_files = [os.path.join(outputFolder, fasta) for fasta in os.listdir(outputFolder) if fasta.endswith(".fasta")]
    for fastagc5 in fastagc5_files:

    # Processamento dos arquivos FASTA
     out_subfoldergc5_fasta = os.path.join(out_subfoldergc5, "processed.fasta")  # Arquivo de saída de FASTA processado
     with open(out_subfoldergc5_fasta, "w") as fasta_arch:
        # Iterar sobre os registros FASTA no arquivo
        for record in SeqIO.parse(fastagc5, "fasta"):
            contig = str(record.id)
            # Verificar se o contig está presente no dicionário
            if contig in dicioGC5:
                numerador = 0
                # Iterar sobre os intervalos para o contig específico
                for interval in sorted(dicioGC5[contig]):
                    inicio, final = interval.begin, interval.end - 1
                    # Extrair a sequência do registro FASTA com base nos intervalos
                    nucleotide_seq = record.seq[inicio:final+1] if final + 1 == len(record.seq) else record.seq[inicio:final+3]
                    if len(nucleotide_seq) % 3 != 0:
                        nucleotide_seq = nucleotide_seq[:-(len(nucleotide_seq) % 3)]  # Truncar para múltiplo de 3
                    # Escrever a sequência no arquivo de saída
                    fasta_arch.write(f">{contig}.{numerador}\n{nucleotide_seq}\n")
                    numerador += 1

    # Plotagem dos gráficos
    for record in SeqIO.parse(fastagc5, "fasta"):
        contig = str(record.id)
        if contig in dicioGC5:
            length = len(record.seq)
            features = []

            for interval in sorted(dicioGC5[contig]):
                inicio, final = interval.begin, interval.end - 1
                features.append(GraphicFeature(start=inicio, end=final, strand=+1, color="#ffd700"))

            record_graphic = GraphicRecord(sequence_length=length, features=features)
            ax, _ = record_graphic.plot(figure_width=5)
            plot_output_path = os.path.join(out_subfoldergc5, f'sequence_{contig}.pdf')  # Caminho de saída específico
            ax.figure.savefig(plot_output_path, bbox_inches='tight')
