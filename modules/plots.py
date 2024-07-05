import os
import pandas as pd
import plotly.express as px
from plotly.offline import plot
from dna_features_viewer import GraphicFeature,GraphicRecord

def scatterPlotBLAST(outputFolder):
    # acha o arquivo processado
    for plot_file in os.listdir(outputFolder):
        if plot_file.endswith("_nonDNA.tsv"):
            inputfile_path = os.path.join(outputFolder, plot_file)
            inputfile = pd.read_csv(inputfile_path, sep='\t')
            break  # Processa apenas o primeiro arquivo encontrado

    # tabela vinda do BLASTn
    for blasttable in os.listdir(outputFolder):
        if blasttable.endswith("_blastn.tsv"):
            inputblast_path = os.path.join(outputFolder,blasttable)
            inputblast = pd.read_csv(inputblast_path, sep='\t')
            break

    # tabela vinda do BLASTx
    for blastxtable in os.listdir(outputFolder):
        if blastxtable.endswith("_blastx.tsv"):
            inputblastx_path = os.path.join(outputFolder,blastxtable)
            inputblastx = pd.read_csv(inputblastx_path, sep='\t')
            break

    # texto principal
    inputfile["MatchSequence"] = inputfile["Species"] + " --> " + inputfile["SubjTitle"]

    inputfile["Family"] = inputfile["Family"].fillna("unknownFamily")
    ################### BLASTn
    inputblast.columns = ['QuerySeq','BLASTn_Cover','BLASTn_Ident','BLASTn_evalue','BLASTn-stitle']
    inputfile = inputfile.merge(inputblast, on='QuerySeq', how='left')
    inputfile["BLASTn"] = (
        inputfile['QuerySeq'] + "| Cover:" + inputfile['BLASTn_Cover'].astype(str) +
        "| Ident:" + inputfile['BLASTn_Ident'].astype(str) + 
        "| E-value:" + inputfile['BLASTn_evalue'].astype(str) + 
        "| Title:" + inputfile['BLASTn-stitle']
    )
    inputfile["BLASTn"] = inputfile["BLASTn"].fillna("No matches on BLASTn")
    ###################

    ################### BLASTx
    inputblastx.columns = ['QuerySeq','BLASTx_Cover','BLASTx_Ident','BLASTx_evalue','BLASTx-stitle']
    inputfile = inputfile.merge(inputblastx, on='QuerySeq', how='left')
    inputfile["BLASTx"] = (
        inputfile['QuerySeq'] + "| Cover:" + inputfile['BLASTx_Cover'].astype(str) +
        "| Ident:" + inputfile['BLASTx_Ident'].astype(str) + 
        "| E-value:" + inputfile['BLASTx_evalue'].astype(str) + 
        "| Title:" + inputfile['BLASTx-stitle']
    )
    inputfile["BLASTx"] = inputfile["BLASTx"].fillna("No matches on BLASTn")
    ###################


    inputfile["Genome.composition"] = inputfile["Genome.composition"].fillna("unknown")
    inputfile["Family_Info"] = inputfile["Family"] + " - " + inputfile["Genome.composition"]

    fig = px.scatter(inputfile, x="QseqLength", y="MatchSequence", size="QCover", color="Family_Info",
                     hover_data=["Pident", "Evalue","QuerySeq","BLASTn","BLASTx"])
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(yaxis={'categoryorder': 'total ascending'}, title='ViewVir interactive scatter plot')

    # Definição do nome do arquivo de saída
    pltname = os.path.join(outputFolder, "scatterPlot.html")
    plot(fig, filename=pltname, auto_open=False)

    # Tabela de curadoria manual
    csv_output_path = os.path.join(outputFolder, "ViewVir-blasts_table.csv")
    inputfile.to_csv(csv_output_path, index=False)

def scatterPlot(outputFolder):
    # Localiza o arquivo correto
    for plot_file in os.listdir(outputFolder):
        if plot_file.endswith("_nonDNA.tsv"):
            inputfile_path = os.path.join(outputFolder, plot_file)
            inputfile = pd.read_csv(inputfile_path, sep='\t')
            break  # Processa apenas o primeiro arquivo encontrado

    # Gráfico Interativo
    inputfile["MatchSequence"] = inputfile["Species"] + " --> " + inputfile["SubjTitle"]
    inputfile["Genome.composition"] = inputfile["Genome.composition"].fillna("NA")

    fig = px.scatter(inputfile, x="QseqLength", y="MatchSequence", size="QCover", color="Genome.composition",
                     hover_data=["Pident", "Evalue","QuerySeq"])
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(yaxis={'categoryorder': 'total ascending'}, title='ViewVir interactive scatter plot')

    # Definição do nome do arquivo de saída
    pltname = os.path.join(outputFolder, "scatterPlot.html")
    plot(fig, filename=pltname, auto_open=False)





