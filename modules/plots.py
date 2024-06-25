import os
import pandas as pd
import plotly.express as px
from plotly.offline import plot
from dna_features_viewer import GraphicFeature,GraphicRecord

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
                     hover_data=["Pident", "Evalue"])
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(yaxis={'categoryorder': 'total ascending'}, title='ViewVir interactive scatter plot')

    # Definição do nome do arquivo de saída
    pltname = os.path.join(outputFolder, "scatterPlot.html")
    plot(fig, filename=pltname, auto_open=False)





