import os
import glob
import csv
import pandas as pd
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord

from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import Label, Range1d
from bokeh.plotting import figure, show
from bokeh.io import output_notebook
from bokeh.models import HoverTool, ColumnDataSource

import plotly.express as px
import plotly.io as pio
from plotly.offline import plot


suffixes = ['_ORFgc1.fasta', '_ORFgc5.fasta', '_ORFgc11.fasta']
out_html = "orf_plots.html" 


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
            try:
                inputblast = pd.read_csv(inputblast_path, sep='\t')
            except pd.errors.EmptyDataError:
                inputblast = pd.DataFrame(columns=['QuerySeq', 'BLASTn_Cover', 'BLASTn_Ident', 'BLASTn_evalue', 'BLASTn_stitle'])
            break

    # tabela vinda do BLASTx
    for blastxtable in os.listdir(outputFolder):
        if blastxtable.endswith("_blastx.tsv"):
            inputblastx_path = os.path.join(outputFolder,blastxtable)
            try:
                inputblastx = pd.read_csv(inputblastx_path, sep='\t')
            except pd.errors.EmptyDataError:
                inputblastx = pd.DataFrame(columns=['QuerySeq', 'BLASTx_Cover', 'BLASTx_Ident', 'BLASTx_evalue', 'BLASTx_stitle'])

            break

    # texto principal
    inputfile["MatchSequence"] = inputfile["Species"] + " --> " + inputfile["SubjTitle"]
    inputfile["Family"] = inputfile["Family"].fillna("unknownFamily")
    ################### BLASTn
    inputblast.columns = ['QuerySeq','BLASTn_Cover','BLASTn_Ident','BLASTn_evalue','BLASTn-stitle']
    inputfile = inputfile.merge(inputblast, on='QuerySeq', how='left')
    inputfile["BLASTn"] = (
        "Cover: " + inputfile['BLASTn_Cover'].astype(str) +
        " | Ident: " + inputfile['BLASTn_Ident'].astype(str) + 
        " | E-value: " + inputfile['BLASTn_evalue'].astype(str) + 
        " | Title: " + inputfile['BLASTn-stitle']
    )
    inputfile["BLASTn"] = inputfile["BLASTn"].fillna("No matches on BLASTn")
    ###################

    ################### BLASTx
    inputblastx.columns = ['QuerySeq','BLASTx_Cover','BLASTx_Ident','BLASTx_evalue','BLASTx-stitle']
    inputfile = inputfile.merge(inputblastx, on='QuerySeq', how='left')
    inputfile["BLASTx"] = (
        "Cover: " + inputfile['BLASTx_Cover'].astype(str) +
        " | Ident: " + inputfile['BLASTx_Ident'].astype(str) + 
        " | E-value: " + inputfile['BLASTx_evalue'].astype(str) + 
        " | Title: " + inputfile['BLASTx-stitle']
    )
    inputfile["BLASTx"] = inputfile["BLASTx"].fillna("No matches on BLASTn")
    ###################


    inputfile["Genome.composition"] = inputfile["Genome.composition"].fillna("unknown")
    inputfile["Family_Info"] = inputfile["Family"] + " - " + inputfile["Genome.composition"]
    inputfile["ID"] = inputfile["QuerySeq"].apply(lambda x: f"contig_{x}")

    fig = px.scatter(inputfile, x="QseqLength", y="MatchSequence", size="QCover", color="Family_Info",
                     hover_data=["Pident", "Evalue", "QuerySeq", "BLASTn", "BLASTx"], custom_data=["ID"])

    fig.update_yaxes(showticklabels=False)
    #fig.update_layout(yaxis={'categoryorder': 'total ascending'}, title='ViewVir interactive scatter plot', template="simple_white")
    fig.update_layout(template="none", xaxis_title="Sequence Lenght (nt)", yaxis_title=None)
    fig.update_traces(textposition='top center')

    # Salva o scatterplot em um arquivo HTML
    pltname = os.path.join(outputFolder, "scatterPlot.html")
    plot(fig, filename=pltname, auto_open=False)

    # Salva a tabela para curadoria manual
    csv_output_path = os.path.join(outputFolder, "ViewVir-blasts_table.csv")
    inputfile.to_csv(csv_output_path, index=False)

    return pltname

def scatterPlot(outputFolder):
    # Localiza o arquivo correto
    for plot_file in os.listdir(outputFolder):
        if plot_file.endswith("_nonDNA.tsv"):
            inputfile_path = os.path.join(outputFolder, plot_file)
            inputfile = pd.read_csv(inputfile_path, sep='\t')
            break  # Processa apenas o primeiro arquivo encontrado

     # texto principal
    inputfile["MatchSequence"] = inputfile["Species"] + " --> " + inputfile["SubjTitle"]
    inputfile["Family"] = inputfile["Family"].fillna("unknownFamily")
    inputfile["Genome.composition"] = inputfile["Genome.composition"].fillna("unknown")
    inputfile["Family_Info"] = inputfile["Family"] + " - " + inputfile["Genome.composition"]
    inputfile["ID"] = inputfile["QuerySeq"].apply(lambda x: f"contig_{x}")

    fig = px.scatter(inputfile, x="QseqLength", y="MatchSequence", size="QCover", color="Family_Info",
                     hover_data=["Pident", "Evalue", "QuerySeq"], custom_data=["ID"])

    fig.update_yaxes(showticklabels=False)
    #fig.update_layout(yaxis={'categoryorder': 'total ascending'}, title='ViewVir interactive scatter plot', template="simple_white")
    fig.update_layout(template="none", xaxis_title="Sequence Lenght (nt)", yaxis_title=None)
    fig.update_traces(textposition='top center')

    # Salva o scatterplot em um arquivo HTML
    pltname = os.path.join(outputFolder, "scatterPlot.html")
    plot(fig, filename=pltname, auto_open=False)


def generate_orf_plots(input_dir, output_file, suffixes):
    # find fasta files with suffix
    def find_orf_files(suffixes):
        orf_fasta_files = []
        for suffix in suffixes:
            files = glob.glob(os.path.join(input_dir, f"*{suffix}"))
            if files:
                orf_fasta_files.append((files[0], f"Genetic Code {suffix.split('gc')[1].split('.')[0]}"))
        return orf_fasta_files

    # parse different genetic codes
    def parse_orf_fastas(file_paths):
        orf_data_by_code = {}
        for file_path, code in file_paths:
            for record in SeqIO.parse(file_path, "fasta"):
                header = record.description
                parts = header.split()
                
                contig_full = parts[0]
                if '_ORF.' in contig_full:
                    contig = contig_full[:contig_full.index('_ORF.')]
                else:
                    contig = contig_full
                
                coordinates_strand = parts[1].replace('[', '').replace(']', '')
                if '(' in coordinates_strand:
                    coordinates, strand = coordinates_strand.split('(')
                    coordinates = coordinates.strip()
                    strand = strand.replace(')', '').strip()
                    start, end = map(int, coordinates.split('-'))
                else:
                    continue

                additional_info = parts[2:]
                additional_info_dict = {info.split(':')[0]: info.split(':')[1] for info in additional_info}

                orf = {
                    'contig_full': contig_full,
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'sequence': str(record.seq),
                    'type': additional_info_dict.get('type'),
                    'length': int(additional_info_dict.get('length', 0)),
                    'frame': additional_info_dict.get('frame'),
                    'start_codon': additional_info_dict.get('start'),
                    'stop_codon': additional_info_dict.get('stop'),
                    'code': code
                }
                
                if contig not in orf_data_by_code:
                    orf_data_by_code[contig] = []
                
                orf_data_by_code[contig].append(orf)
        
        return orf_data_by_code

    
    def parse_nuc_fasta(file_path):
        nuc_data = {}
        for record in SeqIO.parse(file_path, "fasta"):
            nuc_data[record.id] = str(record.seq)
        return nuc_data

    # parse interproscan files
    def parse_cdd_files(input_dir):
        cdd_data = {}
        cdd_suffixes = ['_ORFgc1_CDD.tsv', '_ORFgc5_CDD.tsv', '_ORFgc11_CDD.tsv']
        for suffix in cdd_suffixes:
            files = glob.glob(os.path.join(input_dir, f"*{suffix}"))
            for file in files:
                with open(file, 'r') as f:
                    reader = csv.reader(f, delimiter='\t')
                    for row in reader:
                        orf_name = row[0]
                        database = row[3]
                        domain_code = row[4]
                        domain_name = row[5]
                        e_value = row[8]

                        # remove blank cells
                        if domain_name != '-' and e_value != '-' and database in ['CDD', 'Pfam', 'SUPERFAMILY']:
                            if orf_name not in cdd_data:
                                cdd_data[orf_name] = []
                            cdd_data[orf_name].append({
                                'database': database,
                                'domain_code': domain_code,
                                'domain_name': domain_name,
                                'e_value': e_value
                            })
        return cdd_data

# ORF plot

    def create_graphics(output_file, orf_data_by_code, nuc_data, cdd_data):
        output_file = os.path.join(input_dir, out_html)
        html_content = []

        for contig, orfs in orf_data_by_code.items():
            for suffix in suffixes:
                features = []
                code = f"Genetic Code {suffix.split('gc')[1].split('.')[0]}"
                for orf in orfs:
                    if orf['code'] != code:
                        continue
                    strand = 1 if orf['strand'] == '+' else -1
                    color = "#ffcccc" if strand == 1 else "#ccccff"
                    
                    label = f"{orf['contig_full']}: {orf['start']}-{orf['end']} ({orf['strand']}), Frame ({orf['frame']}), Length = {orf['length']} nt <br>Codons: start {orf['start_codon']}, stop {orf['stop_codon']}"

                    # Adiciona informações CDD ao label se existirem
                    if orf['contig_full'] in cdd_data:
                        cdd_info = cdd_data[orf['contig_full']]
                        for cdd in cdd_info:
                            label += f"<br>{cdd['database']}: {cdd['domain_code']} {cdd['domain_name']} (e-value: {cdd['e_value']})"
                    

                    feature = GraphicFeature(start=orf['start'], end=orf['end'], strand=strand, color=color, label=label)
                    features.append(feature)
            
                if features:
                    sequence_length = len(nuc_data[contig])
                    record = GraphicRecord(sequence_length=sequence_length, features=features)
                    plot = record.plot_with_bokeh(figure_width=15)
                    plot_html = file_html(plot, CDN, f"{contig} - {code}")
                    html_content.append(f"<h3 id='{contig}'>{contig} - {code}</h3>")
                    html_content.append(plot_html)


##################################################################################################
            html_orf = f"""<!DOCTYPE html>
                <html>
                <head>
                    <style>
                        body {{
                            font-family: 'Ubuntu', sans-serif;
                            display: flex;
                            justify-content: center;
                            flex-direction: column;
                            margin: 0;
                            padding: 0;
                            height: auto;
                            overflow-x: hidden;
                        }}
                        .search-bar {{
                            width: 50%;
                            margin: 20px auto;
                            display: flex;
                        }}
                        .search-bar input[type="text"] {{
                            width: 100%;
                            padding: 10px;
                            font-size: 16px;
                            border: 1px solid #ccc;
                            border-radius: 4px;
                            box-sizing: border-box;
                        }}
                        .search-bar input[type="submit"] {{
                            padding: 10px;
                            font-size: 16px;
                            margin-left: 10px;
                            background-color: #031d54;
                            color: white;
                            border: none;
                            border-radius: 4px;
                            cursor: pointer;
                        }}
                        .container {{
                            display: flex;
                            width: 100%;
                            height: auto;
                            justify-content: center;
                            padding: 30px;
                            box-sizing: border-box;
                        }}
                        h3 {{
                            text-align: center; /* Centraliza todos os títulos de contigs */
                        }}
                    </style>
                </head>
                <body>
                    <!-- Barra de Pesquisa -->
                    <div class="search-bar">
                        <input type="text" id="searchInput" placeholder="Search for a contig..." onkeypress="handleEnter(event)">
                        <input type="submit" value="Search" onclick="searchFunction()">
                    </div>

                    <!-- Título menor centralizado acima da primeira contig -->
                    <h3 id="main-title">ORFs Info Visualization</h3>

                    <div class="container">
                        <div class="center">
                            {''.join(html_content[0:])}
                        </div>
                    </div>

                    <!-- Função JavaScript para a pesquisa -->
                    <script>
                        // Função para buscar contigs
                        function searchFunction() {{
                            var input, filter, contigs, i, txtValue;
                            input = document.getElementById('searchInput');
                            filter = input.value.toUpperCase();
                            contigs = document.querySelectorAll('h3'); // Assume que contigs estão em elementos <h3>

                            let found = false;
                            // Itera sobre os títulos de contigs e ORFs
                            for (i = 0; i < contigs.length; i++) {{
                                txtValue = contigs[i].textContent || contigs[i].innerText;
                                // Checa se o texto contém o termo pesquisado
                                if (txtValue.toUpperCase().indexOf(filter) > -1) {{
                                    if (!found) {{
                                        // Rola até o primeiro contig correspondente
                                        contigs[i].scrollIntoView({{ behavior: 'smooth', block: 'start' }});
                                        found = true;
                                    }}
                                }}
                            }}

                            if (!found) {{
                                alert("No contig found.");
                            }}
                        }}

                        // Função para detectar a tecla Enter na barra de pesquisa
                        function handleEnter(event) {{
                            if (event.keyCode === 13) {{ // Código da tecla Enter
                                searchFunction(); // Chama a função de busca
                            }}
                        }}
                    </script>
                </body>
                </html>"""



        
        with open(output_file, "w") as f:
            f.write(html_orf)
##################################################################################################
    # find orfs fasta
    orf_fasta_files = find_orf_files(suffixes)

    # debugs
    for file_path, code in orf_fasta_files:
        if not os.path.isfile(file_path):
            print(f"Error: ORF FASTA file '{file_path}' for code '{code}' not found.")
            return

    
    nuc_fasta_file = None
    for nuc in os.listdir(input_dir):
        if nuc.endswith("_nonDNA.fasta"):
            nuc_fasta_file = os.path.join(input_dir, nuc)
            break

    if not nuc_fasta_file or not os.path.isfile(nuc_fasta_file):
        print(f"Error: Nucleotide FASTA file not found in the specified directory.")
        return

    # parse fasta files of all codes
    orf_data_by_code = parse_orf_fastas(orf_fasta_files)

    # parse CDD files
    cdd_data = parse_cdd_files(input_dir)

    # create out file
    output_file = os.path.join(input_dir,out_html)

    print("Parsing nucleotide FASTA file...")
    nuc_data = parse_nuc_fasta(nuc_fasta_file)
    print(f"Found {len(nuc_data)} nucleotide sequences.")

    print("Creating graphics and writing to HTML file...")
    create_graphics(output_file, orf_data_by_code, nuc_data, cdd_data)
    print(f"Plots saved in {output_file}")


def combine_html(scatterplot_html, orf_plots_html, output_file, image_logo):
    with open(scatterplot_html, 'r') as f:
        scatterplot_content = f.read()
        
    with open(orf_plots_html, 'r') as f:
        orf_plots_content = f.read()

    combined_html = f"""<!DOCTYPE html>
        <html>
        <head>
            <!-- Título da aba do navegador -->
            <title>ViewVir - Viral Diversity Analysis</title>

            <!-- Favicon (logo na aba do navegador) -->
            <link rel="icon" href="{image_logo}" type="image/png">

            <!-- Importa a fonte do Plotly -->
            <link href="https://fonts.googleapis.com/css2?family=Open+Sans:wght@400;600&display=swap" rel="stylesheet">
            <style>
                body {{
                    font-family: 'Open Sans', sans-serif;
                    justify-content: center;
                    display: flex;
                    flex-direction: column;
                    margin: 0;
                    padding: 0;
                }}
                .header {{
                    background-color: #031d54; /* Fundo azul */
                    width: 100%;
                    padding: 20px 40px;
                    box-sizing: border-box;
                    display: flex;
                    align-items: center;
                }}
                .logo {{
                    width: 80px; /* Define o tamanho da logo */
                    height: 80px;
                    margin-right: 20px; /* Espaço entre a logo e o título */
                }}
                h1 {{
                    font-weight: 600; /* Peso da fonte do Plotly */
                    color: white; /* Cor da fonte branca para contraste */
                    margin: 0;
                    text-align: left;
                }}
                h2 {{
                    font-weight: 400; /* Subtítulo com peso menor */
                    color: white;
                    margin: 5px 0 0 0;
                    font-size: 16px;
                    text-align: left;
                }}
                a {{
                    color: white; /* Cor do link no subtítulo */
                    text-decoration: underline;
                }}
                .scatterplot {{
                    width: 100%;
                    height: 50vh;
                    overflow: hidden;
                }}
                .orf-plots {{
                    width: 100%;
                    height: 50vh;
                    justify-content: center;
                    overflow: auto;
                }}
            </style>
        </head>
        <body>
            <!-- Cabeçalho com logo, título e subtítulo -->
            <div class="header">
                <img src="{image_logo}" alt="ViewVir Logo" class="logo">
                <div>
                    <h1>ViewVir - A pipeline for viral diversity analysis</h1>
                    <h2>For further information, please refer to the documentation available on GitHub: <a href="https://github.com/gabrielvpina/ViewVir" target="_blank">https://github.com/gabrielvpina/ViewVir</a></h2>
                </div>
            </div>

            <div class="scatterplot">
                {scatterplot_content}
            </div>
            <div class="orf-plots">
                {orf_plots_content}
            </div>
            <script>
                document.querySelectorAll('.scatterplot .scatter-mark').forEach(function(elem) {{
                    elem.addEventListener('click', function() {{
                        var contigId = elem.getAttribute('text');
                        document.getElementById(contigId).scrollIntoView();
                    }});
                }});
            </script>
        </body>
        </html>
        """



    with open(output_file, 'w') as f:
        f.write(combined_html)







