import os
import glob
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.plotting import figure
from bokeh.layouts import column

# Função para encontrar arquivos FASTA com sufixos específicos
def find_orf_files(suffixes):
    orf_fasta_files = []
    for suffix in suffixes:
        files = glob.glob(f"*{suffix}")
        if files:
            orf_fasta_files.append((files[0], f"Genetic Code {suffix.split('gc')[1].split('.')[0]}"))
    return orf_fasta_files

# Função para parsear os arquivos FASTA de ORFs para diferentes códigos genéticos
def parse_orf_fastas(file_paths):
    orf_data_by_code = {}
    for file_path, code in file_paths:
        for record in SeqIO.parse(file_path, "fasta"):
            header = record.description
            parts = header.split()
            
            contig_full = parts[0]
            # Remover a parte "_ORF.26" do final do identificador da contig
            if '_ORF.' in contig_full:
                contig = contig_full[:contig_full.index('_ORF.')]
            else:
                contig = contig_full
            
            coordinates_strand = parts[1].replace('[', '').replace(']', '')
            
            if '(' in coordinates_strand:  # Verifica se '(' está presente na string
                coordinates, strand = coordinates_strand.split('(')  # Divide na '('
                coordinates = coordinates.strip()  # coordinates será "138-1149"
                strand = strand.replace(')', '').strip()  # strand será "-"
                start, end = map(int, coordinates.split('-'))
            else:
                continue  # Tratamento de erro: se a string não estiver no formato esperado

            # Extraindo informações adicionais
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
                'code': code  # Adicionar o código genético associado à ORF
            }
            
            if code not in orf_data_by_code:
                orf_data_by_code[code] = {}
            
            if contig not in orf_data_by_code[code]:
                orf_data_by_code[code][contig] = []
            
            orf_data_by_code[code][contig].append(orf)
    
    return orf_data_by_code

# Função para parsear o arquivo FASTA de nucleotídeos
def parse_nuc_fasta(file_path):
    nuc_data = {}
    for record in SeqIO.parse(file_path, "fasta"):
        nuc_data[record.id] = str(record.seq)
    return nuc_data

# Função para criar gráficos e gerar arquivo HTML para cada código genético por contig
def create_graphics(output_file, orf_data_by_code, nuc_data):
    html_content = []
    index_content = []
    
    for code_index, (code, orf_data_by_contig) in enumerate(orf_data_by_code.items(), start=1):
        index_content.append(f"<h2>{code}</h2>")
        first_contig_anchor = None
        
        for contig_index, (contig, orfs) in enumerate(orf_data_by_contig.items(), start=1):
            if contig_index == 1:
                first_contig_anchor = f"contig_{code_index}_{contig_index}"
            
            features = []
            for orf in orfs:
                strand = 1 if orf['strand'] == '+' else -1
                color = "#ffcccc" if strand == 1 else "#ccccff"
                label = f"{orf['contig_full']}: {orf['start']}-{orf['end']} ({orf['strand']}) \n Codons: start {orf['start_codon']}, stop {orf['stop_codon']}"
                feature = GraphicFeature(start=orf['start'], end=orf['end'], strand=strand, color=color, label=label)
                features.append(feature)
            
            sequence_length = len(nuc_data[contig])
            record = GraphicRecord(sequence_length=sequence_length, features=features)
            plot = record.plot_with_bokeh(figure_width=15)
            plot_html = file_html(plot, CDN, f"Contig: {contig} - Code: {code}")
            
            # Adicionar o estilo CSS à página HTML gerada
            html_content.append(f"<h3 id='contig_{code_index}_{contig_index}'>Contig: {contig} - Code: {code}</h3>")
            html_content.append(plot_html)

            # Debug: Print the HTML content length for each contig and code
            print(f"Generated HTML for contig {contig}, code {code}, length: {len(plot_html)}")
        
        index_content.append(f"<li><a href='#contig_{code_index}_1'>{code} Contigs</a></li>")
      
    with open(output_file, "w") as f:
        # Incluir estilo CSS e centralizar o conteúdo na página
        f.write(f"""<!DOCTYPE html>
                    <html>
                    <head>
                    <style>
                    body {{
                        font-family: 'Ubuntu', sans-serif;
                        display: flex;
                        justify-content: center;
                        align-items: center;
                        height: 100vh;
                    }}
                    .container {{
                        max-width: 800px;
                        padding: 20px;
                        text-align: center;
                    }}
                    .index {{
                        text-align: left;
                        margin-bottom: 20px;
                    }}
                    </style>
                    </head>
                    <body>
                    <div class="container">
                    <h1>Índice de Contigs</h1>
                    <ul class="index">{"".join(index_content)}</ul>
                    <br><br>
                    {''.join(html_content)}
                    </div>
                    </body>
                    </html>""")

# Função principal
def main(nuc_fasta_file, output_file, suffixes):
    # Encontrar arquivos ORF FASTA com os sufixos especificados
    orf_fasta_files = find_orf_files(suffixes)

    # Verificar se os arquivos fasta de ORFs existem
    for file_path, code in orf_fasta_files:
        if not os.path.isfile(file_path):
            print(f"Error: ORF FASTA file '{file_path}' for code '{code}' not found.")
            return
    if not os.path.isfile(nuc_fasta_file):
        print(f"Error: Nucleotide FASTA file '{nuc_fasta_file}' not found.")
        return

    # Parsear arquivos fasta de ORFs para diferentes códigos genéticos e contigs
    orf_data_by_code = parse_orf_fastas(orf_fasta_files)

    print("Parsing nucleotide FASTA file...")
    nuc_data = parse_nuc_fasta(nuc_fasta_file)
    print(f"Found {len(nuc_data)} nucleotide sequences.")

    print("Creating graphics and writing to HTML file...")
    create_graphics(output_file, orf_data_by_code, nuc_data)
    print(f"Plots saved in {output_file}")

# Sufixos dos arquivos FASTA de ORFs e o arquivo HTML de saída
suffixes = ['_ORFgc1.fasta', '_ORFgc5.fasta', '_ORFgc11.fasta']
nuc_fasta_file = 'teste_nonDNA.fasta'
output_file = 'orf_plots.html'

# Executa a função principal
main(nuc_fasta_file, output_file, suffixes)
