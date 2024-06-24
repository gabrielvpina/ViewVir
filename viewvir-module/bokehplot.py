import os
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.plotting import figure, show, output_file
from bokeh.layouts import column

# Função para parsear o arquivo FASTA de ORFs
def parse_orf_fasta(file_path):
    orf_data = []
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

        orf_data.append({
            'contig': contig,
            'start': start,
            'end': end,
            'strand': strand,
            'sequence': str(record.seq),
            'type': additional_info_dict.get('type'),
            'length': int(additional_info_dict.get('length', 0)),
            'frame': additional_info_dict.get('frame'),
            'start_codon': additional_info_dict.get('start'),
            'stop_codon': additional_info_dict.get('stop')
        })
    return orf_data

def parse_nuc_fasta(file_path):
    nuc_data = {}
    for record in SeqIO.parse(file_path, "fasta"):
        nuc_data[record.id] = str(record.seq)
    return nuc_data

def create_graphics(output_file, orfs_by_contig, nuc_data):
    html_content = []
    for contig, orfs in orfs_by_contig.items():
        if contig not in nuc_data:
            print(f"Contig {contig} not found in nucleotide data, skipping...")
            continue
        
        features = []
        for orf in orfs:
            strand = 1 if orf['strand'] == '+' else -1
            color = "#ffcccc" if strand == 1 else "#ccccff"
            label = f"{orf['start']}-{orf['end']} ({orf['strand']})"
            feature = GraphicFeature(start=orf['start'], end=orf['end'], strand=strand, color=color, label=label)
            features.append(feature)
        
        sequence_length = len(nuc_data[contig])
        record = GraphicRecord(sequence_length=sequence_length, features=features)
        plot = record.plot_with_bokeh(figure_width=15)
        plot_html = file_html(plot, CDN, f"Contig: {contig}")
        html_content.append(plot_html)

        # Debug: Print the HTML content length for each contig
        print(f"Generated HTML for contig {contig}, length: {len(plot_html)}")
      
    with open(output_file, "w") as f:
        f.write("<br><br>".join(html_content))
        
def main(orf_fasta_file, nuc_fasta_file, output_file):
    if not os.path.isfile(orf_fasta_file):
        print(f"Error: ORF FASTA file '{orf_fasta_file}' not found.")
        return
    if not os.path.isfile(nuc_fasta_file):
        print(f"Error: Nucleotide FASTA file '{nuc_fasta_file}' not found.")
        return

    print("Parsing ORF FASTA file...")
    orf_data = parse_orf_fasta(orf_fasta_file)
    print(f"Found {len(orf_data)} ORFs.")

    print("Parsing nucleotide FASTA file...")
    nuc_data = parse_nuc_fasta(nuc_fasta_file)
    print(f"Found {len(nuc_data)} nucleotide sequences.")

    print("Grouping ORFs by contig...")
    orfs_by_contig = {}
    for orf in orf_data:
        if orf['contig'] not in orfs_by_contig:
            orfs_by_contig[orf['contig']] = []
        orfs_by_contig[orf['contig']].append(orf)

    print("Creating graphics and writing to HTML file...")
    create_graphics(output_file, orfs_by_contig, nuc_data)
    print(f"Plots saved in {output_file}")

# Caminhos para os arquivos FASTA e o arquivo HTML de saída
orf_fasta_file = 'teste_ORFgc1.fasta'
nuc_fasta_file = 'teste_nonDNA.fasta'
output_file = 'orf_plots.html'

# Executa a função principal
main(orf_fasta_file, nuc_fasta_file, output_file)
