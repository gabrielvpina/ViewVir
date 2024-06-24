from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import os
import markdown

# Função para parsear o arquivo FASTA de ORFs
def parse_orf_fasta(file_path):
    orf_data = []
    for record in SeqIO.parse(file_path, "fasta"):
        header = record.description
        contig = header.split()[0]
        orf_info = header.split()[1]

        # Extrair coordinates e strand corretamente
        coordinates_strand = orf_info.split(']')[0]  # Pega a parte até o primeiro ']'
        coordinates_strand = coordinates_strand.replace('[', '')  # Remove o '['

        if '(' in coordinates_strand:  # Verifica se '(' está presente na string
            coordinates, strand = coordinates_strand.split('(')  # Divide na '('
            coordinates = coordinates.strip()  # coordinates será "248-497"
            strand = strand.replace(')', '').strip()  # strand será "-"
            start, end = map(int, coordinates.split('-'))
        else:
            continue  # Tratamento de erro: se a string não estiver no formato esperado

        orf_data.append({
            'contig': contig,
            'start': start,
            'end': end,
            'strand': strand,
            'sequence': str(record.seq)
        })
    return orf_data

# Função para parsear o arquivo FASTA de sequências de nucleotídeos
def parse_nuc_fasta(file_path):
    nuc_data = {}
    for record in SeqIO.parse(file_path, "fasta"):
        nuc_data[record.id] = str(record.seq)
    return nuc_data

# Função para plotar as ORFs em uma sequência
def plot_orfs(sequence, orfs, output_file):
    features = []
    for orf in orfs:
        strand = 1 if orf['strand'] == '+' else -1
        color = "#ffcccc" if strand == 1 else "#ccccff"
        label = f"{orf['start']}-{orf['end']} ({orf['strand']})"
        feature = GraphicFeature(start=orf['start'], end=orf['end'], strand=strand, color=color, label=label)
        features.append(feature)

    record = GraphicRecord(sequence_length=len(sequence), features=features)
    ax, _ = record.plot(figure_width=10)
    ax.set_title(f"Contig: {orfs[0]['contig']}")
    plt.savefig(output_file)
    plt.close()

# Caminhos para os arquivos FASTA
orf_fasta_file = 'res/res_ORFgc1.fasta'
nuc_fasta_file = 'res/res_nonDNA.fasta'

# Parsear os arquivos FASTA
orf_data = parse_orf_fasta(orf_fasta_file)
nuc_data = parse_nuc_fasta(nuc_fasta_file)

# Criar diretório para armazenar os gráficos
output_dir = "orf_plots"
os.makedirs(output_dir, exist_ok=True)

# Agrupar ORFs por contig
orfs_by_contig = {}
for orf in orf_data:
    if orf['contig'] not in orfs_by_contig:
        orfs_by_contig[orf['contig']] = []
    orfs_by_contig[orf['contig']].append(orf)

# Plotar e salvar gráficos para cada contig
plot_files = []
for contig, orfs in orfs_by_contig.items():
    if contig in nuc_data:  # Certifique-se de que o contig está no arquivo de sequências nucleotídicas
        output_file = os.path.join(output_dir, f"{contig}.png")
        plot_orfs(nuc_data[contig], orfs, output_file)
        plot_files.append(output_file)

# Criar arquivo markdown com todos os gráficos
markdown_content = "# ORF Plots\n\n"
for plot_file in plot_files:
    markdown_content += f"![{os.path.basename(plot_file)}]({plot_file})\n\n"

with open("orf_plots.md", "w") as md_file:
    md_file.write(markdown_content)

# Converter o arquivo markdown para HTML
html_content = markdown.markdown(markdown_content)
with open("orf_plots.html", "w") as html_file:
    html_file.write(html_content)

print("Plots salvos e reunidos em orf_plots.md e orf_plots.html")
