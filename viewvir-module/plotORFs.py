from Bio import SeqIO
from geneblocks import DiffBlocks, CommonBlocks

def parse_fasta(file_path):
    orf_data = []

    for record in SeqIO.parse(file_path, "fasta"):
        header = record.description
        contig = header.split()[0]
        orf_info = header.split()[1]
        coordinates, strand = orf_info.split(']')[0].replace('[', '').split('(')
        start, end = map(int, coordinates.split('-'))
        strand = strand.replace(')', '')

        orf_data.append({
            'contig': contig,
            'start': start,
            'end': end,
            'strand': strand,
            'sequence': str(record.seq)
        })

    return orf_data

# Caminho para o arquivo fasta - substitua pelo seu caminho real
fasta_file = 'res/res_ORFgc1.fasta'

# Parsear o arquivo fasta
orf_data = parse_fasta(fasta_file)

# Preparar dados para plotagem
sequence_dict = {orf['contig']: orf['sequence'] for orf in orf_data}
features = []
for orf in orf_data:
    strand = '+' if orf['strand'] == '+' else '-'
    feature = {
        'range': (orf['start'], orf['end']),
        'label': f"{orf['start']}-{orf['end']} ({strand})"
    }
    features.append(feature)

# Plotar as sequências e ORFs
plot_sequence(sequence_dict, features)

