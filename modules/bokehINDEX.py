import os
import glob
import csv
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.resources import CDN
from bokeh.embed import file_html

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

    # create ORF graphic for each genetic code
    def create_graphics(output_file, orf_data_by_code, nuc_data, cdd_data):
        output_file = os.path.join(input_dir, output_file)
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
                    
                    label = f"{orf['contig_full']}: {orf['start']}-{orf['end']} ({orf['strand']}) Length = {orf['length']} nt <br>Codons: start {orf['start_codon']}, stop {orf['stop_codon']}"

                    # Adicionar informações CDD ao label se existirem
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
                    plot_html = file_html(plot, CDN, f"Contig: {contig} - Code: {code}")

                    html_content.append(f"<h3>Contig: {contig} - Code: {code}</h3>")
                    html_content.append(plot_html)

                    print(f"Generated HTML for contig {contig}, code {code}, length: {len(plot_html)}")
        
        with open(output_file, "w") as f:
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
                        </style>
                        </head>
                        <body>
                        <div class="container">
                        {''.join(html_content)}
                        </div>
                        </body>
                        </html>""")

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

    print("Parsing nucleotide FASTA file...")
    nuc_data = parse_nuc_fasta(nuc_fasta_file)
    print(f"Found {len(nuc_data)} nucleotide sequences.")

    print("Creating graphics and writing to HTML file...")
    create_graphics(output_file, orf_data_by_code, nuc_data, cdd_data)
    print(f"Plots saved in {output_file}")
