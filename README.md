![alt text]([pipeline.png](https://github.com/gabrielvpina/my_images/blob/main/pipeline.png))

# ViewVir

ViewVir is a pipeline created for the visualization and analysis of viral metagenomes in RNAseq samples. The initial input to run the pipeline is the `unmapped_sequences` folder, which must contain the assembly of the RNAseq libraries of interest. ViewVir is based on Diamond Aligner analysis and uses interactive Plotly graphics to visualize results.

OBS - adicionar processo de importação do refseq viral e criação do banco de dados no Diamon. Add wiki


## Viral data sources - ICTV

- Virus Metadata Feature (VMR)
- Master Species Lists (MSL)

## Processing diamond data

### Basic query:

`./diamond blastx -d reference -q ler.fasta -o correspond.tsv`

### Query used:

`./diamond blastx -d viralDB.dmnd -q SRR7172360_scaffolds.fasta --outfmt '6' qseqid sseqid qlen slen pident evalue bitscore stitle full_qseq --max-target-seqs '1' --out metaviraldmd.tsv`

### Tab-separated columns:

`qseqid sseqid qlen slen pident evalue bitscore stitle full_qseq`

### Columns are added to the first line of the output file.

qseqid - Seq Query - id
sseqid - Subject Seq - id
qlen - Length of query string
slen - Length of subject string
pending - Percentage of identical matches*
evalue - expected value
stitle - Subject Title
full_qseq - Full query sequence

### Columns to be inserted into the file:

QuerySeq SubjectSeq QseqLength SseqLength Pident Evalue SubjTitle FullQueryLength

# ViewVir Output

![alt text](https://github.com/gabrielvpina/my_images/blob/main/viewvir.png)

### 1) Filter

`grep -i 'RdRp\|RNA-dependent\|capsid\|coat\|replicase\|glycoprotein\|replicase\|nucleoprotein\|nucleocapsid' file.tsv > fileFiltered.tsv`

### 2) Inserting new column

- Turning brackets into tabs
- Insertion of new column

QuerySeq SubjectSeq QseqLength SseqLength Pident Evalue bitscore SubjTitle Species FullQueryLength

### 3) Replacing square brackets with tab spaces

`sed -i 's/\[/\t/g; s/\]/\t/g' file.tsv`

### 4) Inserting columns into the table

`echo "QuerySeq SubjectSeq QseqLength SseqLength Pident Evalue SubjTitle Specie FullQueryLength" | cat - metaviraldmnd.tsv > Processmetaviraldmnd.tsv && mv Processmetaviraldmnd.tsv metaviraldmnd.tsv`
