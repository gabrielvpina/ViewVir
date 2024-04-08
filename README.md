![alt text](http://url/to/img.png](https://github.com/gabrielvpina/viralDiversity_analysis/blob/main/viralDiversity_analysis.png)

# Fontes dos dados - ICTV

- Virus Metadata Resource (VMR)
- Master Species Lists (MSL)
# Processando os dados do diamond

## Query basica:

`./diamond blastx -d reference -q reads.fasta -o matches.tsv`

## Query usada: 

`./diamond blastx -d viralDB.dmnd -q SRR7172360_scaffolds.fasta --outfmt '6' qseqid sseqid qlen slen pident evalue bitscore stitle full_qseq --max-target-seqs '1' --out metaviraldmd.tsv`

## Colunas separadas por tab:

`qseqid	sseqid	qlen	slen	pident	evalue	bitscore	stitle	full_qseq`	

## As colunas são adicionadas na primeira linha do arquivo de saída.

qseqid - Query Seq - id
sseqid - Subject Seq - id
qlen - Query sequence length
slen - Subject sequence length
pident - Percentage of identical matches*
evalue - Expect value
stitle - Subject Title
full_qseq - Full query sequence

## Colunas para serem inseridas no arquivo:

QuerySeq	SubjectSeq	QseqLength	SseqLength	Pident	Evalue	SubjTitle	FullQueryLength

### 1) Filtro 

`grep -i 'RdRp\|RNA-dependent\|capsid\|coat\|replicase\|glycoprotein\|replicase\|nucleoprotein\|nucleocapsid' arquivo.tsv > arquivoFiltrado.tsv`

### 2) Inserindo nova coluna

- Transformando os colchetes em tabs
- Inserção de nova coluna

QuerySeq	SubjectSeq	QseqLength	SseqLength	Pident	Evalue	bitscore	SubjTitle	Specie	FullQueryLength

### 3) Substituindo os Colchetes por espaços tabulados

`sed -i 's/\[/\t/g; s/\]/\t/g' arquivo.tsv`

### 4) Inserindo colunas na tabela

`echo "QuerySeq	SubjectSeq	QseqLength	SseqLength	Pident	Evalue	SubjTitle	Specie	FullQueryLength" | cat - metaviraldmnd.tsv > Processmetaviraldmnd.tsv && mv Processmetaviraldmnd.tsv metaviraldmnd.tsv`
