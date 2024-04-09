# ViewVir

ViewVir is a pipeline created for visualization and analysis of viral metagenomes in RNAseq samples. The initial input to run the pipeline is the `nonHost_sequences` folder, which must contain the assembly of the RNAseq libraries of interest. ViewVir is based on Diamond Aligner analysis and uses interactive Plotly graphics to visualize results.

OBS - adicionar processo de importação do refseq viral e criação do banco de dados no Diamon. Add wiki


## Viral data sources - ICTV

- Virus Metadata Feature (VMR)
- Master Species Lists (MSL)

![alt text](https://github.com/gabrielvpina/my_images/blob/main/pipeline_viewvir.png)

## Processing diamond data
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

To get more information, read the [Wiki page](https://github.com/gabrielvpina/ViewVir/wiki).
