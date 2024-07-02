<img src="https://github.com/gabrielvpina/my_images/blob/main/text36800.png" width="200" height="200">

# ViewVir

ViewVir is a pipeline designed to characterize and visualize potential viral contigs of non DNA viruses in RNA-seq samples. The input for the tool is a fasta file of pre-assembled contigs. The analysis performed by ViewVir depends on pre-installed tools such as [CAP3](https://faculty.sites.iastate.edu/xqhuang/cap3-and-pcap-sequence-and-genome-assembly-programs), [ORFipy](https://github.com/urmi-21/orfipy), and [Diamond](https://github.com/bbuchfink/diamond). Additionally, it can be integrated with the local machine's [InterProScan](https://github.com/ebi-pf-team/interproscan).

For more detailed information, please refer to the [Wiki page](https://github.com/gabrielvpina/ViewVir/wiki).

## Viral Data Source
ViewVir utilizes data from the International Committee on Taxonomy of Viruses (ICTV) and NCBI to obtain the genomic composition of each virus. The sources utilized by ViewVir are:
- Virus Metadata Feature [(VMF)](https://ictv.global/vmf)
- Master Species Lists [(MSL)](https://ictv.global/msl)
- [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/)

# Usage
- `--help`: Show all options;
- `-in` or `--input`: Fasta with non-host contigs;
- `-out` or `--outdir`: Name of the output directory (creates a new if doesn't exist);
- `-vir` or `--viralDB`: RefSeq Viral Release in formatted Diamond database (.dmnd file) for blastx in diamond;
- `-scan` or `--interproscan`: Interproscan executable path (/path/to/interproscan/./interproscan.sh);
- `-cpu`: CPU usage (int);
- `-norf` or `--numORFs`: Number of ORFs selected (int);

## Example of usage

```
python ViewVir.py -in mycontigs.fasta -cpu 4 -norf 3 -vir viralDB.dmnd \
-scan /my/path/interproscan/./interproscan.sh -out test

```
# ViewVir Output
ViewVir output files are located in the `--outdir` folder.
General Outputs
- Fasta file of potential non-DNA viral sequences (nt);
- Fasta files of longest ORFs in genetic code 1,5 and 11 (AA);
- Formatted diamond Table (.tsv);
- Interproscan tables of potential non-DNA viral ORFs Conserved Domains (.tsv);
- Interactive Scatter Plot of Diamond results (HTML);
- Interactive ORF viewer (HTML) with conserved domains;

![alt text](https://github.com/gabrielvpina/my_images/blob/main/vvscreen.png)

