# ViewVir

ViewVir is a pipeline designed to characterize and visualize potential viral contigs of non DNA viruses in RNA-seq samples. The input for the tool is a fasta file of pre-assembled contigs. The analysis performed by ViewVir depends on pre-installed tools such as CAP3, ORFipy, and Diamond. Additionally, it can be integrated with the local machine's InterProScan.

For more detailed information, please refer to the [Wiki page](https://github.com/gabrielvpina/ViewVir/wiki).

## Viral Data Source
ViewVir utilizes data from the International Committee on Taxonomy of Viruses (ICTV) and NCBI to obtain the genomic composition of each virus. However, this information may be outdated. The ICTV sources utilized by ViewVir are:
- Virus Metadata Feature [(VMF)](https://ictv.global/vmf)
- Master Species Lists [(MSL)](https://ictv.global/msl)
- [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/

# Usage
- `--help`: Show all options;
- `-in` or `--input`: Fasta with non-host contigs;
- `-out` or `--outdir`: Name of the output directory (creates a new if doesn't exist);
- `-vir` or `--viralDB`: Diamond database (.dmnd file) for blastx in diamond;
- `-scan` or `--interproscan`: Interproscan executable path (/path/to/interproscan/./interproscan.sh);
- `-cpu`: CPU usage <int>;
- `-norf` or `--numORFs`: <int> Number of ORFs selected;

## Example of usage

```
python ViewVir.py -in mycontigs.fasta -cpu 4 -norf 3 -vir viralDB.dmnd \
-scan /my/path/interproscan/./interproscan.sh -out teste

```
# ViewVir Output
ViewVir output files are located in the `--output` folder.
General Outputs
- Fasta of potential non-DNA viral sequences;
- Fasta files of longest ORFs in genetic code 1,5 and 11;
- Formatted diamond Table (.tsv)
- Interproscan tables of potential non-DNA viral ORFs;
- Interactive Bubble Plot (HTML)
- Interactive ORF viewer;

