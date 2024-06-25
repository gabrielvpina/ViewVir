# ViewVir

ViewVir is a simple pipeline created for the visualization and analysis of viral sequences in metatranscriptomics data. The initial input required to run the pipeline is the `non-Host contig`. ViewVir relies on Diamond Aligner and InterproScan analysis and utilizes interactive Plotly and DNA Features Viewer graphics for result visualization.

For more detailed information, please refer to the [Wiki page](https://github.com/gabrielvpina/ViewVir/wiki).

## Viral Data Source
ViewVir utilizes data from the International Committee on Taxonomy of Viruses (ICTV) and NCBI to obtain the genomic composition of each virus. However, this information may be outdated. The ICTV sources utilized by ViewVir are:
- Virus Metadata Feature [(VMF)](https://ictv.global/vmf)
- Master Species Lists [(MSL)](https://ictv.global/msl)
- [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/


# ViewVir Output


ViewVir output files are located in the `ViewVir-results` folder. Within this folder, each submitted sample is represented by a folder named based in the input sample.

General Outputs
- Fasta of potential non-DNA viral sequences;
- Fasta files of longest ORFs in genetic code 1,5 and 11;
- Formatted diamond Table (.tsv)
- Interproscan tables of potential non-DNA viral ORFs;
- Interactive Bubble Plot (HTML)
- Interactive ORF viewer;

