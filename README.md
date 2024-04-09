# ViewVir

ViewVir is a simple pipeline created for the visualization and analysis of viral metagenomes in RNAseq samples. The initial input required to run the pipeline is the `nonHost_sequences` folder, which should contain the assemblies of the RNAseq libraries of interest. ViewVir relies on Diamond Aligner analysis and utilizes interactive Plotly graphics for result visualization.

For more detailed information, please refer to the [Wiki page](https://github.com/gabrielvpina/ViewVir/wiki).

## Viral Data Sources - ICTV
ViewVir utilizes data from the International Committee on Taxonomy of Viruses (ICTV) to obtain the genomic composition of each virus. However, it's worth noting that this information may be outdated. The ICTV sources utilized by ViewVir are:
- Virus Metadata Feature [(VMF)](https://ictv.global/vmf)
- Master Species Lists [(MSL)](https://ictv.global/msl)


![alt text](https://github.com/gabrielvpina/my_images/blob/main/viewvir-pipe.png)


# ViewVir Output


ViewVir output files are located in the `ViewVir-results` folder. Within this folder, each submitted sample is represented by a folder named after the input sample.

There are three main outputs for each sample:
- Formatted Table (.tsv)
- Interactive Bubble Plot (HTML)
- Interactive Data Table (HTML)

![alt text](https://github.com/gabrielvpina/my_images/blob/main/viewvir-print.png)


