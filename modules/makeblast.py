import os
import glob
import subprocess

def blastn(outputFolder,database,cpu):

    fileFasta = os.path.join(outputFolder, "*_nonDNA.fasta")
    nonDNA_files = glob.glob(fileFasta)

    for nonDNA_file in nonDNA_files:
        infile = nonDNA_file
        sample = os.path.basename(nonDNA_file).replace("_nonDNA.fasta", "")
        outfile = os.path.join(outputFolder,sample)

        blastn_input = ["blastn","-query",infile,"-db",database,"-out",f"{outfile}_blastn.tsv","-outfmt",
        "6 qseqid qcovs pident evalue stitle","-num_threads",str(cpu),"-max_target_seqs","1"]

        subprocess.run(blastn_input, check=True)

def blastx(outputFolder,database,cpu):

    fileFasta = os.path.join(outputFolder, "*_nonDNA.fasta")
    nonDNA_files = glob.glob(fileFasta)

    for nonDNA_file in nonDNA_files:
        infile = nonDNA_file
        sample = os.path.basename(nonDNA_file).replace("_nonDNA.fasta", "")
        outfile = os.path.join(outputFolder,sample)

        blastx_input = ["blastx","-query",infile,"-db",database,"-out",f"{outfile}_blastx.tsv","-outfmt",
        "6 qseqid qcovs pident evalue stitle","-num_threads",str(cpu),"-max_target_seqs","1"]

        subprocess.run(blastx_input, check=True)

