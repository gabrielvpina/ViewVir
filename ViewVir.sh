#!/bin/bash

chmod +x scripts/*.sh


scripts/./cap3.sh

scripts/./proc-cap3-out.sh

scripts/./diamondScript.sh && \

scripts/./processing-dmndTables.sh && \

Rscript scripts/results.R && \

scripts/./criando_fasta.sh && \

rm -r diamond-processed
