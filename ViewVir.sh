#!/bin/bash

chmod +x scripts/*.sh

scripts/./diamondScript.sh && \

scripts/./processing-dmndTables.sh && \

Rscript scripts/results.R && \

rm -r diamond-processed
