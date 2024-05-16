#!/bin/bash

chmod +x scripts/*.sh

scripts/./cap3.sh

scripts/./proc-cap3-out.sh

scripts/./diamondScript.sh && \

scripts/./processing-dmndTables.sh && \

python3 scripts/resultsNEW.py && \ 
