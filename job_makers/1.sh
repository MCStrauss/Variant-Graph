#!/bin/bash

java -Xmx15g -jar /hb/groups/corbettlab/DGN/gatk/gatk-package-4.0.10.1-local.jar GenomicsDBImport \
--sample-name-map /hb/groups/corbettlab/DGN/sample_map.txt \
--reader-threads 20 \
--genomicsdb-workspace-path /hb/groups/corbettlab/DGN/NT_033779_5/genomic_database_1-1000000 \
--intervals NT_033779.5:1-1000000 \
