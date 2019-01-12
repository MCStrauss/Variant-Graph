#!/bin/bash

### REQUIRED
#SBATCH -p 128x24                  ### everything we do can be on partition instruction
#SBATCH -N 1                            ### number of nodes ot use
#SBATCH -n 1                            ### number of CPUs to use
#SBATCH --mem=16500                      ### amount of memory to use
#SBATCH -t 0                           ### run time in minutes
#SBATCH -J genotype_test
#SBATCH -o job%.j.out

java -Xmx15g -jar /hb/groups/corbettlab/DGN/gatk/gatk-package-4.0.10.1-local.jar GenotypeGVCFs \
-R /hb/groups/corbettlab/DGN/genome/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna \
-V gendb://genomic_database_2000001-3000000 \
-O test_output_2000001-3000000.vcf