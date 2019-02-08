# VG

#Goal
Working on combining VCFS for drosophila chromosomes to build a VG graph.


#Directories
    1.) job_makers contains scripts to generate jobs for the hummingbird cluster
    2.) vcf_filtering contains scripts to filter vcf's made by job_makers

# Directions for running vcf_parser

python3 vcf_parser.py (see Arguments)

Arguments:

    Required
    -g (global minor allelic frequency filter)  
    Keep the row if the global minor allele frequency is above this float value

    -af (local minor allelic frequency filter)
        Keep the row if the local minor allele frequency is above this float value, where local refers 
    to the population. I.E. North America, South America, ect...

    -n (name of vcf file(s) to filter)
    Accepts a variable number of vcf files to filter

    Optional
    -an (called Chromosomes) if the percentage of chromosomes called is below this value then the line is filtered out
    -fs (Fischer strand) if the fs value is above the argument the line is filtered out
    -p if p is specified then a text file is outputted with all the lines that were not filtered, with the following fields, deliminated by tabs
        (CHROM) (POSITION) (global minor allele frequency) (Population) (# alleles for reference) (# alleles for alternative) (population minor allele frequency)
        
        Example:
        
        NT_033779.5	2000097	global_frequency = 0.05250481695568401	North America	ref = 242	 alt = 7	 min_freq = 0.028112449799196786




#Team
Alex Petrusca, Cade Mirchandani, and Marcus Strauss working under professor Corbett Detig