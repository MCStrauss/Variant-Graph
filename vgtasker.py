command = \
"""#!/bin/bash

java -Xmx15g -jar /hb/groups/corbettlab/DGN/gatk/gatk-package-4.0.10.1-local.jar GenomicsDBImport \\
--sample-name-map /hb/groups/corbettlab/DGN/sample_map.txt \\
--reader-threads 20 \\
--genomicsdb-workspace-path /hb/groups/corbettlab/DGN/{}/genomic_database_{} \\
--intervals {}:{}
"""

# filename = NT_033779_5
# range = 1-1000000
def create_job_command(filename, range):
    filename2 = filename[0:9] + '_' + filename[-1]
    slurm_job = command.format(filename2, range, filename, range)
    return slurm_job

for line in open('intervals.txt'):
    split = line.split()
    filename = split[0]
    length = int(split[1])
    with open('out.txt', 'w') as out:
        val = 0
        while val != length:
            newVal = min(val + 1000000, length)
            range = '{}-{}'.format(val + 1, newVal)
            val = newVal
            out.write(create_job_command(filename, range))
