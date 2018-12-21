import os
import stat
import shutil

command = """
#!/bin/bash

java -Xmx15g -jar /hb/groups/corbettlab/DGN/gatk/gatk-package-4.0.10.1-local.jar GenotypeGVCFs \
-R /hb/groups/corbettlab/DGN/genome/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna \
-V gendb://genomic_database_{} \
-O output_{}.vcf
"""


# to use docstring type in help(create_job_command)
def create_job_command(file_counter, directory, range):
    """
    DOCTSTRING
    :param chromosome: chromosome name
    :param range: steps of 1Mb starting from 0. I.e on first iteration it will be 0 up to 1Mb
    :return: the slurm job format
    """
    file_name = directory + '/' + str(file_counter) + ".sh"
    slurm_job = command.format(range, range)
    with open(file_name, 'w') as out:
        out.write(slurm_job)
        st = os.stat(file_name)
        os.chmod(file_name, st.st_mode | stat.S_IEXEC)


for line in open('intervals.txt'):
    split = line.split()
    length = int(split[1])
    output_dir = split[:-1]
    chromosome = split[0]

    directory = chromosome[0:9] + '_' + chromosome[10]
    #if os.path.exists(directory):  #should be commented out when we put on hummmingbird or github
     #   shutil.rmtree(directory)  #should be commented out
    #os.mkdir(directory)            #should be commented out

    lower_range = 0
    file_counter = 1
    while lower_range != length:
        upper_range = min(lower_range + 1000000, length)
        range = '{}-{}'.format(lower_range + 1, upper_range)
        lower_range = upper_range
        create_job_command(file_counter, directory,range)
        #create_job_command(file_counter,chromosome directory, range)
        file_counter += 1
