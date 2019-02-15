import argparse
import os
import subprocess
'''
Will only work on vcfs with the same header
'''
parser = argparse.ArgumentParser(description = 'takes in a directory in your current working directory '
                                               'and merges all vcf files present')

parser.add_argument('-d', '--directory',
                    help =  'directory in current working directory that contains vcf files to be merged',
                    required = True)

args = parser.parse_args()

class Merger:
    '''
    Merger will only work on vcfs with the same header
    '''
    def __init__(self, directory):
        self.directory = directory
        self.vcfs = []
        self.cmd = 'grep -v \# {} >> temp_merged.vcf'

    def merge(self):
        os.chdir(self.directory)
        self.vcfs = [vcf for vcf in subprocess.getoutput('ls').split() if vcf.endswith('.vcf')]

        for vcf in self.vcfs:
            os.system(self.cmd.format(vcf))

        self.get_header()
        os.system(f'cat header.txt temp_merged.vcf > merged_all.vcf')
        os.system('rm header.txt temp_merged.vcf')

    def get_header(self):
        with open('header.txt', 'w') as header, open(self.vcfs[0]) as vcf:
            info_line = vcf.readline()
            header.write(info_line)

if __name__== '__main__':
    merge_vcfs = Merger(args.directory)
    merge_vcfs.merge()
