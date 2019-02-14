import vcf
import argparse
from count import Frequency
from populations import populations
import os
import csv
from itertools import dropwhile

parser = argparse.ArgumentParser(description  = 'Arguments for VCF filter script')

parser.add_argument('-an', '--an_cutoff',
                    help = 'AN count must represent certain %% of total',
                    type = float) #how good the site is overall

parser.add_argument('-af', '--allelic_frequency',
                    help = '%% to filter for local allelic frequency. '
                    'If a read has no sample with greater than x%% then '
                    'it is removed \n from the vcf. ',
                    type = float, required = True)

parser.add_argument('-p', '--population',
                    help = 'If -p is specified an output file with population and '
                            'their minor allelic frequencies is outputed',
                    action = 'store_true')

parser.add_argument('-g', '--global_frequency',
                    help = 'gives the cutoff for the global frequency of a position on the chromosome,'
                            'if it is less than the cutoff it will not be included in the filtered vcf',
                    type = float,
                    required = True)

parser.add_argument('-fs', '--fischer_strand',
                    help = 'value to filter fischer strand on, if the fischer strand on an '
                            'chromosome is larger than the fs, it will be filtered out',
                    type  = float)

parser.add_argument('-n', '--name',
                    help = 'Name of VCF files to be filtered, can input multiple files',
                    required = True,
                    nargs = '+')

args = parser.parse_args()

fieldnames = ['CHROM', 'POSITION', 'Population', 'Global Minor Allelic Frequency', 'Reference Count',
                          'Alternative Count', 'Local Minor Allelic Frequency']

header = ''

def check():
    if args.an_cutoff: assert 0 <= args.an_cutoff < 1, 'an_cutoff must be atleast 0 and greater than 1'
    if args.allelic_frequency: assert 0 <= args.allelic_frequency < 1, 'allelic_frequency must be atleast 0 and greater than 1'
    if args.fischer_strand: assert args.fischer_strand > 0 , 'fischer strand filter must be greater than 0'
    if args.global_frequency: assert 0 < args.global_frequency <1, 'global frequency must be between 0 and 1'

class Parser:

    def __init__(self, vcf_reader, name, direct):
        '''
        :param vcf_reader: the vcf file we will be iterating over
        :param f_aan: cutoff for mininum number of alleles
        :param f_maf: cutoff for minor allelic frequency
        '''
        self.vcf_reader = vcf_reader
        self.dB = {} #key = populations val = Frequency object
        self.an_cutoff = .5
        self.f_maf = args.allelic_frequency# threshold to filter for minor allele
        self.fs_filter = 10  # default filter value for fischer strand
        self.glob_frequency_filter = args.global_frequency # threshold the filter a line in referenece to the minor allelic frequency fo the entire line
        self.directory = direct

        self.read_from = name #read from this file
        self.output = f'filtered_{name}'  #name of generated VCF
        self.path  = direct + '/' + self.output





    def get_arguments(self):

        if args.fischer_strand:
            self.fs_filter = args.fischer_strand

        if args.fischer_strand:
            self.fs_filter = args.fischer_strand

        if args.an_cutoff:
            self.an_cutoff = args.an_cutoff

    @staticmethod
    def parse(data):
        '''
        :param data: data will be in the form '#/#'
        :return: if there is atleast one '.' the read is bad and do not increment either the ref or alt.
        if there is atleast one 1, increment the alt by 1. if 0/0 increment the ref by 1. Similarly a read of 1/1 will
        only increment the alt by 1.
        '''
        left, right = data.split('/')
        if left == '.' or right == '.':
            return (0, 0)

        left, right = eval(left), eval(right)

        if left or right:
            return (0, 1)
        else:
            return (1, 0)

    def filter_AN_FS(self, line):
        '''
        :param line: a row in a vcf file
        :return:  if the an doesn't have the desired number of samples we want remove the row from our new vcf file.
        '''

        fs = eval(''.join([x for x in line.split(';')[5:7] if x.startswith('FS')]).split('=')[1]) #grabs fischer strand
        an = eval(line.split(';')[2].split('=')[1])  # grabs the an from the line

        if (an / 2) / 1070 < self.an_cutoff or fs > self.fs_filter: #1070 is the total number of samples, an represents the number of chromosomes in the record
            return False
        return True

    def filter_line(self, record):
        '''
        :param chrom: given the (chrom, position) check everywhere it occured and its minor allelic frequency.
        :return: If it has a minor allelic frequency of atleast our cutoff keep it in, else remove the line from our filtered VCF.
        '''
        if max(record.aaf) > .5:  # global minor allele is uniform across populations
            glob_freq = 1 - max(record.aaf)
        else:
            glob_freq = max(record.aaf)

        if glob_freq >= self.glob_frequency_filter:
            return True

        for population in self.dB:
            loc_freq = self.dB[population].minor_freq # minor allelic frequency relative to population on chrom, position
            if loc_freq >= self.f_maf:
                return True
        return False

    def gen_record_data(self, record):
        '''
        :return: returns nothing, this generates our dB,
        '''
        for sample in record.samples:

            x, y = self.parse(sample['GT'])

            for func, pop in populations.items(): #populations is a dictionary of functions and values where the values are our samples and the functions determine which sample corresponds to the sample header from the vcf.
                if func(sample.sample):
                    if pop in self.dB:
                        self.dB[pop].update(ref = x, alt = y) #already in the dB, update it
                        break
                    else:
                        self.dB[pop] = Frequency(ref = x, alt = y) #needs to be initialized in the dB
                        break

        return

    def check_right_maf(self):
        '''
        checks that each minor allelic frequency is less than 50%, if not reverses it by subtracting it from 1
        :return: None
        '''

        for pop in self.dB:
            if self.dB[pop].minor_freq > .5:
                major_allelic = self.dB[pop].minor_freq
                self.dB[pop].minor_freq = 1 - major_allelic

        return

    def write_data(self, record):
        if max(record.aaf) > .5:  # global minor allele is uniform across populations
            glob_freq = 1 - max(record.aaf)
        else:
            glob_freq = max(record.aaf)
        out = self.directory + '/' + 'populations_' + self.read_from.replace('.vcf', '.csv')

        with open(out, mode='a') as csv_file:
            for population in self.dB:
                population_writer = csv.DictWriter(csv_file, fieldnames = fieldnames, delimiter='\t')
                population_writer.writerow({'CHROM': record.CHROM, 'POSITION': record.POS, 'Population': population,
                                            'Global Minor Allelic Frequency': str(glob_freq), 'Reference Count': str(self.dB[population].ref),
                                            'Alternative Count': str(self.dB[population].alt),
                                            'Local Minor Allelic Frequency': str(self.dB[population].minor_freq),})


    def generate_filtered_vcf(self):
        '''
        :param name: name of the new vcf.txt
        :return: will generate a filtered vcf if it fits aan parameters and atleast 1 chromosome has an minor alleliccd
        frequency above .05% or the argument the user specified, relies on two helpers function filter_local_frequency and filter_AAN.
        '''
        with open(self.read_from, 'r') as input, open(self.path, 'w') as out:

            for line in dropwhile(lambda x: x.startswith('##'), input):

                if line.startswith('#'):
                    out.write(f'{line}')
                    global header
                    header = line
                else:
                    split = line.split('\t') #fischer strand

                    record = next(self.vcf_reader) #todo stop iteration

                    self.gen_record_data(record) #generates the data
                    self.check_right_maf()  # makes sure the minor allelic freq is actually the minor allelic frequency

                    if self.filter_AN_FS(split[7]) and self.filter_line(record): #gets chrom name as string and position as int
                        out.write(f'{line}')

                        if args.population:
                            self.write_data(record)

                    self.dB.clear()



if __name__ == '__main__':
    check() #checks to make sure args are valid if there are args
    direct = f'filtered_maf-{args.allelic_frequency}_global-{args.global_frequency}'

    if not os.path.exists(direct):
        os.mkdir(direct)

    for name in args.name:
        vcf_reader = vcf.Reader(open(name, 'r'))
        parse = Parser(vcf_reader, name, direct)
        parse.get_arguments()

        if args.population:
            csv_file = direct + '/' + 'populations_'  + name.replace('.vcf','.csv')
            f = open(csv_file, 'w')
            population_writer = csv.DictWriter(f, fieldnames= fieldnames, delimiter='\t')
            population_writer.writeheader()
            f.close()

        parse.generate_filtered_vcf()

    os.system(f'chmod -R g+w {direct}') #gives  read write to directory for group





