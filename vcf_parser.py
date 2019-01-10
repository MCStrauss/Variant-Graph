import vcf
import argparse
from collections import namedtuple, defaultdict
from count import Frequency
from populations import populations
chromosome = namedtuple('chromosome', ['chrom', 'position'])

parser = argparse.ArgumentParser(description  = 'Arguments for VCF filter script')
parser.add_argument('-an', '--an_cutoff', help = 'AN count must represent certain %% of total', type = float)

parser.add_argument('-af', '--allelic_frequency', help = '%% to filter for local allelic frequency. '
                    'If a read has no sample with greater than x%% then '
                    'it is removed \n from the vcf. ', type = float)

parser.add_argument('-p', '--population', help = 'If -p is specified an outpout file with population and '
                                                 'their minor allelic frequencies is outputed', action = 'store_true')
args = parser.parse_args()

def check():
    if args.an_cutoff: assert 0 <= args.an_cutoff < 1, f'an_cutoff must be atleast 0 and greater than 1'
    if args.allelic_frequency: assert 0 <= args.allelic_frequency < 1, f'allelic_frequency must be atleast 0 and greater than 1'

class Parser:

    def __init__(self, vcf_reader, f_aan = .5, f_maf = .05):
        self.vcf_reader = vcf_reader
        self.dB = defaultdict(dict)
        self.f_aan = f_aan
        self.f_maf = f_maf

    @staticmethod
    def parse(data):
        left, right = data.split('/')
        if left == '.' or right == '.':
            return (0, 0)

        left, right = eval(left), eval(right)

        if left or right:
            return (0, 1)
        else:
            return (1, 0)

    def gen_pop_data(self):
        for record in vcf_reader:
            chrom = chromosome(record.CHROM, record.POS)

            for sample in record.samples:
                x, y = self.parse(sample['GT'])

                for func, key in populations.items():
                    if func(sample.sample):
                        if chrom in self.dB[key]:
                            self.dB[key][chrom].update(ref = x, alt = y) #already in the dB
                            break
                        else:
                            self.dB[key][chrom] = Frequency(ref = x, alt = y) #needs to be initialized in the dB
                            break

    def write_pop_data(self, name = 'check_break.txt'):
        with open(name, 'w') as out:
            for key, value in self.dB.items():
                out.write(f'{key}, {value}\n')

    def filter_ANN(self, line):
        #todo insert f_aan here and also figure out percentages
        an = eval(line.split(';')[2].split('=')[1])  # nasty
        if an / 2 < 535:
            return False
        return True

    def filter_local_frequency(self,chrom):
        for population in self.dB:
            if self.dB[population][chrom].freq >= self.f_maf:
                return True
        return False

    def generate_filtered_vcf(self, name = 'new_vcf.txt'):
        '''

        :param name: name of the new vcf.txt
        :return: will generate a filtered vcf if it fits aan parameters and atleast 1 chromosome has an minor allelic
        frequency above .05% or the argument the user specified, relies on two helpers function filter_local_frequency and filter_AAN.
        '''
        with open('vcf.txt') as input, open(name, 'w') as out:
            for line in input:
                if line.startswith('#'):
                    out.write(f'{line}')
                else:
                    split = line.split('\t')
                    if self.filter_ANN(split[7]) and self.filter_local_frequency((split[0],eval(split[1]))): #gets chrom name as string and position as int
                        out.write(f'{line}')

if __name__ == '__main__':
    check() #checks to make sure args are valid if there are args
    vcf_reader = vcf.Reader(open('vcf.txt', 'r'))
    parse = Parser(vcf_reader, args.an_cutoff if args.an_cutoff else .5, args.allelic_frequency if args.allelic_frequency else .05)
    parse.gen_pop_data()

    if args.population: parse.write_pop_data()
    parse.generate_filtered_vcf()