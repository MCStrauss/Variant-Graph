import vcf
import argparse
from collections import namedtuple, defaultdict
from count import Frequency
from populations import populations
#filter on fs
chromosome = namedtuple('chromosome', ['chrom', 'position'])

parser = argparse.ArgumentParser(description  = 'Arguments for VCF filter script')
parser.add_argument('-an', '--an_cutoff', help = 'AN count must represent certain %% of total', type = float)

parser.add_argument('-af', '--allelic_frequency', help = '%% to filter for local allelic frequency. '
                    'If a read has no sample with greater than x%% then '
                    'it is removed \n from the vcf. ', type = float)

parser.add_argument('-p', '--population', help = 'If -p is specified an output file with population and '
                                                 'their minor allelic frequencies is outputed', action = 'store_true')
parser.add_argument('-g', '--global_frequency', help = 'gives the cutoff for the global frequency of a position on the chromosome,'
                                                       'if it is less than the cutoff it will not be included in the filtered vcf')


args = parser.parse_args()

def check():
    if args.an_cutoff: assert 0 <= args.an_cutoff < 1, 'an_cutoff must be atleast 0 and greater than 1'
    if args.allelic_frequency: assert 0 <= args.allelic_frequency < 1, 'allelic_frequency must be atleast 0 and greater than 1'

class Parser:
    def __init__(self, vcf_reader, f_aan = .5, f_maf = .05):
        '''
        :param vcf_reader: the vcf file we will be iterating over
        :param f_aan: cutoff for mininum number of alleles
        :param f_maf: cutoff for minor allelic frequency
        '''
        self.vcf_reader = vcf_reader
        self.dB = defaultdict(dict) #our nested dictionary, it will be represented  as follows
        # key = population: val = (chrom_name, positiont): val = Frequency(#alleles_refernce, #alleles_alternate, #minor allelic frequency)
        self.f_aan = f_aan
        self.f_maf = f_maf

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

    def gen_pop_data(self):
        '''
        :return: returns nothing, this generates our dB,
        '''
        for record in vcf_reader:

            chrom = chromosome(record.CHROM, record.POS)
            aaf = max(record.aaf)

            for sample in record.samples:
                x, y = self.parse(sample['GT'])

                for func, key in populations.items(): #populations is a dictionary of functions and values where the values are our samples and the functions determine which sample corresponds to the sample header from the vcf.
                    if func(sample.sample):
                        if chrom in self.dB[key]:
                            self.dB[key][chrom].update(ref = x, alt = y) #already in the dB
                            break
                        else:
                            self.dB[key][chrom] = Frequency(aaf, ref = x, alt = y) #needs to be initialized in the dB
                            break
        return

    def write_pop_data(self, name = 'data_populations.txt'):
        '''
        :param name: name for the text file the user wants the population data in
        :return: A text file that displays the population: and the minor alellic frequency at each chromosome.
        '''
        with open(name, 'w') as out:
            for key, value in self.dB.items():
                out.write(f'{key}\t')
                for chrom, freq in value.items():
                    out.write(f'{chrom.chrom}\t{chrom.position}\t{freq}\t' )
                out.write('\n')


    def filter_AN_FS(self, line):
        '''
        :param line: a row in a vcf file
        :return:  if the an doesn't have the desired number of samples we want remove the row from our new vcf file.
        '''
        #todo filter on fs
        fs = eval(''.join([x for x in line.split(';')[5:7] if x.startswith('FS')]).split('=')[1]) #grabs fischer strand
        an = eval(line.split(';')[2].split('=')[1])  # grabs the an from the line

        if (an / 2) / 1070 < self.f_aan: #1070 is the total number of samples, an represents the number of chromosomes in the record
            return False
        return True

    def filter_frequency(self, chrom):
        '''
        :param chrom: given the (chrom, position) check everywhere it occured and its minor allelic frequency.
        :return: If it has a minor allelic frequency of atleast our cutoff keep it in, else remove the line from our filtered VCF.
        '''
        for population in self.dB:
            if chrom in self.dB[population]: # haven't seen this yet but its here to prevent edge case where the a population
                # doesn't have a read on a chromosome, remove this line and we will get a KeyError otherwise
                loc_freq, glob_freq = self.dB[population][chrom].freq, self.dB[population][chrom].aaf
                if loc_freq >= self.f_maf or glob_freq >= self.f_maf :
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
                    split = line.split('\t') #fischer strand
                    chrom = split[0], eval(split[1])

                    if self.filter_AN_FS(split[7]) and self.filter_frequency(chrom): #gets chrom name as string and position as int
                        out.write(f'{line}')

if __name__ == '__main__':
    check() #checks to make sure args are valid if there are args
    vcf_reader = vcf.Reader(open('vcf.txt', 'r'))
    parse = Parser(vcf_reader, args.an_cutoff if args.an_cutoff else .5, args.allelic_frequency if args.allelic_frequency else .05)
    parse.gen_pop_data()
    parse.generate_filtered_vcf()
    if args.population: parse.write_pop_data()

    '''
    Analysis of runtime: generating our database based on populations takes O(n * m * 18) where 
    n = number of row and m = number of columns. The constant 18 comes from populations.py, at worst to find the populations
    our for loop will have to execute the entire way through (however this will almost never happen).
    
    Generating the new vcf requires us to make another pass through the original vcf and compare it to the data we generated. So we must again 
    go through each line of the vcf (n) and then walk through every chromosome in our database by population (m). So another n*m.
    
    Our approximate runtime is therefore O(2n * 2m) = O(n*m).
    
    Storage:
    Since we are using generators to iterate over the file our memory is being consumed by our database which hold every population and chromosome
    so O(n*m) storage. When generating the new vcf we also take into one line to check the an value however it is constant so O(1) storage.
     
    '''