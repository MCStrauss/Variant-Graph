from itertools import dropwhile
from chromosome import chromosomeValue
import argparse
db={} #storage (chro. name, position): reference count, mutation count stored as an object in chromosome.py
#accessed as val=chromsomeValue(): val.Mcount or val.Rcount

parser=argparse.ArgumentParser(description = 'Find the allele frequencies of '
                                             'drosophilia flies from vcf files '
                                             'Example run: python3 vcf_parser.py .5')

parser.add_argument('filter', type = float, metavar='', help = 'filter: Cutoff to see if a variant is '
                                                               'significant or not Please supply filter as a float'
                                                               'for instance .5')
args = parser.parse_args()
assert 0 <= args.filter <= 1,  ('Filter value should not be below 0 or exceed 1')


def main():
    with open('vcf.txt') as f:
        for line in dropwhile(lambda x: x.startswith('#'),f):
            split=line.split()
            name, position = split[0], split[1]
            key = (name, position)

            if key not in db.keys():
                value=chromosomeValue()
                db[key]=value

            for sample in split[9:]:
                if len(sample.split('|'))<2:
                    sample=sample.split('/')
                else:
                    sample=sample.split('|')
                c1,c2=eval(sample[0]), eval(sample[1][0])
                db[key].find_count(c1, c2) #accumulating the count of reference and mutation alleles

            if db[key].total()<args.filter: #filtering step
                del db[key]

if __name__=='__main__':
    main()
    for item in db.keys():
        print(item,db[item].total(),sep=': ')





