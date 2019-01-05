import vcf,  vcf.filters
from collections import namedtuple, deque
import argparse

#parser = argparse.ArgumentParser(description = 'filter for VCF file')
#parser.add_argument('siteQ', type = int, help = 'integer value to filter bad sites')

vcf_reader = vcf.Reader(open('vcf.txt', 'r'))

Allfreq = namedtuple('Allfreq',['DP', 'AD', 'AF'])

alt_freq = {} #holds sample: namedtuple Allfreq

history  = deque(maxlen = 2)

DP = AD = 0


for record in vcf_reader:
    for call in record.samples:
        if call.sample[0].isdigit():
            sample = call.sample[:2]
        else:
            sample = call.sample[0]

        history.append(sample)
        if sample != history[0]:
            try:
                alt_freq[sample] = Allfreq(DP, AD, (DP - AD) / DP)
            except ZeroDivisionError:
                alt_freq[sample] = Allfreq(DP, AD, 0)
            DP = AD = 0

        ifDP = call['DP']

        if ifDP:
            DP += ifDP
            AD += int(call['AD'][0])

with open('output.txt','w') as out:
    for key, value in alt_freq.items():
        out.write(f'sample_name - {key}:  {value}\n')


def main():
    pass



if __name__ == '__main__':
    main()


