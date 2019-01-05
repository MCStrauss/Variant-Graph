import vcf,  vcf.filters
from collections import deque
import argparse
from frequency import Frequency
#parser = argparse.ArgumentParser(description = 'filter for VCF file')
#parser.add_argument('siteQ', type = int, help = 'integer value to filter bad sites')

vcf_reader = vcf.Reader(open('vcf.txt', 'r'))

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
                if sample not in alt_freq:
                    alt_freq[sample] = Frequency(DP, AD, (DP - AD) / DP)
                else:
                    alt_freq[sample].update(DP,AD)

            except ZeroDivisionError:
                pass
                #alt_freq[sample] = Frequency(DP, AD, 0)

            DP = AD = 0

        ifDP = call['DP']

        if ifDP:
            DP += ifDP
            AD += int(call['AD'][0])

with open('other.txt','w') as out:
    for key, value in alt_freq.items():
        out.write(f'sample_name - {key}:  {value}\n')


def main():
    pass

if __name__ == '__main__':
    main()


