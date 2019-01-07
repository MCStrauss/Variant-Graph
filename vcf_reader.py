import vcf,  vcf.filters
from collections import deque, defaultdict
import argparse
from frequency import Frequency
#parser = argparse.ArgumentParser(description = 'filter for VCF file')
#parser.add_argument('siteQ', type = int, help = 'integer value to filter bad sites')

vcf_reader = vcf.Reader(open('vcf.txt', 'r'))

alt_freq = {} #holds sample: namedtuple Allfreq

history  = deque(maxlen = 2)

def main():
    check = 0
    DP = AD = 0
    for record in vcf_reader:
        key = (record.CHROM, record.POS)
        alt_freq[key] = {}
        for call in record.samples:

            if call.sample[0].isdigit():
                sample = call.sample[:2]
            else:
                sample = call.sample[0]
            history.append(sample)

            if sample != history[0]:
                DP = AD = 0
            check += 1
            try:
                DP += call['DP']
                AD += int(call['AD'][0])
                #print(key, sample, DP, AD)
                if DP and AD:
                    if sample not in alt_freq[key]:
                        alt_freq[key][sample] = Frequency(DP, AD, (DP - AD) / DP)
                    else:
                        alt_freq[key][sample].update(DP, AD)
                        #print(key, sample, alt_freq[key][sample].DP, alt_freq[key][sample].AD)
            except TypeError: #avoids getting edgecase where DP = None and thus adding int + Nonetype throws error
                DP = AD = 0
            #if check>8:
             #   break


def check():
    vcf_reader = vcf.Reader(open('vcf.txt', 'r'))
    with open('reference.txt', 'w') as check:
        for record in vcf_reader:
            check.write(f'{record.CHROM} {record.POS} {record.aaf}\n')

    with open('our_data.txt', 'w') as data:
        for chrom, samples in alt_freq.items():
            DP = AD = 0
            for val in samples.values():
                DP += val.DP
                AD += val.AD

            data.write(f'{chrom} {(DP - AD) / DP}\n')



if __name__ == '__main__':
    main()
    with open('output.txt', 'w') as out:
        for key, value in alt_freq.items():
            out.write(f'sample_name - {key}')
            for val in value.items():
                out.write(f'{val}'.center(10,' '))
            out.write('\n')
    check()



'''
 history.append(sample)

            if sample != history[0]:
                try:
                    if sample not in alt_freq[(record.CHROM, record.POS)]:
                        alt_freq[(record.CHROM, record.POS)][sample] = Frequency(DP, AD, (DP - AD) / DP)
                    else:
                        alt_freq[(record.CHROM, record.POS)][sample].update(DP, AD)

                except ZeroDivisionError:
                    pass
                    # alt_freq[sample] = Frequency(DP, AD, 0)

                DP = AD = 0

'''