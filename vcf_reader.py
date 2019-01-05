import vcf,  vcf.filters
from itertools import takewhile
from collections import namedtuple, deque


vcf_reader = vcf.Reader(open('vcf.txt', 'r'))

Allfreq = namedtuple('Allfreq',['DP', 'AD', 'Complement_AD'])

alt_freq = {}



#divide minor_allele_freq =  (AD[0] for AD[1])   / DP

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
            alt_freq[sample] = Allfreq(DP, AD, DP - AD)
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



    #print(record, record.aaf, record.call_rate)


#for i, sample in enumerate(record.samples):
    #print(i, sample)

#for i, item in enumerate(vcf_reader):
    #print(i, item)
#print(vcf_reader.samples)








#test_output_2000001-3000000.vcf
#parser = argparse.ArgumentParser(description = 'filtering VCF files for a certain quality')
#parser.add_argument('siteQ', type = int, help = 'filter sites below this quality')
#args = parser.parse_args()
#vcf.filters.SiteQuality.customize_parser(parser)
#vcf_filter.py --local-script vcf_reader.py

'''
def place_in_dict(sample):
    split = sample['GT'].split('/')
    if len(split) < 2:
        split = sample['GT'].split('|')
    zero = split.count('0')
    one  = split.count('1')
    dot = split.count('.')
    if dot:
        #this means we need to throw it out
        return
    elif zero and not one:
        alt_freq[sample[:2]].Reference += 1
    elif one and not zero:
        alt_freq[sample[:2]].Alternative += 1
'''
