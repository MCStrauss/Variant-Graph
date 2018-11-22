from itertools import dropwhile
db={}
class chromosomeValue:
    Rcount=Mcount=0
    def init(self):
        self.Rcount=0  #reference count
        self.Mcount=0  #Mutation count

    def __str__(self):
        return 'reference_count: {}, mutation_count: {}'.format(self.Rcount, self.Mcount)

    def __iter__(self):
        return iter(self.Rcount,self.Mcount)

    def __getitem__(self, item):
        return (self.Rcount,self.Mcount)[item]

    def find_count(self,c1,c2):
        '''
        GT genotype, encoded as alleles values separated by
        either of ”/” or “|”, e.g. The allele values are 0 for the reference allele (what is in the reference sequence),
        1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on.
        :param c1: 0,1 or 2
        :param c2: The count of a mutation or a reference.
        :return:
        '''
        if c1 and c2:
            print('case1',c1,c2,sep=',')
            self.Mcount+=2
        elif not c1 and not c2:
            print('case2',c1, c2, sep=',')
            self.Rcount += 2
        else:
            print('case3',c1, c2, sep=',')
            self.Rcount+=1
            self.Mcount+=1
    def total(self):
        #print(self.Mcount,self.Rcount,sep=',')
        return 'Mutation rate = {}'.format(self.Mcount/(self.Rcount+self.Mcount))

def main():
    with open('vcf.txt') as f:
        for line in dropwhile(lambda x: x.startswith('#'),f):
            split=line.split()
            name, position=split[0], split[1]
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
                db[key].find_count(c1, c2)

if __name__=='__main__':
    main()
    for item in db.keys():
        print(item,db[item].total(),sep=': ')
        



