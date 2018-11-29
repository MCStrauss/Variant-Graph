from itertools import dropwhile
db={}
cc_count=0

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
            self.Mcount+=1 #increment mutation count by 2
        elif not c1 and not c2:
            self.Rcount += 1 #increment reference count by 2
        else:
            self.Rcount+=1 #increment reference count by 1
            self.Mcount+=1 #increment mutation count by 1
    def total(self):
        return 'Mutation rate = {}'.format(self.Mcount/(self.Rcount+self.Mcount))

def main():

    with open('vcf.txt') as f:
        for line in dropwhile(lambda x: x.startswith('#'),f):
            split=line.split()
            global cc_count
            if split[2]=='.': cc_count+=1
            name, position=split[0], split[1]
            key = (name, position)

            if key not in db.keys():
                value=chromosomeValue()
                db[key]=value

                print(cc_count)
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




