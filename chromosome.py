class chromosomeValue:
    num_instances = 0 #keeps track of the number of instances of our class
    def __init__ (self):
        self.Rcount = 0  #reference count
        self.Mcount = 0  #Mutation count
        self.__class__.num_instances += 1

    def __str__(self):
        return 'reference_count: {}, mutation_count: {}'.format(self.Rcount, self.Mcount)

    def __iter__(self):
        return iter(( self.Rcount, self.Mcount))

    def __getitem__(self, item):
        return (self.Rcount, self.Mcount)[item]

    def find_count(self, c1 ,c2):
        '''
        GT genotype, encoded as alleles values separated by
        either of ”/” or “|”, e.g. The allele values are 0 for the reference allele (what is in the reference sequence),
        1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on.
        :param c1: 0,1 or 2
        :param c2: The count of a mutation or a reference.
        :return:
        '''
        if c1 and c2:
            self.Mcount += 1 #increment mutation count by 2
        elif not c1 and not c2:
            self.Rcount += 1 #increment reference count by 2
        else:
            self.Rcount += 1 #increment reference count by 1
            self.Mcount += 1 #increment mutation count by 1

    def total(self):
        return self.Mcount/(self.Rcount+self.Mcount)