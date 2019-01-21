#Frequency object will be used to store the reference count, alternative count, and the minor allelic frequency
class Frequency:

    def __init__(self, aaf, ref = 0, alt = 0):
        self.ref = ref #number of reads in support of reference
        self.alt = alt #number of reads in support of alternative
        self.get_freq()
        self.aaf = aaf  #aaf = alternate allele frequency


    def update(self, ref = 0, alt = 0):
        self.ref += ref
        self.alt += alt
        self.get_freq()

    def get_freq(self):
        try:
            self.freq = self.alt / (self.ref + self.alt)
        except ZeroDivisionError:
            self.freq = 0

    def __repr__(self):
        return f'{self.__class__.__name__}(ref {self.ref}, alt {self.alt}, maf = {self.freq}, gaf {self.aaf})'