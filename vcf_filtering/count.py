#Frequency object will be used to store the reference count, alternative count, and the minor allelic frequency
class Frequency:

    def __init__(self,  ref = 0, alt = 0):
        self.ref = ref #number of reads in support of reference
        self.alt = alt #number of reads in support of alternative
        self.get_freq()


    def update(self, ref = 0, alt = 0):
        self.ref += ref
        self.alt += alt
        self.get_freq()

    def get_freq(self):
        try:
            self.minor_freq = self.alt / (self.ref + self.alt)
        except ZeroDivisionError:
            self.minor_freq = 0 # local minor allelic frequency
    def __repr__(self):
        #for debugging purposes
        return f'{self.__class__.__name__}(ref {self.ref}, alt {self.alt}, maf = {self.minor_freq}\t'

    def __str__(self):
        #for looking pretty and making it easier to parse
        return f'ref = {self.ref}\talt = {self.alt}\tmaf = {self.minor_freq}\t'