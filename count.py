class Frequency:

    def __init__(self, ref = 0, alt = 0):
        self.ref = ref
        self.alt = alt
        self.get_freq()

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
        return f'{self.__class__.__name__}(reference {self.ref}, alternative {self.alt}, minor allele frequency = {self.freq})'