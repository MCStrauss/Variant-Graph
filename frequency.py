class Frequency:
    def __init__(self, DP, AD, AF):
        self.DP = DP
        self.AD = AD
        self.AF = AF
    def __str__(self):
        return f'DP: {self.DP} AD: {self.AD} AF: {self.AF}'

    def __repr__(self):
        return f'{self.__class__.__name__}({self.DP}, {self.AD}, {self.AF})'

    def update(self, DP, AD):
        self.DP += DP
        self.AD += AD
        self.AF = (self.DP - self.AD) / self.DP


