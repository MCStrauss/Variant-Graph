#dictionary  used to support generating our population data
#todo double check which samples we want grouped together and make sure that every sample has a function
#exclude Cairo 'EG'
#include GA under georgia North America
#Gabon goes to Africa
#include lower_t and lower_w in global but exclude from populations
#uganda Tanzania <- Africa
#Zimbabwe is other Africa
#make Europe a population group

populations = { lambda x: x[0].isdigit(): 'North America',
                lambda x: x.startswith('B'): 'China',
                lambda x: x.startswith('C'): 'Other Africa',
                lambda x: x.startswith('E'): 'Ethiopa',
                lambda x: x.startswith('F'): 'Europe',
                lambda x: x.startswith('G'): 'North America',
                lambda x: x.startswith('H'): 'North America',
                lambda x: x.startswith('I'): 'North America',
                lambda x: x.startswith('K'): 'Other Africa',
                lambda x: x.startswith('N'): 'Europe',
                lambda x: x.startswith('RAL'): 'North America',
                lambda x: x.startswith('RC'): 'Other Africa',
                lambda x: x.startswith('RA'): 'Other Africa',
                lambda x: x.startswith('S'): 'Southern Africa',
                #lambda x: x.startswith('T'): 'Australia',
                lambda x: x.startswith('W'): 'North America',
                lambda x: x.startswith('ZH') : 'Other Africa',
                lambda x: x.startswith('ZK'): 'Other Africa',
                lambda x: x.startswith('ZS'): 'Other Africa',
                lambda x: x.startswith('ZW'): 'Other Africa',
                lambda x: x.startswith('ZI') : 'Southern Africa',
                lambda x: x.startswith('ZL'): 'Southern Africa',
                lambda x: x.startswith('ZO'): 'Southern Africa',
                lambda x: x.startswith('b'): 'China',
                #lambda x: x.startswith('t'): 'lower_t',
                #lambda x: x.startswith('w'): 'lower_w',
}

