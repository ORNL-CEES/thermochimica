import propertyOfMixing

property = 'entropy'
phase = 'BCCN'
temperature = 1000
tunit = 'C'
endpoints = [{'Pd':1},{'Mo':1}]
ncalc = 10
mixtures = [i/ncalc for i in range(1,ncalc)]
database = 'data/Kaye_NobleMetals.dat'

values = propertyOfMixing.propertyOfMixing(property, phase, temperature, endpoints, mixtures, database, tunit = tunit)

for i in range(len(mixtures)):
    print(f'{mixtures[i]} {values[i]}')
