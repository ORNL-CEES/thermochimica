import subprocess

filename = 'inputs/pythonInput.ti'

temperature = 900
pressure = 2
elements = [6,8]
masses = [1,1.7]
tunit = 'K'
punit = 'atm'
munit = 'moles'
datafile = 'data/C-O.dat'

with open(filename, 'w') as inputFile:
    inputFile.write('! Python-generated input file for Thermochimica\n')
    inputFile.write('temperature          = ' + str(temperature) + '\n')
    inputFile.write('pressure          = ' + str(pressure) + '\n')
    for i in range(len(masses)):
        inputFile.write('mass(' + str(elements[i]) + ')           = ' + str(masses[i]) + '\n')
    inputFile.write('temperature unit         = ' + tunit + '\n')
    inputFile.write('pressure unit          = ' + punit + '\n')
    inputFile.write('mass unit          = ' + munit + '\n')
    inputFile.write('data file         = ' + datafile + '\n')
    inputFile.write('print mode        = 2\n')
    inputFile.write('debug mode        = .FALSE.\n')

subprocess.call(['./bin/ThermochimicaInputScriptMode',filename])
