import thermoTools
import subprocess
import json

def propertyOfMixing(property, phase, temperature, endpoints, mixtures, database,
                     thermochimica_path = 'bin/RunCalculationList',
                     pressure = 1,
                     tunit = 'K',
                     punit = 'atm',
                     munit = 'moles'):
    # Check that property is one of the recognized options
    propList = ['integral Gibbs energy','enthalpy','entropy']
    if property not in propList:
        print(f'Property {property} not recognized, please use one of {", ".join(propList)}')

    elements = set()
    for endpoint in endpoints:
        elements.update(endpoint.keys())
    # If an element is not in an endpoint, add it (0 concentration)
    for endpoint in endpoints:
        for element in elements:
            if element not in endpoint.keys():
                endpoint[element] = 0

    inputFileName = f'inputs/{property.replace(" ","_")}OfMixing-{phase}-{temperature}{tunit}.ti'
    with open(inputFileName, 'w') as inputFile:
        inputFile.write(f'data file              = {database}\n')
        inputFile.write(f'temperature unit       = {tunit}\n')
        inputFile.write(f'pressure unit          = {punit}\n')
        inputFile.write(f'mass unit              = {munit}\n')
        inputFile.write(f'nEl                    = {len(elements)} \n')
        inputFile.write(f'iEl                    = {" ".join([str(thermoTools.atomic_number_map.index(element)+1) for element in elements])}\n')
        inputFile.write(f'number excluded except = 1\n')
        inputFile.write(f'phases excluded except = {phase}\n')
        inputFile.write(f'heat capacity          = {".TRUE." if property in ["enthalpy","entropy"] else ".FALSE."}\n')
        inputFile.write(f'nCalc                  = {len(endpoints) + len(mixtures)}\n')
        # Write endpoint calculations
        for endpoint in endpoints:
            inputFile.write(f'{temperature} {pressure} {" ".join([str(endpoint[element]) for element in elements])}\n')
        # Write mixture calculations
        for mixture in mixtures:
            inputFile.write(f'{temperature} {pressure} {" ".join([str((1-mixture)*endpoints[0][element] + (mixture)*endpoints[1][element]) for element in elements])}\n')
    # Run calculation
    subprocess.run([thermochimica_path,inputFileName])

    # Process output
    f = open('thermoout.json',)
    try:
        data = json.load(f)
        f.close()
    except:
        f.close()
        print('Data load failed, aborting phase diagram update')
        return

    endpointProp = []
    i = 0
    for endpoint in endpoints:
        i += 1
        endpointProp.append(data[str(i)][property])

    print()
    mixtureProp = []
    for mixture in mixtures:
        i += 1
        prop = data[str(i)][property] - ((1-mixture)*endpointProp[0] + (mixture)*endpointProp[1])
        mixtureProp.append(prop)
        print(f'{mixture} {prop}')


property = 'entropy'
phase = 'MSFL'
temperature = 1000
tunit = 'C'
# endpoints = [dict([('Li',1),('F',1)]),dict([('Pu',1),('F',3)])]
endpoints = [dict([('Na',1),('F',1)]),dict([('Ca',1),('F',2)])]
mixtures = [i/10 for i in range(1,10)]
database = 'mstdb/Models and Documentation/MSTDB-TC_V1.3_Fluorides_8-0.dat'

propertyOfMixing(property, phase, temperature, endpoints, mixtures, database, tunit = tunit)
