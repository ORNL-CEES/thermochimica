import thermoTools
import subprocess
import json

def propertyOfMixing(property, phase, temperature, endpoints, mixtures, database,
                     thermochimica_path = '.',
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

    # Write and run endpoint calculations
    inputFileName = f'{thermochimica_path}/inputs/{property.replace(" ","_")}OfMixing-{phase}-{temperature}{tunit}.ti'
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
        inputFile.write(f'nCalc                  = {len(endpoints)}\n')
        # Write endpoint calculations
        for endpoint in endpoints:
            inputFile.write(f'{temperature} {pressure} {" ".join([str(endpoint[element]) for element in elements])}\n')

    # Run calculation
    subprocess.check_output([thermochimica_path+'/bin/RunCalculationList',inputFileName])

    # Process output
    f = open(thermochimica_path+'/thermoout.json',)
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

    # Write and run mixture calculations
    inputFileName = f'{thermochimica_path}/inputs/{property.replace(" ","_")}OfMixing-{phase}-{temperature}{tunit}.ti'
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
        inputFile.write(f'nCalc                  = {len(mixtures)}\n')
        # Write mixture calculations
        for mixture in mixtures:
            inputFile.write(f'{temperature} {pressure} {" ".join([str((1-mixture)*endpoints[0][element] + (mixture)*endpoints[1][element]) for element in elements])}\n')
    # Run calculation
    subprocess.check_output([thermochimica_path+'/bin/RunCalculationList',inputFileName])

    # Process output
    f = open(thermochimica_path+'/thermoout.json',)
    try:
        data = json.load(f)
        f.close()
    except:
        f.close()
        print('Data load failed, aborting phase diagram update')
        return

    print()
    mixtureProp = []
    i = 0
    for mixture in mixtures:
        i += 1
        prop = data[str(i)][property] - ((1-mixture)*endpointProp[0] + (mixture)*endpointProp[1])
        mixtureProp.append(prop)
        print(f'{mixture} {prop}')
        data[str(i)][f'{property} of mixing'] = prop

    # Save data
    with open('thermoout.json', 'w') as outfile:
        json.dump(data, outfile, indent=4)
    subprocess.run(['cp','thermoout.json',f'{property.replace(" ","_")}OfMixing-{phase}-{temperature}{tunit}.json'])

    return mixtureProp
