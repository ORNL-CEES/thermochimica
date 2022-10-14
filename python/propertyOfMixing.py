import thermoTools
import subprocess
import json
import shutil

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

    # Check if heat capacity should be turned on
    heatCapacity = True if property in ["enthalpy","entropy"] else False

    elements = set()
    for endpoint in endpoints:
        elements.update(endpoint.keys())
    # If an element is not in an endpoint, add it (0 concentration)
    for endpoint in endpoints:
        for element in elements:
            if element not in endpoint.keys():
                endpoint[element] = 0

    # Write and run endpoint calculations
    inputFileName = f'{thermochimica_path}/inputs/{property.replace(" ","_")}OfMixing-{phase}-{temperature}{tunit}-endpoints.ti'
    # Set up endpoint calculations
    calcList = []
    for endpoint in endpoints:
        calc = [temperature, pressure]
        calc.extend([endpoint[element] for element in elements])
        calcList.append(calc)
    thermoTools.WriteRunCalculationList(inputFileName,database,elements,calcList,tunit=tunit,punit=punit,munit=munit,printMode=0,heatCapacity=heatCapacity,excludePhasesExcept=[phase])

    # Run calculation
    thermoTools.RunRunCalculationList(inputFileName)

    # Process output
    f = open(thermochimica_path+'/outputs/thermoout.json',)
    try:
        data = json.load(f)
        f.close()
    except:
        f.close()
        print('Data load failed, aborting mixture calculation')
        return

    endpointProp = []
    i = 0
    for endpoint in endpoints:
        i += 1
        endpointProp.append(data[str(i)][property])

    # Write and run mixture calculations
    inputFileName = f'{thermochimica_path}/inputs/{property.replace(" ","_")}OfMixing-{phase}-{temperature}{tunit}.ti'
    # Set up mixture calculations
    calcList = []
    for mixture in mixtures:
        calc = [temperature, pressure]
        calc.extend([(1-mixture)*endpoints[0][element] + (mixture)*endpoints[1][element] for element in elements])
        calcList.append(calc)
    thermoTools.WriteRunCalculationList(inputFileName,database,elements,calcList,tunit=tunit,punit=punit,munit=munit,printMode=0,heatCapacity=heatCapacity,excludePhasesExcept=[phase])

    # Run calculation
    thermoTools.RunRunCalculationList(inputFileName)

    # Process output
    f = open(thermochimica_path+'/outputs/thermoout.json',)
    try:
        data = json.load(f)
        f.close()
    except:
        f.close()
        print('Data load failed, aborting mixture calculation')
        return

    mixtureProp = []
    i = 0
    for mixture in mixtures:
        i += 1
        prop = data[str(i)][property] - ((1-mixture)*endpointProp[0] + (mixture)*endpointProp[1])
        mixtureProp.append(prop)
        data[str(i)][f'{property} of mixing'] = prop

    # Save data
    with open('outputs/thermoout.json', 'w') as outfile:
        json.dump(data, outfile, indent=4)
    shutil.copy2('outputs/thermoout.json', f'outputs/{property.replace(" ","_")}OfMixing-{phase}-{temperature}{tunit}.json')

    return mixtureProp
