import PySimpleGUI as sg
import subprocess
import math
import os

atomic_number_map = [
    'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P',
    'S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
    'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh',
    'Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd',
    'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re',
    'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
    'Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db',
    'Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts', 'Og'
]

file_list_column = [
    [
        sg.Text("Database Folder"),
        sg.In(size=(25, 1), enable_events=True, key="-FOLDER-"),
        sg.FolderBrowse(),
    ],
    [
        sg.Listbox(
            values=[], enable_events=True, size=(40, 20), key="-FILE LIST-"
        )
    ],
]

dataWindow = sg.Window('Thermochimica database selection', file_list_column,finalize=True)
folder = os.getcwd()+'/data'
# print(folder)
try:
    file_list = os.listdir(folder)
except:
    file_list = []

fnames = [
    f
    for f in file_list
    if os.path.isfile(os.path.join(folder, f))
    and f.lower().endswith((".dat", ".DAT"))
]
dataWindow["-FILE LIST-"].update(fnames)

timeout = 50

while True:
    event, values = dataWindow.read()
    # print(event, values)
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    elif event == "-FOLDER-":
        folder = values["-FOLDER-"]
        try:
            file_list = os.listdir(folder)
        except:
            file_list = []

        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(folder, f))
            and f.lower().endswith((".dat", ".DAT"))
        ]
        dataWindow["-FILE LIST-"].update(fnames)
    elif event == "-FILE LIST-":  # A file was chosen from the listbox
        try:
            datafile = os.path.join(
                folder, values["-FILE LIST-"][0]
            )
            with open(datafile) as f:
                f.readline() # read comment line
                line = f.readline() # read first data line (# elements, # phases, n*# species)
                nElements = int(line[1:5])
                nSoln = int(line[6:10])
                elements = ['']*nElements
                elIt = 0
                for i in range(math.ceil((nSoln+3)/15)-1):
                    f.readline() # read the rest of the # species but don't need them)
                for i in range(math.ceil(nElements/3)):
                    els = f.readline() # read a line of elements (3 per line)
                    elLen = 25 # formatted 25 wide
                    for j in range(3):
                        elements[elIt] = els[1+j*elLen:(1+j)*elLen].strip()
                        elIt = elIt + 1
                        if elIt == nElements:
                            break # breaks inner
                    else:
                        continue # if inner not broken, continue inner
                    break # breaks outer if inner broken
            i = 0
            while i < nElements:
                try:
                    index = atomic_number_map.index(elements[i])+1 # get element indices in PT (i.e. # of protons)
                    # print(elements[i],index)
                    i = i + 1
                except ValueError:
                    print(elements[i]+' not in list') # if the name is bogus (or e(phase)), discard
                    elements.remove(elements[i])
                    # print(elements)
                    nElements = nElements - 1
            tempLayout = [[sg.Text('Temperature')],[sg.Input(key='-temperature-',size=(16,1))],
                          [sg.Text('End Temperature',key='-endtemperaturelabel-')],[sg.Input(key='-endtemperature-',size=(16,1))],
                          [sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstep-',size=(8,1))],
                          [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')],
                          [sg.Text('Temperature range:')],
                          [sg.Radio('Disabled', 'trange', default=True, enable_events=True, key='-tdis-')],
                          [sg.Radio('Enabled', 'trange', default=False, enable_events=True, key='-ten-')]]
            presLayout = [[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(16,1))],
                          [sg.Text('End Pressure',key='-endpressurelabel-')],[sg.Input(key='-endpressure-',size=(16,1))],
                          [sg.Text('# of steps',key='-psteplabel-')],[sg.Input(key='-pstep-',size=(8,1))],
                          [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')],
                          [sg.Radio('Disabled', 'prange', default=True, enable_events=True, key='-pdis-')],
                          [sg.Radio('Enabled', 'prange', default=False, enable_events=True, key='-pen-')],
                          [sg.Radio('Enabled, step with temperature', 'prange', default=False, enable_events=True, key='-pent-')]]
            elemLayout = []
            for i in range(nElements):
                elemLayout.append([sg.Text(elements[i])])
                elemLayout.append([sg.Input(key='-'+elements[i]+'-',size=(16,1))])
            calcLayout = [tempLayout,
                          presLayout,
                          elemLayout,
                          [sg.Text('Mass unit')],
                          [sg.Combo(['moles', 'kg', 'atoms', 'g'],default_value='moles',key='-munit-')],
                          [sg.Button('Run'), sg.Exit()]]
            calcWindow = sg.Window('Thermochimica calculation', calcLayout, finalize=True)
            calcWindow.Element('-endtemperature-').Update(visible = False)
            calcWindow.Element('-endtemperaturelabel-').Update(visible = False)
            calcWindow.Element('-ntstep-').Update(visible = False)
            calcWindow.Element('-tsteplabel-').Update(visible = False)
            calcWindow.Element('-endpressure-').Update(visible = False)
            calcWindow.Element('-endpressurelabel-').Update(visible = False)
            calcWindow.Element('-pstep-').Update(visible = False)
            calcWindow.Element('-psteplabel-').Update(visible = False)
            calcWindow.Element('-pent-').Update(disabled = True)
            while True:
                event, values = calcWindow.read(timeout=timeout)
                eventd, valuesd = dataWindow.read(timeout=timeout)
                # print(event, values)
                if event == sg.WIN_CLOSED or event == 'Exit' or eventd == sg.WIN_CLOSED or eventd == 'Exit':
                    break
                elif event == '-tdis-':
                    calcWindow.Element('-endtemperature-').Update(visible = False)
                    calcWindow.Element('-endtemperaturelabel-').Update(visible = False)
                    calcWindow.Element('-ntstep-').Update(visible = False)
                    calcWindow.Element('-tsteplabel-').Update(visible = False)
                    calcWindow.Element('-pent-').Update(disabled = True)
                    if values['-pent-']:
                        calcWindow.Element('-pdis-').Update(value = True)
                elif event == '-ten-':
                    calcWindow.Element('-endtemperature-').Update(visible = True)
                    calcWindow.Element('-endtemperaturelabel-').Update(visible = True)
                    calcWindow.Element('-ntstep-').Update(visible = True)
                    calcWindow.Element('-tsteplabel-').Update(visible = True)
                    calcWindow.Element('-pent-').Update(disabled = False)
                elif event == '-pdis-':
                    calcWindow.Element('-endpressure-').Update(visible = False)
                    calcWindow.Element('-endpressurelabel-').Update(visible = False)
                    calcWindow.Element('-pstep-').Update(visible = False)
                    calcWindow.Element('-psteplabel-').Update(visible = False)
                elif event == '-pen-':
                    calcWindow.Element('-endpressure-').Update(visible = True)
                    calcWindow.Element('-endpressurelabel-').Update(visible = True)
                    calcWindow.Element('-pstep-').Update(visible = True)
                    calcWindow.Element('-psteplabel-').Update(visible = True)
                elif event == '-pent-':
                    calcWindow.Element('-endpressure-').Update(visible = True)
                    calcWindow.Element('-endpressurelabel-').Update(visible = True)
                    calcWindow.Element('-pstep-').Update(visible = False)
                    calcWindow.Element('-psteplabel-').Update(visible = False)
                elif event =='Run':
                    temperature = values['-temperature-']
                    pressure = values['-pressure-']
                    filename = 'inputs/pythonInput.ti'
                    masses = [0.0]*nElements
                    for i in range(nElements):
                        masses[i] = values['-'+elements[i]+'-']
                    tunit = values['-tunit-']
                    punit = values['-punit-']
                    munit = values['-munit-']
                    if values['-ten-']:
                        tend = values['-endtemperature-']
                        ntstep = values['-ntstep-']
                    else:
                        ntstep = 0
                    if values['-pen-']:
                        pend = values['-endpressure-']
                        npstep = values['-pstep-']
                    elif values['-pent-']:
                        pend = values['-endpressure-']
                        npstep = ntstep
                    else:
                        npstep = 0
                    with open(filename, 'w') as inputFile:
                        inputFile.write('! Python-generated input file for Thermochimica\n')
                        if float(ntstep) > 0:
                            tstep = (float(tend)-float(temperature))/float(ntstep)
                            inputFile.write('temperature          = ' + str(temperature) + ':' + str(tend) + ':' + str(tstep) + '\n')
                        else:
                            inputFile.write('temperature          = ' + str(temperature) + '\n')
                        if float(npstep) > 0:
                            pstep = (float(pend)-float(pressure))/float(npstep)
                            inputFile.write('pressure          = ' + str(pressure) + ':' + str(pend) + ':' + str(pstep) + '\n')
                        else:
                            inputFile.write('pressure          = ' + str(pressure) + '\n')
                        for i in range(nElements):
                            inputFile.write('mass(' + str(atomic_number_map.index(elements[i])+1) + ')           = ' + str(masses[i]) + '\n')
                        inputFile.write('temperature unit         = ' + tunit + '\n')
                        inputFile.write('pressure unit          = ' + punit + '\n')
                        inputFile.write('mass unit          = ' + munit + '\n')
                        inputFile.write('data file         = ' + datafile + '\n')
                        inputFile.write('print mode        = 2\n')
                        inputFile.write('debug mode        = .FALSE.\n')
                    thermoOut=subprocess.check_output(['./bin/ThermochimicaInputScriptMode',filename]).decode("utf-8")
                    resultOutput = [[sg.Text(thermoOut)]]
                    outWindow = sg.Window('Thermochimica output',resultOutput)
                    while True:
                        event, values = outWindow.read(timeout=timeout)
                        eventc, valuesc = calcWindow.read(timeout=timeout)
                        eventd, valuesd = dataWindow.read(timeout=timeout)
                        # print(event, values)
                        if event == sg.WIN_CLOSED or event == 'Exit' or eventd == sg.WIN_CLOSED or eventd == 'Exit'or eventc == sg.WIN_CLOSED or eventc == 'Exit':
                            break
                    outWindow.close()
            calcWindow.close()
        except:
            pass

dataWindow.close()
