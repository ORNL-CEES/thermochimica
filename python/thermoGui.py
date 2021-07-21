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

dataWindow = sg.Window('Thermochimica database selection', file_list_column)

calcIter = 0
while True:
    event, values = dataWindow.read()
    print(event, values)
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
                values["-FOLDER-"], values["-FILE LIST-"][0]
            )
            with open(datafile) as f:
                f.readline()
                line = f.readline()
                nElements = int(line[1:5])
                nSoln = int(line[6:10])
                elements = ['']*nElements
                elIt = 0
                for i in range(math.ceil((nSoln+3)/15)-1):
                    f.readline()
                for i in range(math.ceil(nElements/3)):
                    els = f.readline()
                    elements[elIt] = els[1:25].strip()
                    elIt = elIt + 1
                    if elIt == nElements:
                        break
                    elements[elIt] = els[26:50].strip()
                    elIt = elIt + 1
                    if elIt == nElements:
                        break
                    elements[elIt] = els[51:75].strip()
                    elIt = elIt + 1
                    if elIt == nElements:
                        break

            i = 0
            while i < nElements:
                try:
                    index = atomic_number_map.index(elements[i])+1
                    print(elements[i],index)
                    i = i + 1
                except ValueError:
                    print(elements[i]+' not in list')
                    elements.remove(elements[i])
                    print(elements)
                    nElements = nElements - 1
            tempLayout = [[sg.Text('Temperature')],[sg.Input(key='-temperature-',size=(16,1))],
                          [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')]]
            presLayout = [[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(16,1))],
                          [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')]]
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
            calcWindow = sg.Window('Thermochimica calculation', calcLayout)
            while True:
                event, values = calcWindow.read()
                print(event, values)
                if event == sg.WIN_CLOSED or event == 'Exit':
                    break
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
                    with open(filename, 'w') as inputFile:
                        inputFile.write('! Python-generated input file for Thermochimica\n')
                        inputFile.write('temperature          = ' + str(temperature) + '\n')
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
                        event, values = outWindow.read()
                        print(event, values)
                        if event == sg.WIN_CLOSED or event == 'Exit':
                            break
                    outWindow.close()
            calcWindow.close()
        except:
            pass

dataWindow.close()
