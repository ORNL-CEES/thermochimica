import PySimpleGUI as sg
import subprocess

layout = [[sg.Text('Temperature')],
          [sg.Input(key='-temperature-')],
          [sg.Text('Temperature unit')],
          [sg.Combo(['K', 'C', 'F'],key='-tunit-')],
          [sg.Text('Pressure')],
          [sg.Input(key='-pressure-')],
          [sg.Text('Pressure unit')],
          [sg.Combo(['atm', 'Pa', 'bar'],key='-punit-')],
          [sg.Text('Mass unit')],
          [sg.Combo(['moles', 'kg', 'atoms', 'g'],key='-munit-')],
          [sg.Button('Run'), sg.Exit()]]

window = sg.Window('Thermochimica', layout)

while True:                             # The Event Loop
    event, values = window.read()
    print(event, values)
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    elif event =='Run':
        temperature = values['-temperature-']
        pressure = values['-pressure-']
        filename = 'inputs/pythonInput.ti'
        elements = [6,8]
        masses = [1,1.7]
        tunit = values['-tunit-']
        punit = values['-punit-']
        munit = values['-munit-']
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


window.close()
