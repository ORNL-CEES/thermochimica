import PySimpleGUI as sg
import subprocess
import os
import sys
import numpy as np
import shutil
import thermoTools
import thermoToolsGUI

class CalculationWindow:
    def __init__(self, parent, datafile, nElements, elements, active):
        windowList.append(self)
        self.datafile = datafile
        self.nElements = nElements
        self.elements = elements
        self.makeLayout()
        self.sgw = sg.Window(f'Thermochimica calculation: {os.path.basename(self.datafile)}', self.layout, location = [400,0], finalize=True)
        self.children = []
        self.exportFileName = 'thermoout'
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event == '-tdis-':
            self.sgw.Element('-endtemperature-').Update(visible = False)
            self.sgw.Element('-endtemperaturelabel-').Update(visible = False)
            self.sgw.Element('-ntstep-').Update(visible = False)
            self.sgw.Element('-tsteplabel-').Update(visible = False)
            self.sgw.Element('-pent-').Update(disabled = True)
            if values['-pent-']:
                self.sgw.Element('-pdis-').Update(value = True)
                self.sgw.Element('-endpressure-').Update(visible = False)
                self.sgw.Element('-endpressurelabel-').Update(visible = False)
        elif event == '-ten-':
            self.sgw.Element('-endtemperature-').Update(visible = True)
            self.sgw.Element('-endtemperaturelabel-').Update(visible = True)
            self.sgw.Element('-ntstep-').Update(visible = True)
            self.sgw.Element('-tsteplabel-').Update(visible = True)
            self.sgw.Element('-pent-').Update(disabled = False)
        elif event == '-pdis-':
            self.sgw.Element('-endpressure-').Update(visible = False)
            self.sgw.Element('-endpressurelabel-').Update(visible = False)
            self.sgw.Element('-pstep-').Update(visible = False)
            self.sgw.Element('-psteplabel-').Update(visible = False)
        elif event == '-pen-':
            self.sgw.Element('-endpressure-').Update(visible = True)
            self.sgw.Element('-endpressurelabel-').Update(visible = True)
            self.sgw.Element('-pstep-').Update(visible = True)
            self.sgw.Element('-psteplabel-').Update(visible = True)
        elif event == '-pent-':
            self.sgw.Element('-endpressure-').Update(visible = True)
            self.sgw.Element('-endpressurelabel-').Update(visible = True)
            self.sgw.Element('-pstep-').Update(visible = False)
            self.sgw.Element('-psteplabel-').Update(visible = False)
        elif event == '-mdis-':
            self.sgw.Element('-composition2-').Update(visible = False)
            self.sgw.Element('-xstepcol-').Update(visible = False)
            self.sgw.Element('-ten-').Update(disabled = False)
            self.sgw.Element('-pen-').Update(disabled = False)
        elif event == '-men-':
            self.sgw.Element('-composition2-').Update(visible = True)
            self.sgw.Element('-xstepcol-').Update(visible = True)
            self.sgw.Element('-ten-').Update(disabled = True)
            self.sgw.Element('-endtemperature-').Update(visible = False)
            self.sgw.Element('-endtemperaturelabel-').Update(visible = False)
            self.sgw.Element('-ntstep-').Update(visible = False)
            self.sgw.Element('-tsteplabel-').Update(visible = False)
            self.sgw.Element('-pen-').Update(disabled = True)
            self.sgw.Element('-pent-').Update(disabled = True)
            self.sgw.Element('-endpressure-').Update(visible = False)
            self.sgw.Element('-endpressurelabel-').Update(visible = False)
            self.sgw.Element('-pstep-').Update(visible = False)
            self.sgw.Element('-psteplabel-').Update(visible = False)
            self.sgw.Element('-tdis-').Update(value = True)
            self.sgw.Element('-pdis-').Update(value = True)
        elif event == 'Run':
                temperature = 300
                try:
                    templo = float(values['-temperature-'])
                    if 295 <= templo <= 6000:
                        temperature = templo
                except:
                    pass
                pressure = 1
                try:
                    tempPress = float(values['-pressure-'])
                    if 1e-6 < tempPress < 1e6:
                        pressure = float(values['-pressure-'])
                except:
                    pass
                filename = 'inputs/pythonInput.ti'
                masses1 = [0.0]*self.nElements
                masses2 = [0.0]*self.nElements
                for i in range(self.nElements):
                    try:
                        tempMass = float(values[f'-{self.elements[i]}1-'])
                        if 0 < tempMass:
                            masses1[i] = tempMass
                    except:
                        pass
                    try:
                        tempMass = float(values[f'-{self.elements[i]}2-'])
                        if 0 < tempMass:
                            masses2[i] = tempMass
                    except:
                        pass
                tunit = values['-tunit-']
                punit = values['-punit-']
                munit = values['-munit-']
                if values['-ten-']:
                    tend = values['-endtemperature-']
                    ntstep = values['-ntstep-']
                    try:
                        ntstep = float(ntstep)
                    except ValueError:
                        ntstep = 1
                else:
                    tend = None
                    ntstep = 0
                if values['-pen-']:
                    pend = values['-endpressure-']
                    npstep = values['-pstep-']
                    try:
                        npstep = float(npstep)
                    except ValueError:
                        npstep = 1
                elif values['-pent-']:
                    pend = values['-endpressure-']
                    npstep = ntstep
                else:
                    pend = None
                    npstep = 0
                if values['-men-']:
                    nxstep = 10
                    try:
                        tempstep = int(values['-nxstep-'])
                        if tempstep >= 0:
                            nxstep = tempstep
                    except:
                        pass

                    xs = np.array([np.linspace(masses1[i],masses2[i],nxstep) for i in range(self.nElements)]).T
                    calcList = []
                    for x in xs:
                        calc = [temperature,pressure]
                        calc.extend(x)
                        calcList.append(calc)

                    thermoTools.WriteRunCalculationList(filename,self.datafile,self.elements,calcList,tunit=tunit,punit=punit,munit=munit,heatCapacity=values["-cp_h_s-"],writeJson=values["-json-"])
                    thermoOut = subprocess.check_output(['./bin/RunCalculationList',filename]).decode("utf-8")
                else:
                    thermoTools.WriteInputScript(filename,self.datafile,self.elements,temperature,tend,ntstep,pressure,pend,npstep,masses1,tunit=tunit,punit=punit,munit=munit,heatCapacity=values["-cp_h_s-"],writeJson=values["-json-"],stepTogether=values["-pent-"])
                    thermoOut = subprocess.check_output(['./bin/InputScriptMode',filename]).decode("utf-8")
                nLines = thermoOut.count('\n')
                if (nLines < 5000):
                    resultOutput = [[sg.Column([[sg.Multiline(thermoOut, size = (65, nLines),font='TkFixedFont')]], size = (None, 800), scrollable = True, vertical_scroll_only = True)]]
                else:
                    resultOutput = [[sg.Text('Output is too large to display')]]
                resultWindow = ResultWindow(resultOutput)
                self.children.append(resultWindow)
                if values['-json-']:
                    try:
                        shutil.copy2('thermoout.json', f'{self.exportFileName}.json')
                    except:
                        pass
        elif event == 'Set name':
            setNameLayout = [[sg.Input(key='-jsonname-',size=(thermoToolsGUI.inputSize,1)),sg.Text('.json')],[sg.Button('Accept'), sg.Button('Cancel')]]
            setNameWindow = sg.Window('Set JSON name', setNameLayout, location = [400,0], finalize=True)
            while True:
                event, values = setNameWindow.read(timeout=thermoToolsGUI.timeout)
                if event == sg.WIN_CLOSED or event == 'Cancel':
                    break
                elif event == 'Accept':
                    try:
                        if str(values['-jsonname-']) != '':
                            self.exportFileName = str(values['-jsonname-'])
                    except:
                        pass
                    break
            setNameWindow.close()
    def makeLayout(self):
        tempLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperature-',size=(thermoToolsGUI.inputSize,1))],
                      [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')]],vertical_alignment='t'),
                      sg.Column([[sg.Text('End Temperature',key='-endtemperaturelabel-',visible=False)],
                      [sg.Input(key='-endtemperature-',size=(thermoToolsGUI.inputSize,1),visible=False)],
                      [sg.Text('Temperature range:')],
                      [sg.Radio('Disabled', 'trange', default=True, enable_events=True, key='-tdis-')],
                      [sg.Radio('Enabled', 'trange', default=False, enable_events=True, key='-ten-')]],vertical_alignment='t'),
                      sg.Column([[sg.Text('# of steps',key='-tsteplabel-',visible=False)],[sg.Input(key='-ntstep-',size=(8,1),visible=False)]],vertical_alignment='t')]
        presLayout = [sg.Column([[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(thermoToolsGUI.inputSize,1))],
                      [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')]],vertical_alignment='t'),
                      sg.Column([[sg.Text('End Pressure',key='-endpressurelabel-',visible=False)],[sg.Input(key='-endpressure-',size=(thermoToolsGUI.inputSize,1),visible=False)],
                      [sg.Text('Pressure range:')],
                      [sg.Radio('Disabled', 'prange', default=True, enable_events=True, key='-pdis-')],
                      [sg.Radio('Enabled', 'prange', default=False, enable_events=True, key='-pen-')],
                      [sg.Radio('Enabled, step\nwith temperature', 'prange', default=False, disabled=True, enable_events=True, key='-pent-')]],vertical_alignment='t'),
                      sg.Column([[sg.Text('# of steps',key='-psteplabel-',visible=False)],[sg.Input(key='-pstep-',size=(8,1),visible=False)]],vertical_alignment='t')
                      ]
        elem1Layout = [[sg.Text('Composition 1')]]
        elem2Layout = [[sg.Text('Composition 2')]]
        for el in self.elements:
            elem1Layout.append([sg.Text(el)])
            elem1Layout.append([sg.Input(key=f'-{el}1-',size=(thermoToolsGUI.inputSize,1))])
        for el in self.elements:
            elem2Layout.append([sg.Text(el)])
            elem2Layout.append([sg.Input(key=f'-{el}2-',size=(thermoToolsGUI.inputSize,1))])
        massLayout = [sg.Column([
                        [sg.Text('Mass unit')],
                        [sg.Combo(['moles', 'kg', 'atoms', 'g'],default_value='moles',key='-munit-')]
                        ],vertical_alignment='t'),
                      sg.Column([
                        [sg.Text('Composition range:')],
                        [sg.Radio('Disabled', 'mrange', default=True, enable_events=True, key='-mdis-')],
                        [sg.Radio('Enabled', 'mrange', default=False, enable_events=True, key='-men-')]
                        ],vertical_alignment='t'),
                        sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstep-',size=(8,1))]],key='-xstepcol-',visible=False,vertical_alignment='t')]
        if (self.nElements < 8):
            self.layout = [tempLayout,
                          presLayout,
                          [sg.Column(elem1Layout,vertical_alignment='t'),
                           sg.Column(elem2Layout,key='-composition2-',visible=False,vertical_alignment='t')],
                          massLayout,
                          [sg.Checkbox('Save JSON',key='-json-'), sg.Button('Set name')],
                          [sg.Checkbox('Calculate heat capacity, entropy, and enthalpy',key='-cp_h_s-')],
                          [sg.Button('Run'), sg.Exit()]]
        else:
            self.layout = [tempLayout,
                          presLayout,
                          [sg.Column(elem1Layout,vertical_alignment='t', scrollable = True, vertical_scroll_only = True, expand_y = True),
                           sg.Column(elem2Layout,vertical_alignment='t', scrollable = True, vertical_scroll_only = True, expand_y = True,key='-composition2-',visible=False)],
                          massLayout,
                          [sg.Checkbox('Save JSON',key='-json-'), sg.Button('Set name')],
                          [sg.Checkbox('Calculate heat capacity, entropy, and enthalpy',key='-cp_h_s-')],
                          [sg.Button('Run'), sg.Exit()]]

class ResultWindow:
    def __init__(self, layout):
        windowList.append(self)
        self.sgw = sg.Window('Thermochimica output',layout, location = [825,0], finalize=True, resizable=True)
    def close(self):
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()

if not(os.path.isfile('bin/InputScriptMode')):
    errorLayout = [[sg.Text('No Thermochimica executable available.')],
                   [sg.Text('Either Thermochimica has not been built (run make),')],
                   [sg.Text('or this script was not executed from Thermochimica root directory.')],
                   [sg.Button('Exit')]]
    errorWindow = sg.Window('Thermochimica Error Message', errorLayout, location = [0,0], finalize=True)
    while True:
        event, values = errorWindow.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            break
    errorWindow.close()
    sys.exit()

windowList = []
dataWindow = thermoToolsGUI.DataWindow(windowList,CalculationWindow,thermoToolsGUI.DatFileParse)
while len(windowList) > 0:
    for window in windowList:
        window.read()
