import pseudoBinaryPhaseDiagramFunctions
import PySimpleGUI as sg
import matplotlib.pyplot as plt
import os
import sys
import shutil
import copy
import thermoToolsGUI

# For boundaries of phase regions where both sides have (# phases) < (# elements), only plot points within phaseFractionTol of the boundary
phaseFractionTol = 1e-2
# Below this tolerance, set phase fraction = 0
phaseIncludeTol = 1e-8

class CalculationWindow:
    def __init__(self, parent, datafile, nElements, elements, active):
        self.parent = parent
        self.datafile = datafile
        self.nElements = nElements
        self.elements = elements
        self.active = active
        if self.active:
            self.makeLayout()
            self.sgw = sg.Window(f'Phase Diagram Setup: {os.path.basename(self.datafile)}', self.layout, location = [400,0], finalize=True)
            windowList.append(self)
        self.children = []
        self.calculation = pseudoBinaryPhaseDiagramFunctions.diagram(self.datafile, True, True)
        self.macro = []
        self.showExperiment = True
        self.backup = []
        self.macroSaveName = 'macroPhaseDiagram.py'
    def close(self):
        for child in self.children:
            child.close()
        for fig in self.calculation.figureList:
            plt.close(fig=fig)
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event =='Run':
                self.calculation.makeBackup()
                tlo = 300
                try:
                    templo = float(values['-temperature-'])
                    if 295 <= templo <= 6000:
                        tlo = templo
                except:
                    pass
                thi = 1000
                try:
                    temphi = float(values['-endtemperature-'])
                    if 295 <= temphi <= 6000:
                        thi = temphi
                except:
                        pass
                ntstep = 10
                try:
                    tempstep = int(values['-ntstep-'])
                    if tempstep >= 0:
                        ntstep = tempstep
                except:
                    pass
                nxstep = 10
                try:
                    tempstep = int(values['-nxstep-'])
                    if tempstep >= 0:
                        nxstep = tempstep
                except:
                    pass
                pressure = 1
                try:
                    tempPress = float(values['-pressure-'])
                    if 1e-6 < tempPress < 1e6:
                        pressure = float(values['-pressure-'])
                except:
                    pass
                if values['-pressure-'] != '':
                    pressure = values['-pressure-']
                masses1 = [0.0]*self.nElements
                masses2 = [0.0]*self.nElements
                elementsUsed = []
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
                    if (masses1[i] > 0.0) or (masses2[i] > 0.0):
                        elementsUsed.append(self.elements[i])
                for i in reversed(range(self.nElements)):
                    if not self.elements[i] in elementsUsed:
                        del masses1[i]
                        del masses2[i]
                nElementsUsed = len(elementsUsed)
                massLabels = ['','']
                for i in range(nElementsUsed):
                    if masses1[i] > 0:
                        massLabels[0] += elementsUsed[i]
                        if masses1[i] != 1:
                            if int(masses1[i]) == masses1[i]:
                                massLabels[0] += f'$_{ {int(masses1[i])} }$'
                            else:
                                massLabels[0] += f'$_{ {masses1[i]} }$'
                    if masses2[i] > 0:
                        massLabels[1] += elementsUsed[i]
                        if masses2[i] != 1:
                            if int(masses2[i]) == masses2[i]:
                                massLabels[1] += f'$_{ {int(masses2[i])} }$'
                            else:
                                massLabels[1] += f'$_{ {masses2[i]} }$'
                sum1 = sum(masses1)
                sum2 = sum(masses2)
                if (sum1 == 0) or (sum2 == 0) or (min(masses1) < 0) or (min(masses2) < 0):
                    alertLayout = [[sg.Text('One or more input masses are invalid.')],[sg.Button('Continue')]]
                    alertWindow = sg.Window('Invalid Mass Alert', alertLayout, location = [400,0], finalize=True)
                    while True:
                        event, values = alertWindow.read(timeout=thermoToolsGUI.timeout)
                        if event == sg.WIN_CLOSED or event == 'Continue':
                            break
                    alertWindow.close()
                    return
                masses1 = [mass / sum1 for mass in masses1]
                masses2 = [mass / sum2 for mass in masses2]
                plane = [masses1,masses2]
                tunit = values['-tunit-']
                # Check temperature unit for shift
                if tunit == 'K':
                    tshift = 0
                elif tunit == 'C':
                    tshift = 273.15
                # Initially reverse shift (opposite at plot)
                mint = tlo + tshift
                maxt = thi + tshift
                punit = values['-punit-']
                munit = values['-munit-']
                self.calculation.initRun(pressure,tunit,punit,plane,sum1,sum2,mint,maxt,elementsUsed,massLabels,munit,tshift)
                self.macro.append(f'macroPD.initRun({pressure},\'{tunit}\',\'{punit}\',{plane},{sum1},{sum2},{mint},{maxt},{elementsUsed},{massLabels},\'{munit}\',{tshift})')
                self.calculation.runCalc(0,1,nxstep,tlo,thi,ntstep)
                self.macro.append(f'macroPD.runCalc(0,1,{nxstep},{tlo},{thi},{ntstep})')
                self.calculation.processPhaseDiagramData()
                self.macro.append(f'macroPD.processPhaseDiagramData()')
                self.calculation.makePlot()
                self.sgw.Element('Refine').Update(disabled = False)
                self.sgw.Element('Add Label').Update(disabled = False)
                self.sgw.Element('Plot').Update(disabled = False)
                self.sgw.Element('Export Plot').Update(disabled = False)
                self.sgw.Element('Undo').Update(disabled = False)
                self.sgw.Element('Add Data').Update(disabled = False)
        elif event =='Refine':
            refineWindow = RefineWindow(self)
            self.children.append(refineWindow)
        elif event =='Add Label':
            labelWindow = LabelWindow(self)
            self.children.append(labelWindow)
        elif event =='Remove Label':
            removeWindow = RemoveWindow(self)
            self.children.append(removeWindow)
        elif event =='Plot':
            self.calculation.makePlot()
            self.macro.append(f'macroPD.makePlot()')
        elif event =='Export Plot':
            exportStatus = self.calculation.exportPlot()
            if exportStatus:
                errorLayout = [[sg.Text('The export failed, try changing plot settings.')],[sg.Button('Continue'), sg.Button('Cancel')]]
                errorWindow = sg.Window('Plot export failed', errorLayout, location = [400,0], finalize=True, keep_on_top = True)
                while True:
                    event, values = errorWindow.read(timeout=thermoToolsGUI.timeout)
                    if event == sg.WIN_CLOSED or event == 'Continue':
                        break
                errorWindow.close()
            else:
                self.macro.append(f'macroPD.makePlot()')
                self.macro.append(f'macroPD.exportPlot()')
        elif event =='Plot Settings':
            settingsWindow = SettingsWindow(self)
            self.children.append(settingsWindow)
        elif event =='Undo':
            for fig in self.calculation.figureList:
                plt.close(fig=fig)
            backup = copy.deepcopy(self.calculation.backup)
            self.calculation = self.calculation.backup
            self.calculation.backup = backup
            self.macro.append('backup = copy.deepcopy(macroPD.backup)')
            self.macro.append('macroPD = macroPD.backup')
            self.macro.append('macroPD.backup = backup')
            self.calculation.makePlot()
            self.sgw.Element('Refine').Update(disabled = False)
            self.sgw.Element('Add Label').Update(disabled = False)
            self.sgw.Element('Plot').Update(disabled = False)
            self.sgw.Element('Export Plot').Update(disabled = False)
            if len(self.calculation.labels) > 0:
                self.sgw.Element('Remove Label').Update(disabled = False)
        elif event =='Clear Macro':
            self.macro = []
        elif event =='Export Macro':
            with open('python/' + self.macroSaveName, 'w') as f:
                f.write('import pseudoBinaryPhaseDiagramFunctions\n')
                f.write('import copy\n')
                f.write('import numpy as np\n')
                f.write(f'macroPD = pseudoBinaryPhaseDiagramFunctions.diagram("{self.datafile}", False, False)\n')
                for command in self.macro:
                    f.write(f'{command}\n')
                f.write('macroPD.makePlot()\n')
        elif event =='Run Macro':
            if 'macroPhaseDiagram' in sys.modules:
                del sys.modules['macroPhaseDiagram']
            import macroPhaseDiagram
            self.calculation = macroPhaseDiagram.macroPD
            self.calculation.active = True
            self.calculation.interactivePlot = True
            self.sgw.Element('Refine').Update(disabled = False)
            self.sgw.Element('Add Label').Update(disabled = False)
            self.sgw.Element('Plot').Update(disabled = False)
            self.sgw.Element('Export Plot').Update(disabled = False)
            self.sgw.Element('Undo').Update(disabled = False)
            self.sgw.Element('Add Data').Update(disabled = False)
        elif event =='Add Data':
            self.calculation.makeBackup()
            addDataWindow = AddDataWindow(self)
            self.children.append(addDataWindow)
        elif event =='Macro Settings':
            macroSettingsWindow = MacroSettingsWindow(self)
            self.children.append(macroSettingsWindow)
    def makeLayout(self):
        tempLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperature-',size=(thermoToolsGUI.inputSize,1))],
                      [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')]],
                      vertical_alignment='t'),
                      sg.Column([[sg.Text('End Temperature',key='-endtemperaturelabel-')],
                      [sg.Input(key='-endtemperature-',size=(thermoToolsGUI.inputSize,1))]],
                      vertical_alignment='t'),
                      sg.Column([[sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstep-',size=(8,1))]],
                      vertical_alignment='t')]
        presLayout = [sg.Column([[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(thermoToolsGUI.inputSize,1))],
                      [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')]],
                      vertical_alignment='t')]
        elem1Layout = [[sg.Text('Composition 1')]]
        elem2Layout = [[sg.Text('Composition 2')]]
        for el in self.elements:
            elem1Layout.append([sg.Text(el)])
            elem1Layout.append([sg.Input(key=f'-{el}1-',size=(thermoToolsGUI.inputSize,1))])
        for el in self.elements:
            elem2Layout.append([sg.Text(el)])
            elem2Layout.append([sg.Input(key=f'-{el}2-',size=(thermoToolsGUI.inputSize,1))])
        if (self.nElements < 8):
            elemLayout = [sg.Column(elem1Layout),sg.Column(elem2Layout)]
        else:
            elemLayout = [sg.Column(elem1Layout,vertical_alignment='t', scrollable = True, vertical_scroll_only = True, expand_y = True),
                          sg.Column(elem2Layout,vertical_alignment='t', scrollable = True, vertical_scroll_only = True, expand_y = True)]
        self.layout = [tempLayout,
                       presLayout,
                       elemLayout,
                       [sg.Text('# of steps')],[sg.Input(key='-nxstep-',size=(8,1))],
                       [sg.Text('Mass unit')],
                       [sg.Combo(['moles'],default_value='moles',key='-munit-')],
                       [
            sg.Column([[sg.Button('Run', size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Undo', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Exit(size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Add Data', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Macro Settings', size = thermoToolsGUI.buttonSize)]
                      ],vertical_alignment='t'),
            sg.Column([[sg.Button('Refine', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Auto Refine', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Auto Smoothen', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Inspect', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Run Macro', size = thermoToolsGUI.buttonSize)]],vertical_alignment='t'),
            sg.Column([[sg.Button('Add Label', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Auto Label', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Remove Label', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Load Diagram', size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Export Macro', size = thermoToolsGUI.buttonSize)]],vertical_alignment='t'),
            sg.Column([[sg.Button('Plot', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Export Plot', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Plot Settings', size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Export Diagram Data', disabled = True, size = thermoToolsGUI.buttonSize)],
                       [sg.Button('Clear Macro', size = thermoToolsGUI.buttonSize)]],vertical_alignment='t')
            ]]

class RefineWindow:
    def __init__(self, parent):
        self.parent = parent
        xRefLayout    = [sg.Column([[sg.Text('Start Concentration')],[sg.Input(key='-xlor-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t'),
                          sg.Column([[sg.Text('End Concentration')],[sg.Input(key='-xhir-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t'),
                          sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstepr-',size=(8,1))]],vertical_alignment='t')]
        tempRefLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperaturer-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t'),
                        sg.Column([[sg.Text('End Temperature')],[sg.Input(key='-endtemperaturer-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t'),
                        sg.Column([[sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstepr-',size=(8,1))]],vertical_alignment='t')]
        refineLayout = [xRefLayout,tempRefLayout,[sg.Button('Refine'), sg.Button('Cancel')]]
        self.sgw = sg.Window('Phase diagram refinement', refineLayout, location = [400,0], finalize=True)
        windowList.append(self)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            self.close()
        elif event =='Refine':
            cancelRun = False
            ntstep = 10
            try:
                tempstep = int(values['-ntstepr-'])
                if tempstep >= 0:
                    ntstep = tempstep
            except:
                pass
            nxstep = 10
            try:
                tempstep = int(values['-nxstepr-'])
                if tempstep >= 0:
                    nxstep = tempstep
            except:
                pass
            if (float(ntstep) * float(nxstep)) > 50000:
                cancelRun = True
                confirmLayout = [[sg.Text('The selected calculation is large and may take some time.')],[sg.Button('Continue'), sg.Button('Cancel')]]
                confirmWindow = sg.Window('Large calculation confirmation', confirmLayout, location = [400,0], finalize=True)
                while True:
                    event, values = confirmWindow.read(timeout=thermoToolsGUI.timeout)
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Continue':
                        cancelRun = False
                        break
                confirmWindow.close()
            xlo = 0
            try:
                templo = float(values['-xlor-'])
                if 0 <= templo <= 1:
                    xlo = templo
            except:
                pass
            xhi = 1
            try:
                temphi = float(values['-xhir-'])
                if 0 <= temphi <= 1:
                    xhi = temphi
            except:
                pass
            tlo = 300
            try:
                templo = float(values['-temperaturer-'])
                if 295 <= templo <= 6000:
                    tlo = templo
            except:
                pass
            thi = 1000
            try:
                temphi = float(values['-endtemperaturer-'])
                if 295 <= temphi <= 6000:
                    thi = temphi
            except:
                    pass

            # refine x-coords are going to come in scaled to axis
            xlo = xlo
            xhi = xhi
            if not self.parent.calculation.normalizeX:
                if xlo > 0:
                    xlo = 1/(1+((1-xlo)/xlo)*(self.parent.calculation.sum1/self.parent.calculation.sum2))
                if xhi > 0:
                    xhi = 1/(1+((1-xhi)/xhi)*(self.parent.calculation.sum1/self.parent.calculation.sum2))
            if not cancelRun:
                self.parent.calculation.makeBackup()
                self.parent.macro.append(f'macroPD.makeBackup()')
                self.parent.sgw.Element('Undo').Update(disabled = False)
                self.parent.calculation.runCalc(xlo,xhi,nxstep,tlo,thi,ntstep)
                self.parent.macro.append(f'macroPD.runCalc({xlo},{xhi},{nxstep},{tlo},{thi},{ntstep})')
                self.parent.calculation.processPhaseDiagramData()
                self.parent.macro.append(f'macroPD.processPhaseDiagramData()')
                self.parent.calculation.makePlot()

class LabelWindow:
    def __init__(self, parent):
        self.parent = parent
        xLabLayout  = [[sg.Text(f'{self.parent.calculation.massLabels[1].translate({ord(i):None for i in "{}_$"})} Concentration')],[sg.Input(key='-xlab-',size=(thermoToolsGUI.inputSize,1))]]
        tLabLayout  = [[sg.Text('Temperature')],[sg.Input(key='-tlab-',size=(thermoToolsGUI.inputSize,1))]]
        labelLayout = [xLabLayout,tLabLayout,[sg.Button('Add Label'), sg.Button('Cancel')]]
        self.sgw = sg.Window('Add phase label', labelLayout, location = [400,0], finalize=True)
        windowList.append(self)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            self.close()
        elif event =='Add Label':
            try:
                try:
                    xlab = float(values['-xlab-'])
                except ValueError:
                    num, den = values['-xlab-'].split('/')
                    xlab = float(num)/float(den)
                tlab = float(values['-tlab-'])
                if (0 <= xlab <= 1) and (295 <= tlab <= 6000):
                    self.parent.calculation.makeBackup()
                    self.parent.macro.append(f'macroPD.makeBackup()')
                    self.parent.sgw.Element('Undo').Update(disabled = False)
                    self.parent.calculation.addLabel(xlab,tlab)
                    self.parent.macro.append(f'macroPD.addLabel({xlab},{tlab})')
                    self.parent.calculation.processPhaseDiagramData()
                    self.parent.macro.append(f'macroPD.processPhaseDiagramData()')
                    self.parent.calculation.makePlot()
                    self.parent.sgw.Element('Remove Label').Update(disabled = False)
            except:
                pass

class RemoveWindow:
    def __init__(self, parent):
        self.parent = parent
        headingsLayout = [[sg.Text('Label Text',   size = [55,1],justification='left'),
                               sg.Text('Concentration',size = [15,1],justification='center'),
                               sg.Text('Temperature',  size = [15,1],justification='center'),
                               sg.Text('Remove Label?',size = [15,1])]]
        labelListLayout = []
        for i in range(len(self.parent.calculation.labels)):
            labelListLayout.append([[sg.Text(self.parent.calculation.labels[i][1],size = [55,1],justification='left'),
                                        sg.Text("{:.3f}".format(self.parent.calculation.labels[i][0][0]),size = [15,1],justification='center'),
                                        sg.Text("{:.0f}".format(self.parent.calculation.labels[i][0][1]),size = [15,1],justification='center'),
                                        sg.Checkbox('',key='-removeLabel'+str(i)+'-',pad=[[40,0],[0,0]])]])
        removeLayout = [headingsLayout,labelListLayout,[sg.Button('Remove Label(s)'), sg.Button('Cancel')]]
        self.sgw = sg.Window('Add phase label', removeLayout, location = [400,0], finalize=True)
        windowList.append(self)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            self.close()
        if event == 'Remove Label(s)':
            self.parent.calculation.makeBackup()
            self.parent.macro.append(f'macroPD.makeBackup()')
            self.parent.sgw.Element('Undo').Update(disabled = False)
            tempLength = len(self.parent.calculation.labels)
            for i in reversed(range(tempLength)):
                if values['-removeLabel'+str(i)+'-']:
                    del self.parent.calculation.labels[i]
                    self.parent.macro.append(f'del macroPD.labels[{i}]')
            if len(self.parent.calculation.labels) == 0:
                self.parent.sgw.Element('Remove Label').Update(disabled = True)
            self.parent.calculation.makePlot()
            self.close()

class SettingsWindow:
    def __init__(self, parent):
        self.parent = parent
        windowList.append(self)
        if self.parent.calculation.plotMarker == '-':
            line  = True
            point = False
            both  = False
        elif self.parent.calculation.plotMarker == '.':
            line  = False
            point = True
            both  = False
        else:
            line  = False
            point = False
            both  = True
        if self.parent.calculation.plotColor == 'colorful':
            colorful = True
            bland    = False
        else:
            colorful = False
            bland    = True
        if self.parent.calculation.experimentColor == 'colorful':
            expcolorful = True
            expbland    = False
        else:
            expcolorful = False
            expbland    = True
        settingsLayout = [[sg.Text('Marker Style:')],
                          [sg.Radio('Lines', 'mstyle', default=line,  enable_events=True, key='-mline-')],
                          [sg.Radio('Points','mstyle', default=point, enable_events=True, key='-mpoint-')],
                          [sg.Radio('Both',  'mstyle', default=both,  enable_events=True, key='-mboth-')],
                          [sg.Text('Plot Colors:')],
                          [sg.Radio('Colorful', 'mcolor', default=colorful, enable_events=True, key='-mcolorful-')],
                          [sg.Radio('Black',    'mcolor', default=bland,    enable_events=True, key='-mbland-')],
                          [sg.Text('Experimental Data Colors:')],
                          [sg.Radio('Colorful', 'mexpcolor', default=expcolorful, enable_events=True, key='-mexpcolorful-')],
                          [sg.Radio('Black',    'mexpcolor', default=expbland,    enable_events=True, key='-mexpbland-')],
                          [sg.Text('Show:')],
                          [sg.Checkbox('Experimental Data', default=self.parent.calculation.showExperiment, key='-showExperiment-'),
                           sg.Checkbox('Data Legend', default=self.parent.calculation.showExperimentLegend, key='-showExperimentLegend-')
                        #   ,sg.Checkbox('Loaded Diagram', default=self.parent.calculation.showLoaded, key='-showLoaded-')
                          ],
                        #   [sg.Text('Auto-Label Settings:')],
                        #   [sg.Checkbox('1-Phase Regions', default=self.parent.calculation.label1phase, key='-label1phase-'),
                        #    sg.Checkbox('2-Phase Regions', default=self.parent.calculation.label2phase, key='-label2phase-')],
                          [sg.Text('Export Filename'),sg.Input(key='-filename-',size=(thermoToolsGUI.inputSize,1))],
                          [sg.Text('Export Format'),sg.Combo(['png', 'pdf', 'ps', 'eps', 'svg'],default_value='png',key='-format-')],
                          [sg.Text('Export DPI'),sg.Input(key='-dpi-',size=(thermoToolsGUI.inputSize,1))],
                          [sg.Button('Accept')]]
        self.sgw = sg.Window('Plot Settings', settingsLayout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        # Log settings in macro before closing
        self.parent.macro.append(f'macroPD.plotMarker = "{self.parent.calculation.plotMarker}"')
        self.parent.macro.append(f'macroPD.plotColor = "{self.parent.calculation.plotColor}"')
        self.parent.macro.append(f'macroPD.experimentColor = "{self.parent.calculation.experimentColor}"')
        self.parent.macro.append(f'macroPD.showExperiment = {self.parent.calculation.showExperiment}')
        self.parent.macro.append(f'macroPD.exportFileName = "{self.parent.calculation.exportFileName}"')
        self.parent.macro.append(f'macroPD.exportFormat = "{self.parent.calculation.exportFormat}"')
        self.parent.macro.append(f'macroPD.exportDPI = {self.parent.calculation.exportDPI}')
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED:
            self.close()
        elif event == '-mline-':
            self.parent.calculation.plotMarker = '-'
        elif event =='-mpoint-':
            self.parent.calculation.plotMarker = '.'
        elif event =='-mboth-':
            self.parent.calculation.plotMarker = '.-'
        elif event =='-mcolorful-':
            self.parent.calculation.plotColor = 'colorful'
        elif event =='-mbland-':
            self.parent.calculation.plotColor = 'bland'
        elif event =='-mexpcolorful-':
            self.parent.calculation.experimentColor = 'colorful'
        elif event =='-mexpbland-':
            self.parent.calculation.experimentColor = 'bland'
        elif event =='Accept':
            self.parent.calculation.showExperiment = values['-showExperiment-']
            self.parent.calculation.showExperimentLegend = values['-showExperimentLegend-']
            # self.parent.calculation.showLoaded = values['-showLoaded-']
            # self.parent.calculation.label1phase = values['-label1phase-']
            # self.parent.calculation.label2phase = values['-label2phase-']
            try:
                if str(values['-filename-']) != '':
                    self.parent.calculation.exportFileName = str(values['-filename-'])
            except:
                pass
            self.parent.calculation.exportFormat = values['-format-']
            try:
                tempDPI = int(values['-dpi-'])
                if tempDPI > 0 > 10000:
                    self.parent.calculation.exportDPI = int(values['-dpi-'])
            except:
                pass
            self.parent.calculation.makePlot()
            self.close()

class AddDataWindow:
    def __init__(self,parent):
        self.parent = parent
        windowList.append(self)
        file_list_column = [
            [
                sg.Text("Experimental Data Folder"),
                sg.In(size=(25, 1), enable_events=True, key="-FOLDER-"),
                sg.FolderBrowse(),
            ],
            [
                sg.Listbox(
                    values=[], enable_events=True, size=(40, 20), key="-FILE LIST-"
                )
            ],
        ]
        self.folder = os.getcwd()
        try:
            file_list = os.listdir(self.folder)
        except:
            file_list = []
        buttonLayout = [sg.Button('Add Data'), sg.Button('Exit')]
        addDataLayout = [file_list_column,buttonLayout]
        self.sgw = sg.Window('Experimental data selection', addDataLayout, location = [0,0], finalize=True)
        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(self.folder, f))
            and f.lower().endswith((".csv"))
        ]
        fnames = sorted(fnames, key=str.lower)
        self.sgw["-FILE LIST-"].update(fnames)
        self.children = []
        self.filename = ''
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
        elif event == "-FOLDER-":
            self.filename = ''
            self.folder = values["-FOLDER-"]
            try:
                file_list = os.listdir(self.folder)
            except:
                file_list = []

            fnames = [
                f
                for f in file_list
                if os.path.isfile(os.path.join(self.folder, f))
                and f.lower().endswith((".csv"))
            ]
            fnames = sorted(fnames, key=str.lower)
            self.sgw["-FILE LIST-"].update(fnames)
        elif event == "-FILE LIST-":  # A file was chosen from the listbox
            self.filename = values["-FILE LIST-"][0]
        elif event == 'Add Data':
            if not self.filename:
                return
            datafile = os.path.join(self.folder, self.filename)
            expName = self.filename.split('.',1)[0]
            self.parent.calculation.addData(datafile,expName)
            self.parent.macro.append(f'macroPD.addData("{datafile}","{expName}")')

class MacroSettingsWindow:
    def __init__(self,parent):
        self.parent = parent
        windowList.append(self)
        file_list_column = [
            [
                sg.Text("Macro Folder"),
                sg.In(size=(25, 1), enable_events=True, key="-FOLDER-"),
                sg.FolderBrowse(),
            ],
            [
                sg.Listbox(
                    values=[], enable_events=True, size=(40, 20), key="-FILE LIST-"
                )
            ],
        ]
        self.folder = os.getcwd() + '/python'
        try:
            file_list = os.listdir(self.folder)
        except:
            file_list = []
        buttonLayout = [sg.Button('Select Macro', size = thermoToolsGUI.buttonSize)]
        inputNameLayout = [[sg.Text('Macro File Save Name:')],[sg.Input(key='-macroSaveName-',size=(thermoToolsGUI.inputSize,1)),sg.Text('.py')],[sg.Button('Set Save Name', size = thermoToolsGUI.buttonSize), sg.Button('Exit', size = thermoToolsGUI.buttonSize)]]
        addDataLayout = [file_list_column,buttonLayout,inputNameLayout]
        self.sgw = sg.Window('Macro file', addDataLayout, location = [0,0], finalize=True)
        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(self.folder, f))
            and f.lower().endswith((".py"))
        ]
        fnames = sorted(fnames, key=str.lower)
        self.sgw["-FILE LIST-"].update(fnames)
        self.children = []
        self.filename = ''
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
        elif event == "-FOLDER-":
            self.filename = ''
            self.folder = values["-FOLDER-"]
            try:
                file_list = os.listdir(self.folder)
            except:
                file_list = []

            fnames = [
                f
                for f in file_list
                if os.path.isfile(os.path.join(self.folder, f))
                and f.lower().endswith((".py"))
            ]
            fnames = sorted(fnames, key=str.lower)
            self.sgw["-FILE LIST-"].update(fnames)
        elif event == "-FILE LIST-":  # A file was chosen from the listbox
            self.filename = values["-FILE LIST-"][0]
        elif event == 'Select Macro':
            if not self.filename:
                return
            datafile = os.path.join(self.folder, self.filename)
            # I don't want to mess with import logic, so rename/overwrite file in one spot instead
            shutil.copy(datafile,os.getcwd() + '/python/' + 'macroPhaseDiagram.py')
            self.close()
        elif event == 'Set Save Name':
            if not values["-macroSaveName-"]:
                return
            self.parent.macroSaveName = values["-macroSaveName-"] + '.py'

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
dataWindow = thermoToolsGUI.DataWindow(windowList,CalculationWindow)
while len(windowList) > 0:
    for window in windowList:
        window.read()
