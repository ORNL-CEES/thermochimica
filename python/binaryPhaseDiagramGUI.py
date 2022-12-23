import binaryPhaseDiagramFunctions
import PySimpleGUI as sg
import os
import sys
import pickle
import copy
import matplotlib.pyplot as plt
import thermoToolsGUI
from phaseDiagramCommon import *

class CalculationWindow:
    def __init__(self, parent, datafile, nElements, elements, active):
        self.parent = parent
        self.datafile = datafile
        self.nElements = nElements
        self.elements = elements
        self.makeLayout()
        self.sgw = sg.Window(f'Phase Diagram Setup: {os.path.basename(self.datafile)}', self.layout, location = [400,0], finalize=True)
        windowList.append(self)
        self.children = []
        self.calculation = binaryPhaseDiagramFunctions.diagram(self.datafile, True, True)
        self.macro = []
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
            cancelRun = False
            grid_density = 10
            try:
                tempstep = int(values['-grid_density-'])
                if tempstep >= 0:
                    grid_density = tempstep
            except:
                pass
            if grid_density > 200:
                cancelRun = True
                confirmLayout = [[sg.Text('The selected calculation is large and may take some time.')],[sg.Button('Continue'), sg.Button('Cancel')]]
                confirmWindow = sg.Window('Large calculation confirmation', confirmLayout, location = [400,0], finalize=True, keep_on_top = True)
                while True:
                    event, values = confirmWindow.read(timeout=thermoToolsGUI.timeout)
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Continue':
                        cancelRun = False
                        break
                confirmWindow.close()
            self.calculation.makeBackup()
            pressure = 1
            try:
                tempPress = float(values['-pressure-'])
                if 1e-6 < tempPress < 1e6:
                    pressure = float(values['-pressure-'])
            except:
                pass
            tunit = values['-tunit-']
            punit = values['-punit-']
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
            el1 = values['-el1-']
            el2 = values['-el2-']
            try:
                if (str(el1) == str(el2)) or (float(tlo) == float(thi)):
                    cancelRun = True
                    repeatLayout = [[sg.Text('Values cannot be equal.')],[sg.Button('Cancel')]]
                    repeatWindow = sg.Window('Repeat value notification', repeatLayout, location = [400,0], finalize=True, keep_on_top = True)
                    while True:
                        event, values = repeatWindow.read(timeout=thermoToolsGUI.timeout)
                        if event == sg.WIN_CLOSED or event == 'Cancel':
                            break
                    repeatWindow.close()
                    return
            except ValueError:
                errorLayout = [[sg.Text('Invalid value detected.')],[sg.Button('Cancel')]]
                errorWindow = sg.Window('Invalid value notification', errorLayout, location = [400,0], finalize=True, keep_on_top = True)
                while True:
                    event, values = errorWindow.read(timeout=thermoToolsGUI.timeout)
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                errorWindow.close()
                return
            if not cancelRun:
                self.calculation.run(grid_density,grid_density,pressure,tunit,punit,0,1,tlo,thi,el1,el2,'moles',fuzzy=values["-fuzzy-"])
                self.calculation.makePlot()
                enableButtons(self)
                self.macro.append(f'macroPD.run({grid_density},{grid_density},{pressure},"{tunit}","{punit}",{0},{1},{tlo},{thi},"{el1}","{el2}","moles",fuzzy={values["-fuzzy-"]})')
        elif event =='Refine':
            refineWindow = RefineWindow(self,windowList)
            self.children.append(refineWindow)
        elif event =='Auto Refine':
            self.calculation.makeBackup()
            self.calculation.refinery()
            self.calculation.makePlot()
            self.macro.append('macroPD.makeBackup()')
            self.macro.append('macroPD.refinery()')
        elif event =='Auto Smoothen':
            self.calculation.makeBackup()
            self.sgw.Element('Undo').Update(disabled = False)
            self.calculation.autoSmooth()
            self.calculation.makePlot()
            self.macro.append('macroPD.makeBackup()')
            self.macro.append('macroPD.autoSmooth()')
        elif event =='Add Label':
            labelWindow = LabelWindow(self, windowList)
            self.children.append(labelWindow)
        elif event =='Auto Label':
            self.calculation.makeBackup()
            self.calculation.autoLabel()
            self.calculation.makePlot()
            self.sgw.Element('Remove Label').Update(disabled = False)
            self.macro.append('macroPD.makeBackup()')
            self.macro.append('macroPD.autoLabel()')
        elif event =='Remove Label':
            removeWindow = RemoveWindow(self, windowList)
            self.children.append(removeWindow)
        elif event =='Plot':
            self.calculation.makePlot()
        elif event == 'Close Plots':
            for fig in self.calculation.figureList:
                plt.close(fig=fig)
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
            settingsWindow = SettingsWindow(self, windowList)
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
            enableButtons(self)
            if len(self.calculation.labels) > 0:
                self.sgw.Element('Remove Label').Update(disabled = False)
        elif event =='Add Data':
            self.calculation.makeBackup()
            addDataWindow = thermoToolsGUI.PhaseDiagramAddDataWindow(self,windowList)
            self.children.append(addDataWindow)
        elif event =='Inspect':
            self.calculation.makeBackup()
            inspectWindow = InspectWindow(self,endMember2=self.calculation.el2,phases=self.calculation.phases,windowList=windowList)
            self.children.append(inspectWindow)
        elif event =='Export Diagram Data':
            saveDataWindow = SaveDataWindow(self)
            self.children.append(saveDataWindow)
        elif event =='Load Diagram':
            self.calculation.makeBackup()
            loadDataWindow = LoadDataWindow(self)
            self.children.append(loadDataWindow)
        elif event =='Clear Macro':
            self.macro = []
        elif event =='Export Macro':
            with open('python/' + self.macroSaveName, 'w') as f:
                f.write('import binaryPhaseDiagramFunctions\n')
                f.write('import copy\n')
                f.write(f'macroPD = binaryPhaseDiagramFunctions.diagram("{self.datafile}", False, False)\n')
                for command in self.macro:
                    f.write(f'{command}\n')
                f.write('macroPD.makePlot()\n')
        elif event =='Run Macro':
            runMacro(self)
        elif event =='Macro Settings':
            macroSettingsWindow = thermoToolsGUI.PhaseDiagramMacroSettingsWindow(self,windowList)
            self.children.append(macroSettingsWindow)
    def makeLayout(self):
        elSelectLayout = [sg.Column([[sg.Text('Element 1')],[sg.Combo(self.elements[:self.nElements],default_value=self.elements[0],key='-el1-')]],vertical_alignment='t'),
                          sg.Column([[sg.Text('Element 2')],[sg.Combo(self.elements[:self.nElements],default_value=self.elements[1],key='-el2-')]],vertical_alignment='t')]
        # xLayout        = [sg.Column([[sg.Text('Start Element 2 Concentration')],[sg.Input(key='-xlo-',size=(thermoToolsGUI.inputSize,1))],
        #                              [sg.Text('Concentration unit')],[sg.Combo(['mole fraction'],default_value='mole fraction',key='-munit-')]],vertical_alignment='t'),
        #                   sg.Column([[sg.Text('End Element 2 Concentration')],[sg.Input(key='-xhi-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t'),
        #                   sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstep-',size=(8,1))]],vertical_alignment='t')]
        tempLayout     = [sg.Column([[sg.Text('Minimum Temperature')],[sg.Input(key='-temperature-',size=(thermoToolsGUI.inputSize,1))],
                                     [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')]],vertical_alignment='t'),
                          sg.Column([[sg.Text('Maximum Temperature')],[sg.Input(key='-endtemperature-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t')]
        presLayout     = [sg.Column([[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(thermoToolsGUI.inputSize,1))],
                                     [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')]],vertical_alignment='t')]
        densityLayout  = [sg.Column([[sg.Text('Initial grid density')],[sg.Input(key='-grid_density-',size=(8,1))],
                                     [sg.Checkbox('Use fuzzy stoichiometry',key='-fuzzy-')]],vertical_alignment='t')]
        buttonLayout   = [
                            sg.Column([[sg.Button('Run', size = thermoToolsGUI.buttonSize)],
                                    [sg.Button('Undo', disabled = True, size = thermoToolsGUI.buttonSize)],
                                    [sg.Button('Close Plots', size = thermoToolsGUI.buttonSize)],
                                    [sg.Button('Add Data', size = thermoToolsGUI.buttonSize)],
                                    [sg.Button('Macro Settings', size = thermoToolsGUI.buttonSize)]],vertical_alignment='t'),
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
                         ]
        self.layout = [elSelectLayout,tempLayout,presLayout,densityLayout,buttonLayout]

class SaveData(object):
    def __init__(self,pdPoints,boundaries,phases,b,x0data,x1data,mint,maxt):
        self.pdPoints = pdPoints
        self.boundaries = boundaries
        self.phases = phases
        self.b = b
        self.x0data = x0data
        self.x1data = x1data
        self.mint = mint
        self.maxt = maxt

class SaveDataWindow:
    def __init__(self, parent):
        self.parent = parent
        windowList.append(self)
        self.children = []
        layout = [[sg.Input(key='-saveName-',size=(thermoToolsGUI.inputSize,1)), sg.Text('.pkl')],
                  [sg.Button('Save'), sg.Button('Cancel')]]
        self.sgw = sg.Window('Save Diagram Data', layout, location = [400,0], finalize=True)
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
        elif event =='Save':
            try:
                tempName = str(values['-saveName-'])
                if not tempName == '':
                    self.parent.calculation.saveDataName = tempName
            except:
                pass
            saveData = SaveData(self.parent.calculation.pdPoints,
                                self.parent.calculation.boundaries,
                                self.parent.calculation.phases,
                                self.parent.calculation.b,
                                self.parent.calculation.x0data,
                                self.parent.calculation.x1data,
                                self.parent.calculation.mint,
                                self.parent.calculation.maxt)
            with open('outputs/'+self.parent.calculation.saveDataName+'.pkl','wb') as outp:
                pickle.dump(saveData, outp, pickle.HIGHEST_PROTOCOL)
            self.close()

class LoadDataWindow:
    def __init__(self,parent):
        self.parent = parent
        windowList.append(self)
        file_list_column = [
            [
                sg.Text("Phase Diagram Data Folder"),
                sg.In(size=(25, 1), enable_events=True, key="-FOLDER-"),
                sg.FolderBrowse(),
            ],
            [
                sg.Listbox(
                    values=[], enable_events=True, size=(40, 20), key="-FILE LIST-"
                )
            ],
        ]
        self.folder = os.getcwd() + '/outputs'
        try:
            file_list = os.listdir(self.folder)
        except:
            file_list = []
        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(self.folder, f))
            and f.lower().endswith((".pkl"))
        ]
        fnames = sorted(fnames, key=str.lower)
        self.sgw = sg.Window('Phase diagram data selection', file_list_column, location = [0,0], finalize=True)
        self.sgw["-FILE LIST-"].update(fnames)
        self.children = []
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
            self.folder = values["-FOLDER-"]
            try:
                file_list = os.listdir(self.folder)
            except:
                file_list = []
            fnames = [
                f
                for f in file_list
                if os.path.isfile(os.path.join(self.folder, f))
                and f.lower().endswith((".pkl"))
            ]
            fnames = sorted(fnames, key=str.lower)
            self.sgw["-FILE LIST-"].update(fnames)
        elif event == "-FILE LIST-":  # A file was chosen from the listbox
            filename = values["-FILE LIST-"][0]
            datafile = os.path.join(self.folder, filename)
            with open(datafile, 'rb') as inp:
                self.parent.calculation.loadedDiagram = pickle.load(inp)
                self.parent.calculation.loaded = True
            self.close()

if not(os.path.isfile('bin/InputScriptMode')):
    errorLayout = [[sg.Text('No Thermochimica executable available.')],
                   [sg.Text('Either Thermochimica has not been built (run make),')],
                   [sg.Text('or this script was not executed from Thermochimica root directory.')],
                   [sg.Button('Exit')]]
    errorWindow = sg.Window('Thermochimica Error Message', errorLayout, location = [0,0], finalize=True, keep_on_top = True)
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
