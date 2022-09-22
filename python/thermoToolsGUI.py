import PySimpleGUI as sg
import os
import math
import thermoTools
import shutil

timeout = 50
inputSize = 20
buttonSize = 12
popupLocation = [300,0]

futureBlue = '#003C71'
simcoeBlue = '#0077CA'
techTangerine = '#E75D2A'
coolGrey = '#A7A8AA'
sg.theme_add_new('OntarioTech', {'BACKGROUND': futureBlue,
                                 'TEXT': 'white',
                                 'INPUT': 'white',
                                 'TEXT_INPUT': 'black',
                                 'SCROLL': coolGrey,
                                 'BUTTON': ('white', techTangerine),
                                 'PROGRESS': ('#01826B', '#D0D0D0'),
                                 'BORDER': 1,
                                 'SLIDER_DEPTH': 0,
                                 'PROGRESS_DEPTH': 0})
sg.theme('OntarioTech')

class DataWindow:
    def __init__(self,windowList,calc,parser,ext='.dat',rootDir='/data'):
        self.windowList = windowList
        self.calc = calc
        self.parser = parser
        self.ext = ext.lower()
        windowList.append(self)
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
        self.folder = os.getcwd()+rootDir
        try:
            file_list = os.listdir(self.folder)
        except:
            file_list = []
        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(self.folder, f))
            and f.lower().endswith(self.ext)
        ]
        fnames = sorted(fnames, key=str.lower)
        self.sgw = sg.Window('Thermochimica database selection', file_list_column, location = [0,0], finalize=True)
        self.sgw["-FILE LIST-"].update(fnames)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in self.windowList:
            self.windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
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
                and f.lower().endswith(self.ext)
            ]
            fnames = sorted(fnames, key=str.lower)
            self.sgw["-FILE LIST-"].update(fnames)
        elif event == "-FILE LIST-":  # A file was chosen from the listbox
            try:
                datafile = os.path.join(self.folder, values["-FILE LIST-"][0])
            except:
                return
            # Run selected parser
            self.parser(self,datafile)

class PhaseDiagramAddDataWindow:
    def __init__(self,parent,windowList):
        self.parent = parent
        self.windowList = windowList
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
        if self in self.windowList:
            self.windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
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

class PhaseDiagramMacroSettingsWindow:
    def __init__(self,parent,windowList):
        self.parent = parent
        self.windowList = windowList
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
        buttonLayout = [sg.Button('Select Macro', size = buttonSize)]
        inputNameLayout = [[sg.Text('Macro File Save Name:')],[sg.Input(key='-macroSaveName-',size=(inputSize,1)),sg.Text('.py')],[sg.Button('Set Save Name', size = buttonSize), sg.Button('Exit', size = buttonSize)]]
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
        if self in self.windowList:
            self.windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
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

def DatFileParse(parent,datafile):
    try:
        with open(datafile) as f:
            f.readline() # read comment line
            line = f.readline() # read first data line (# elements, # phases, n*# species)
            nElements = int(line[1:5])
            nSoln = int(line[6:10])
            elements = []
            while True:
                line = f.readline() # read the rest of the # species but don't need them)
                if any(c.isalpha() for c in line):
                    break
            elLen = 25 # element names are formatted 25 wide
            els = line # get the first line with letters in it
            for i in range(math.ceil(nElements/3)):
                for j in range(3):
                    elements.append(els[1+j*elLen:(1+j)*elLen].strip())
                els = f.readline() # read a line of elements (3 per line)
                # It doesn't matter now, but this reads one more line than required
    except:
        return
    # Loop over elements and check for bad names
    for el in elements:
        try:
            _ = thermoTools.atomic_number_map.index(el)+1 # get element indices in PT (i.e. # of protons)
        except ValueError:
            if len(el) > 0:
                if el[0] != 'e':
                    print(el+' not in list') # if the name is bogus (or e(phase)), discard
            elements = list(filter(lambda a: a != el, elements))
    nElements = len(elements)
    if nElements == 0:
        return
    calcWindow = parent.calc(parent,datafile,nElements,elements,True)
    parent.children.append(calcWindow)

def JSONParse(parent,datafile):
    try:
        plotWindow = parent.calc(datafile)
        parent.children.append(plotWindow)
    except:
        return