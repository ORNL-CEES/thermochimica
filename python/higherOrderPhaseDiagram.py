import PySimpleGUI as sg
import json
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import subprocess
import scipy.optimize

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

def line_intersection(line1, lines):
    l1 = np.array(line1)
    ls = np.array(lines)
    def diff(mu):
        sum = l1[0] + mu[0]*(l1[1] - l1[0])
        for i in range(nElements - 2):
            sum -= ls[i][0] + mu[i+1]*(ls[i][1] - ls[i][0])
        return sum
    return scipy.optimize.least_squares(diff, [0.5 for i in range(nElements-1)]).x

timeout = 50
inputSize = 20

windowList = []

class DataWindow:
    def __init__(self):
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
        self.folder = os.getcwd()+'/data'
        try:
            file_list = os.listdir(self.folder)
        except:
            file_list = []
        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(self.folder, f))
            and f.lower().endswith((".dat", ".DAT"))
        ]
        self.sgw = sg.Window('Thermochimica database selection', file_list_column, location = [0,0], finalize=True)
        self.sgw["-FILE LIST-"].update(fnames)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
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
                and f.lower().endswith((".dat", ".DAT"))
            ]
            self.sgw["-FILE LIST-"].update(fnames)
        elif event == "-FILE LIST-":  # A file was chosen from the listbox
            try:
                datafile = os.path.join(
                    self.folder, values["-FILE LIST-"][0]
                )
                with open(datafile) as f:
                    f.readline() # read comment line
                    line = f.readline() # read first data line (# elements, # phases, n*# species)
                    nElements = int(line[1:5])
                    nSoln = int(line[6:10])
                    elements = []
                    for i in range(math.ceil((nSoln+3)/15)-1):
                        f.readline() # read the rest of the # species but don't need them)
                    for i in range(math.ceil(nElements/3)):
                        els = f.readline() # read a line of elements (3 per line)
                        elLen = 25 # formatted 25 wide
                        for j in range(3):
                            elements.append(els[1+j*elLen:(1+j)*elLen].strip())
                i = 0
                while i < nElements:
                    try:
                        index = atomic_number_map.index(elements[i])+1 # get element indices in PT (i.e. # of protons)
                        i = i + 1
                    except ValueError:
                        if elements[i][0] != 'e':
                            print(elements[i]+' not in list') # if the name is bogus (or e(phase)), discard
                        elements.remove(elements[i])
                        nElements = nElements - 1
                tempLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperature-',size=(inputSize,1))],
                              [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')]],
                              vertical_alignment='t'),
                              sg.Column([[sg.Text('End Temperature',key='-endtemperaturelabel-')],
                              [sg.Input(key='-endtemperature-',size=(inputSize,1))]],
                              vertical_alignment='t'),
                              sg.Column([[sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstep-',size=(8,1))]],
                              vertical_alignment='t')]
                presLayout = [sg.Column([[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(inputSize,1))],
                              [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')]],
                              vertical_alignment='t')]
                elem1Layout = [[sg.Text('Composition 1')]]
                elem2Layout = [[sg.Text('Composition 2')]]
                for i in range(nElements):
                    elem1Layout.append([sg.Text(elements[i])])
                    elem1Layout.append([sg.Input(key='-'+elements[i]+'1-',size=(inputSize,1))])
                for i in range(nElements):
                    elem2Layout.append([sg.Text(elements[i])])
                    elem2Layout.append([sg.Input(key='-'+elements[i]+'2-',size=(inputSize,1))])
                if (nElements < 8):
                    calcLayout = [tempLayout,
                                  presLayout,
                                  [sg.Column(elem1Layout),sg.Column(elem2Layout)],
                                  [sg.Text('Mass unit')],
                                  [sg.Combo(['moles'],default_value='moles',key='-munit-')],
                                  [sg.Checkbox('Save JSON',key='-json-')],
                                  [sg.Button('Run'), sg.Exit()]]
                else:
                    calcLayout = [tempLayout,
                                  presLayout,
                                  [sg.Column(elem1Layout,vertical_alignment='t', scrollable = True, vertical_scroll_only = True, expand_y = True)],
                                  [sg.Text('Mass unit')],
                                  [sg.Combo(['moles'],default_value='moles',key='-munit-')],
                                  [sg.Checkbox('Save JSON',key='-json-')],
                                  [sg.Button('Run'), sg.Exit()]]
                calcWindow = CalculationWindow(calcLayout,datafile,nElements,elements)
                self.children.append(calcWindow)
            except:
                pass

class CalculationWindow:
    def __init__(self, windowLayout, datafile, nElements, elements):
        windowList.append(self)
        self.datafile = datafile
        self.nElements = nElements
        self.elements = elements
        self.sgw = sg.Window(f'Thermochimica calculation: {os.path.basename(datafile)}', windowLayout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event =='Run':
                temperature = values['-temperature-']
                pressure = values['-pressure-']
                filename = 'inputs/pythonInput.ti'
                masses1 = [0.0]*self.nElements
                masses2 = [0.0]*self.nElements
                elementsUsed = []
                for i in range(self.nElements):
                    if values['-'+self.elements[i]+'1-'] != '':
                        masses1[i] = values['-'+self.elements[i]+'1-']
                    if values['-'+self.elements[i]+'2-'] != '':
                        masses2[i] = values['-'+self.elements[i]+'2-']
                    if masses1[i] or masses2[i]:
                        elementsUsed.append(self.elements[i])
                tunit = values['-tunit-']
                punit = values['-punit-']
                munit = values['-munit-']
                tend = values['-endtemperature-']
                ntstep = values['-ntstep-']
                with open(filename, 'w') as inputFile:
                    inputFile.write('! Python-generated input file for Thermochimica\n')
                    inputFile.write(f'data file         = {datafile}\n')
                    inputFile.write(f'temperature unit         = {tunit}\n')
                    inputFile.write(f'pressure unit          = {punit}\n')
                    inputFile.write(f'mass unit          = \'{munit}\'\n')
                    inputFile.write(f'nEl         = {nElements} \n')
                    inputFile.write(f'iEl         = {atomic_number_map.index(el1)+1} {atomic_number_map.index(el2)+1} {atomic_number_map.index(el3)+1}\n')
                    inputFile.write(f'nCalc       = {len(xs)*len(temperature)}\n')
                    for t in temperature:
                        for x in xs:
                            inputFile.write(f'{t} {pressure} {(1-x3)*(1-x)} {(1-x3)*x} {x3}\n')
                subprocess.run(['./bin/RunCalculationList',filename])

dataWindow = DataWindow()
while len(windowList) > 0:
    for window in windowList:
        window.read()
