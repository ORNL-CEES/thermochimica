import PySimpleGUI as sg
import json
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import math
import os
import sys
import subprocess
import copy
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.geometry import MultiPoint
from shapely.geometry import LineString
from shapely.prepared import prep
from shapely.ops import split
from functools import reduce
import operator

timeout = 50
inputSize = 20
buttonSize = 12

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

phaseIncludeTol = 1e-8

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

def fmt(x,pos=None):
    return '{:.2f}'.format(1-x)

class DataWindow:
    def __init__(self):
        windowList.append(self)
        file_list_column = [
            [
                sg.Text('Database Folder'),
                sg.In(size=(25, 1), enable_events=True, key='-FOLDER-'),
                sg.FolderBrowse(),
            ],
            [
                sg.Listbox(values=[], enable_events=True, size=(40, 20), key='-FILE LIST-')
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
        fnames = sorted(fnames, key=str.lower)
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
            fnames = sorted(fnames, key=str.lower)
            self.sgw["-FILE LIST-"].update(fnames)
        elif event == "-FILE LIST-":  # A file was chosen from the listbox
            try:
                datafile = os.path.join(self.folder, values["-FILE LIST-"][0])
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
            i = 0
            for el in elements:
                try:
                    index = atomic_number_map.index(el)+1 # get element indices in PT (i.e. # of protons)
                except ValueError:
                    if len(el) > 0:
                        if el[0] != 'e':
                            print(el+' not in list') # if the name is bogus (or e(phase)), discard
                    elements = list(filter(lambda a: a != el, elements))
            nElements = len(elements)
            if nElements == 0:
                return
            calcWindow = CalculationWindow(self,datafile,nElements,elements,True)
            self.children.append(calcWindow)

class CalculationWindow:
    def __init__(self, parent, datafile, nElements, elements, active):
        self.parent = parent
        self.datafile = datafile
        self.nElements = nElements
        self.elements = elements
        self.x1 = np.empty([0])
        self.x2 = np.empty([0])
        self.p1 = []
        self.p2 = []
        self.points3 = []
        self.points1 = []
        self.el1 = ''
        self.el2 = ''
        self.el3 = ''
        self.tunit = 'K'
        self.punit = 'atm'
        self.munit = 'moles'
        self.active = active
        if self.active:
            self.makeLayout()
            self.sgw = sg.Window(f'Phase Diagram Setup: {os.path.basename(self.datafile)}', self.layout, location = [400,0], finalize=True)
            windowList.append(self)
        self.children = []
        self.labels = []
        outline = MultiPolygon([])
        self.temperature = 300
        self.pressure = 1
        self.inputFileName = 'inputs/pythonPhaseDiagramInput.ti'
        self.outputFileName = 'thermoout.json'
        self.plotMarker = '-'
        self.plotColor = 'colorful'
        self.backup = []
        self.currentPlot = []
        self.exportFormat = 'png'
        self.exportFileName = 'thermochimicaPhaseDiagram'
        self.exportDPI = 300
        self.resRef = 4
        self.resSmooth = 4
        self.figureList = []
        self.tielines = True
        self.tiegap = 1/27
    def close(self):
        for child in self.children:
            child.close()
        for fig in self.figureList:
            plt.close(fig=fig)
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event =='Run':
            cancelRun = False
            nxstep = 10
            try:
                tempstep = int(values['-nxstep-'])
                if tempstep >= 0:
                    nxstep = tempstep
            except:
                pass
            if (float(nxstep)) > 500:
                cancelRun = True
                confirmLayout = [[sg.Text('The selected calculation is large and may take some time.')],[sg.Button('Continue'), sg.Button('Cancel')]]
                confirmWindow = sg.Window('Large calculation confirmation', confirmLayout, location = [400,0], finalize=True, keep_on_top = True)
                while True:
                    event, values = confirmWindow.read(timeout=timeout)
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Continue':
                        cancelRun = False
                        break
                confirmWindow.close()
            self.makeBackup()
            self.temperature = 300
            try:
                templo = float(values['-temperature-'])
                if 295 <= templo <= 6000:
                    self.temperature = templo
            except:
                pass
            self.pressure = 1
            try:
                tempPress = float(values['-pressure-'])
                if 1e-6 < tempPress < 1e6:
                    self.pressure = float(values['-pressure-'])
            except:
                pass
            self.tunit = values['-tunit-']
            self.punit = values['-punit-']
            self.el1 = values['-el1-']
            self.el2 = values['-el2-']
            self.el3 = values['-el3-']
            try:
                if (str(self.el1) == str(self.el2)) or (str(self.el2) == str(self.el3)) or (str(self.el1) == str(self.el3)):
                    cancelRun = True
                    repeatLayout = [[sg.Text('Elements cannot be identical.')],[sg.Button('Cancel')]]
                    repeatWindow = sg.Window('Repeat value notification', repeatLayout, location = [400,0], finalize=True, keep_on_top = True)
                    while True:
                        event, values = repeatWindow.read(timeout=timeout)
                        if event == sg.WIN_CLOSED or event == 'Cancel':
                            break
                        elif event == 'Continue':
                            cancelRun = False
                            break
                    repeatWindow.close()
            except ValueError:
                errorLayout = [[sg.Text('Invalid value detected.')],[sg.Button('Cancel')]]
                errorWindow = sg.Window('Invalid value notification', errorLayout, location = [400,0], finalize=True, keep_on_top = True)
                while True:
                    event, values = errorWindow.read(timeout=timeout)
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                errorWindow.close()
                return
            if not cancelRun:
                self.writeInputFile(0,1,0,1,nxstep)
                self.x1 = np.empty([0])
                self.x2 = np.empty([0])
                self.p1 = []
                self.p2 = []
                self.points3 = []
                self.points1 = []
                self.labels = []
                self.resRef = 7
                self.resSmooth = 7
                self.runCalc()
                self.makePlot()
                self.sgw.Element('Refine').Update(disabled = False)
                self.sgw.Element('Auto Refine').Update(disabled = False)
                self.sgw.Element('Auto Smoothen').Update(disabled = False)
                self.sgw.Element('Add Label').Update(disabled = False)
                self.sgw.Element('Auto Label').Update(disabled = False)
                self.sgw.Element('Plot').Update(disabled = False)
                self.sgw.Element('Undo').Update(disabled = False)
        elif event =='Refine':
            xRefLayout    = [sg.Column([[sg.Text(f'Start {self.el1} Concentration')],[sg.Input(key='-xlor1-',size=(inputSize,1))]],vertical_alignment='t'),
                          sg.Column([[sg.Text(f'End {self.el1} Concentration')],[sg.Input(key='-xhir1-',size=(inputSize,1))]],vertical_alignment='t'),
                          sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstepr-',size=(8,1))]],vertical_alignment='t')]
            tempRefLayout = [sg.Column([[sg.Text(f'Start {self.el2} Concentration')],[sg.Input(key='-xlor2-',size=(inputSize,1))]],vertical_alignment='t'),
                          sg.Column([[sg.Text(f'End {self.el2} Concentration')],[sg.Input(key='-xhir2-',size=(inputSize,1))]],vertical_alignment='t')]
            refineLayout = [xRefLayout,tempRefLayout,[sg.Button('Refine'), sg.Button('Cancel')]]
            refineWindow = RefineWindow(self, refineLayout)
            self.children.append(refineWindow)
        elif event =='Auto Refine':
            self.makeBackup()
            self.sgw.Element('Undo').Update(disabled = False)
            self.autoRefine(self.resRef**2)
            # self.writeInputFile(0,1,0,1,self.resRef**2)
            # self.runCalc()
            self.resRef += 1
            self.makePlot()
        elif event =='Auto Smoothen':
            self.makeBackup()
            self.sgw.Element('Undo').Update(disabled = False)
            self.autoRefine2Phase(self.resSmooth**2)
            self.makePlot()
            self.resSmooth += 1
        elif event =='Add Label':
            xLabLayout    = [[sg.Text(f'{self.el1} Concentration')],[sg.Input(key='-x1lab-',size=(inputSize,1))]]
            tLabLayout = [[sg.Text(f'{self.el2} Concentration')],[sg.Input(key='-x2lab-',size=(inputSize,1))]]
            labelLayout = [xLabLayout,tLabLayout,[sg.Button('Add Label'), sg.Button('Cancel')]]
            labelWindow = LabelWindow(self,labelLayout)
            self.children.append(labelWindow)
        elif event =='Auto Label':
            self.makeBackup()
            self.autoLabel()
            self.makePlot()
            self.sgw.Element('Remove Label').Update(disabled = False)
        elif event =='Remove Label':
            headingsLayout = [[sg.Text('Label Text',   size = [55,1],justification='left'),
                               sg.Text(f'{self.el1} Concentration',size = [15,1],justification='center'),
                               sg.Text(f'{self.el2} Concentration',  size = [15,1],justification='center'),
                               sg.Text('Remove Label?',size = [15,1])]]
            labelListLayout = []
            for i in range(len(self.labels)):
                labelListLayout.append([[sg.Text(self.labels[i][1],size = [55,1],justification='left'),
                                         sg.Text("{:.3f}".format(float(self.labels[i][0][0])),size = [15,1],justification='center'),
                                         sg.Text("{:.3f}".format(float(self.labels[i][0][1])),size = [15,1],justification='center'),
                                         sg.Checkbox('',key='-removeLabel'+str(i)+'-',pad=[[40,0],[0,0]])]])
            removeLayout = [headingsLayout,labelListLayout,[sg.Button('Remove Label(s)'), sg.Button('Cancel')]]
            removeWindow = RemoveWindow(self, removeLayout)
            self.children.append(removeWindow)
        elif event =='Plot':
            self.makePlot()
        elif event =='Export Plot':
            self.exportPlot()
        elif event =='Plot Settings':
            if self.plotMarker == '-':
                line  = True
                point = False
                both  = False
            elif self.plotMarker == '.':
                line  = False
                point = True
                both  = False
            else:
                line  = False
                point = False
                both  = True
            if self.plotColor == 'colorful':
                colorful = True
                bland    = False
            else:
                colorful = False
                bland    = True
            settingsLayout = [[sg.Text('Marker Style:')],
                             [sg.Radio('Lines', 'mstyle', default=line,  enable_events=True, key='-mline-')],
                             [sg.Radio('Points','mstyle', default=point, enable_events=True, key='-mpoint-')],
                             [sg.Radio('Both',  'mstyle', default=both,  enable_events=True, key='-mboth-')],
                             [sg.Text('Plot Colors:')],
                             [sg.Radio('Colorful', 'mcolor', default=colorful, enable_events=True, key='-mcolorful-')],
                             [sg.Radio('Black',    'mcolor', default=bland,    enable_events=True, key='-mbland-')],
                             [sg.Checkbox('Tielines', default=self.tielines, key='-tielines-'),
                              sg.Text('Density:'),sg.Input(key='-tiedensity-',size=(inputSize,1))],
                             [sg.Text('Export Filename'),sg.Input(key='-filename-',size=(inputSize,1))],
                             [sg.Text('Export Format'),sg.Combo(['png', 'pdf', 'ps', 'eps', 'svg'],default_value='png',key='-format-')],
                             [sg.Text('Export DPI'),sg.Input(key='-dpi-',size=(inputSize,1))],
                             [sg.Button('Accept')]]
            settingsWindow = SettingsWindow(self, settingsLayout)
            self.children.append(settingsWindow)
        elif event =='Undo':
            self.backup.activate()
            self.close()
    def processPhaseDiagramData(self):
        f = open(self.outputFileName,)
        data = json.load(f)
        f.close()
        if list(data.keys())[0] != '1':
            print('Output does not contain data series')
            exit()
        x1 = self.x1.tolist()
        x2 = self.x2.tolist()
        for i in list(data.keys()):
            try:
                nPhases = 0
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                            nPhases += 1
            except:
                continue
            # 1-phase data points (edges only)
            if nPhases == 1:
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                            boundPhase = phaseName
                            tempComps = [0,0,0]
                            if self.el1 in list(data[i][phaseType][phaseName]['elements'].keys()):
                                tempComps[0] = data[i][phaseType][phaseName]['elements'][self.el1]['mole fraction of phase by element']
                            if self.el2 in list(data[i][phaseType][phaseName]['elements'].keys()):
                                tempComps[1] = data[i][phaseType][phaseName]['elements'][self.el2]['mole fraction of phase by element']
                            if self.el3 in list(data[i][phaseType][phaseName]['elements'].keys()):
                                tempComps[2] = data[i][phaseType][phaseName]['elements'][self.el3]['mole fraction of phase by element']
                # only record points on a diagram boundary
                if min(tempComps) > 0:
                    continue
                self.points1.append([[tempComps[0],tempComps[1]],boundPhase])
            # 2-phase data points
            if nPhases == 2:
                boundPhases = []
                boundComps = []
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                            boundPhases.append(phaseName)
                            tempComps = [0,0]
                            if self.el1 in list(data[i][phaseType][phaseName]['elements'].keys()):
                                tempComps[0] = data[i][phaseType][phaseName]['elements'][self.el1]['mole fraction of phase by element']
                            if self.el2 in list(data[i][phaseType][phaseName]['elements'].keys()):
                                tempComps[1] = data[i][phaseType][phaseName]['elements'][self.el2]['mole fraction of phase by element']
                            boundComps.append(tempComps)
                x1.append(boundComps[0])
                x2.append(boundComps[1])
                self.p1.append(boundPhases[0])
                self.p2.append(boundPhases[1])
            # 3-phase data points
            if nPhases == 3:
                boundPhases = []
                boundComps = []
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                            boundPhases.append(phaseName)
                            tempComps = [0,0]
                            if self.el1 in list(data[i][phaseType][phaseName]['elements'].keys()):
                                tempComps[0] = data[i][phaseType][phaseName]['elements'][self.el1]['mole fraction of phase by element']
                            if self.el2 in list(data[i][phaseType][phaseName]['elements'].keys()):
                                tempComps[1] = data[i][phaseType][phaseName]['elements'][self.el2]['mole fraction of phase by element']
                            boundComps.append(tempComps)
                # Record triplet (check values to avoid duplicating)
                if len(self.points3) > 0:
                    mindist = np.sqrt(min([np.linalg.norm(np.array(boundComps)-np.array(p[0])) for p in self.points3]))
                else:
                    mindist = 100
                if mindist > 1e-2:
                    self.points3.append([boundComps,boundPhases])
                    # Add first pair
                    x1.append(boundComps[0])
                    x2.append(boundComps[1])
                    self.p1.append(boundPhases[0])
                    self.p2.append(boundPhases[1])
                    # Add second pair
                    x1.append(boundComps[0])
                    x2.append(boundComps[2])
                    self.p1.append(boundPhases[0])
                    self.p2.append(boundPhases[2])
                    # Add third pair
                    x1.append(boundComps[1])
                    x2.append(boundComps[2])
                    self.p1.append(boundPhases[1])
                    self.p2.append(boundPhases[2])
        self.x1 = np.array(x1)
        self.x2 = np.array(x2)
    def runCalc(self):
        print('Thermochimica calculation initiated.')
        subprocess.run(['./bin/Phase3DiagramDataGen',self.inputFileName])
        print('Thermochimica calculation finished.')
        self.processPhaseDiagramData()
    def makePlot(self):
        boundaries = []
        phases = []
        b = []
        for i in range(len(self.p1)):
            # If a miscibility gap label has been used unnecessarily, remove it
            if self.p1[i].find('#') > 0:
                if not(self.p1[i][0:self.p1[i].find('#')] == self.p2[i]):
                    self.p1[i] = self.p1[i][0:self.p1[i].find('#')]
            if self.p2[i].find('#') > 0:
                if not(self.p2[i][0:self.p2[i].find('#')] == self.p1[i]):
                    self.p2[i] = self.p2[i][0:self.p2[i].find('#')]
            repeat = False
            for j in range(len(boundaries)):
                if (boundaries[j][0] == self.p1[i]) and (boundaries[j][1] == self.p2[i]):
                    b.append(j)
                    repeat = True
            if not(repeat):
                boundaries.append([self.p1[i],self.p2[i]])
                b.append(len(boundaries)-1)

        for i in range(len(boundaries)):
            repeat1 = False
            repeat2 = False
            for j in range(len(phases)):
                if (boundaries[i][0] == phases[j]):
                    repeat1 = True
                if (boundaries[i][1] == phases[j]):
                    repeat2 = True
            if not(repeat1 or boundaries[i][0].find('#') > 0):
                phases.append(boundaries[i][0])
            if not(repeat2 or boundaries[i][1].find('#') > 0):
                phases.append(boundaries[i][1])

        # Start figure
        fig = plt.figure()
        plt.ion()
        ax = fig.add_axes([0.125, 0.1, 0.75, 0.85])

        # plot 2-phase region boundaries
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(boundaries))))
        for j in range(len(boundaries)):
            if self.plotColor == 'colorful':
                c = next(color)
            else:
                c = 'k'
            inds = [i for i, k in enumerate(b) if k == j]
            if len(inds) < 2:
                continue
            v1 = np.array([self.x1[inds[1],0],self.x1[inds[1],1],0])
            v2 = np.array([self.x2[inds[1],0],self.x2[inds[1],1],0])
            normal = np.cross(v1-v2,np.array([0,0,1]))
            order = [inds[i] for i, k in sorted(enumerate(self.x1[inds].tolist()), key=lambda coord: coord[1][0]*normal[0]+coord[1][1]*normal[1])]
            x1s = self.x1[order]
            x2s = self.x2[order]
            # Draw tie lines before adding points and flipping order
            if self.tielines:
                lastline1 = x1s[0] + (x1s[-1]-x1s[0]) * self.tiegap / 2
                lastline2 = x2s[0] + (x2s[-1]-x2s[0]) * self.tiegap / 2
                endline1 = x1s[-1] - (x1s[-1]-x1s[0]) * self.tiegap / 2
                endline2 = x2s[-1] - (x2s[-1]-x2s[0]) * self.tiegap / 2
                for i in range(len(x1s)-1):
                    gap = min(max(np.linalg.norm(x1s[i]-lastline1),np.linalg.norm(x2s[i]-lastline2)),max(np.linalg.norm(x1s[i]-endline1),np.linalg.norm(x2s[i]-endline2)))
                    if gap > self.tiegap:
                        ax.plot([1-(x1s[i,0]+x1s[i,1]/2),1-(x2s[i,0]+x2s[i,1]/2)],[x1s[i,1],x2s[i,1]],'--k')
                        lastline1 = x1s[i]
                        lastline2 = x2s[i]
            # Reverse second half and add end points from opposite side to form box
            x2s = np.flip(x2s,axis=0)
            x1s = np.append(x1s,[x2s[0]],axis=0)
            x2s = np.append(x2s,[x1s[0]],axis=0)
            ax.plot(1-(x1s[:,0]+x1s[:,1]/2),x1s[:,1],self.plotMarker,c=c)
            ax.plot(1-(x2s[:,0]+x2s[:,1]/2),x2s[:,1],self.plotMarker,c=c)

        ax.plot([0,0.5,1,0],[0,1,0,0],'k-')
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_title(f'{self.el1} - {self.el2} - {self.el3} ternary phase diagram')
        ax.set_xlabel(f'Mole fraction {self.el1}')
        ax.set_ylabel(f'Mole fraction {self.el2}')
        for lab in self.labels:
            plt.text(1-(lab[0][0]+lab[0][1]/2),lab[0][1],lab[1], ha='center')
        # reverse x tick labels for ternary plot
        ax.xaxis.set_major_formatter(FuncFormatter(fmt))
        ax2 = ax.twinx()
        ax2.yaxis.set_major_formatter(FuncFormatter(fmt))
        ax2.set_ylabel(f'Mole fraction {self.el3}')
        plt.sca(ax)
        plt.show()
        plt.pause(0.001)
        self.currentPlot = fig
        self.figureList.append(fig)
        self.sgw.Element('Export Plot').Update(disabled = False)
    def writeInputFile(self,xlo1,xhi1,xlo2,xhi2,nxstep):
        if xlo1 > xhi1:
            temp = xlo1
            xlo1 = xhi1
            xhi1 = temp
        if xlo2 > xhi2:
            temp = xlo2
            xlo2 = xhi2
            xhi2 = temp
        with open(self.inputFileName, 'w') as inputFile:
            inputFile.write('! Python-generated input file for Thermochimica\n')
            if float(nxstep) > 0:
                xstep1 = (float(xhi1)-float(xlo1))/float(nxstep)
            else:
                xstep1 = 0
            if float(nxstep) > 0:
                xstep2 = (float(xhi2)-float(xlo2))/float(nxstep)
            else:
                xstep2 = 0
            inputFile.write(f'x1               = {str(xlo1)}:{str(xhi1)}:{str(xstep1)}\n')
            inputFile.write(f'x2               = {str(xlo2)}:{str(xhi2)}:{str(xstep2)}\n')
            inputFile.write(f'temperature      = {str(self.temperature)}\n')
            inputFile.write(f'pressure         = {str(self.pressure)}\n')
            inputFile.write(f'temperature unit = {self.tunit}\n')
            inputFile.write(f'pressure unit    = {self.punit}\n')
            inputFile.write(f'mass unit        = \'{self.munit}\'\n')
            inputFile.write(f'iEl              = {str(atomic_number_map.index(self.el1)+1)} {str(atomic_number_map.index(self.el2)+1)} {str(atomic_number_map.index(self.el3)+1)}\n')
            inputFile.write(f'data file        = {self.datafile}\n')
    def addLabel(self,x1lab,x2lab):
        with open(self.inputFileName, 'w') as inputFile:
            inputFile.write('! Python-generated input file for Thermochimica\n')
            inputFile.write(f'data file         = {self.datafile}\n')
            inputFile.write(f'temperature unit  = {self.tunit}\n')
            inputFile.write(f'pressure unit     = {self.punit}\n')
            inputFile.write(f'mass unit         = \'{self.munit}\'\n')
            inputFile.write( 'nEl               = 3 \n')
            inputFile.write(f'iEl               = {atomic_number_map.index(self.el1)+1} {atomic_number_map.index(self.el2)+1} {atomic_number_map.index(self.el3)+1}\n')
            inputFile.write(f'nCalc             = 1\n')
            inputFile.write(f'{self.temperature} {self.pressure} {x1lab} {x2lab} {1-x1lab-x2lab}\n')
        print('Thermochimica calculation initiated.')
        subprocess.run(['./bin/RunCalculationList',self.inputFileName])
        print('Thermochimica calculation finished.')
        f = open(self.outputFileName,)
        data = json.load(f)
        f.close()
        if list(data.keys())[0] != '1':
            print('Output does not contain data series')
            exit()
        labelName = []
        for phaseName in list(data['1']['solution phases'].keys()):
            if (data['1']['solution phases'][phaseName]['moles'] > 0):
                labelName.append(phaseName)
        for phaseName in list(data['1']['pure condensed phases'].keys()):
            if (data['1']['pure condensed phases'][phaseName]['moles'] > 0):
                labelName.append(phaseName)
        self.labels.append([[x1lab,x2lab],'+'.join(labelName)])
        self.processPhaseDiagramData()
    def autoRefine(self,res):
        outline = Polygon([[0,0],[0,1],[1,0],[0,0]])
        maxArea = 0
        boundaries = []
        b = []
        for i in range(len(self.p1)):
            # If a miscibility gap label has been used unnecessarily, remove it
            if self.p1[i].find('#') > 0:
                if not(self.p1[i][0:self.p1[i].find('#')] == self.p2[i]):
                    self.p1[i] = self.p1[i][0:self.p1[i].find('#')]
            if self.p2[i].find('#') > 0:
                if not(self.p2[i][0:self.p2[i].find('#')] == self.p1[i]):
                    self.p2[i] = self.p2[i][0:self.p2[i].find('#')]
            repeat = False
            for j in range(len(boundaries)):
                if (boundaries[j][0] == self.p1[i]) and (boundaries[j][1] == self.p2[i]):
                    b.append(j)
                    repeat = True
            if not(repeat):
                boundaries.append([self.p1[i],self.p2[i]])
                b.append(len(boundaries)-1)

        # Make list of phases
        phases = []
        for i in range(len(boundaries)):
            repeat1 = False
            repeat2 = False
            for j in range(len(phases)):
                if (boundaries[i][0] == phases[j]):
                    repeat1 = True
                if (boundaries[i][1] == phases[j]):
                    repeat2 = True
            if not(repeat1 or boundaries[i][0].find('#') > 0):
                phases.append(boundaries[i][0])
            if not(repeat2 or boundaries[i][1].find('#') > 0):
                phases.append(boundaries[i][1])

        # find and subtract 1-phase regions
        for phase in phases:
            inds1 = [i for i, k in enumerate(self.p1) if k == phase]
            inds2 = [i for i, k in enumerate(self.p2) if k == phase]
            points = []
            if len(inds1) > 0:
                points.extend(self.x1[inds1].tolist())
            if len(inds2) > 0:
                points.extend(self.x2[inds2].tolist())
            for point in self.points1:
                if point[1] == phase:
                    points.append(point[0])
            if len(point) < 3:
                continue
            average = np.average(np.array(points),axis=0)
            sortpoints = sorted(points, key=lambda coord: (-135 - math.degrees(math.atan2(*tuple(map(operator.sub, coord, average))[::-1]))) % 360)
            phaseOutline = Polygon(sortpoints).buffer(0)
            try:
                outline = outline - phaseOutline
            except:
                continue

        # find and subtract 2-phase regions
        for j in range(len(boundaries)):
            inds = [i for i, k in enumerate(b) if k == j]
            if len(inds) < 2:
                continue
            v1 = np.array([self.x1[inds[1],0],self.x1[inds[1],1],0])
            v2 = np.array([self.x2[inds[1],0],self.x2[inds[1],1],0])
            normal = np.cross(v1-v2,np.array([0,0,1]))
            order = [inds[i] for i, k in sorted(enumerate(self.x1[inds].tolist()), key=lambda coord: coord[1][0]*normal[0]+coord[1][1]*normal[1])]
            x1s = self.x1[order]
            x2s = self.x2[order]
            # Reverse second half and add end points from opposite side to form box
            x2s = np.flip(x2s,axis=0)
            x1s = np.append(x1s,x2s,axis=0)
            x1s = np.append(x1s,[x1s[0]],axis=0)
            phaseOutline = Polygon(x1s).buffer(0)
            try:
                outline = outline - phaseOutline
            except:
                continue

        # find and subtract 3-phase regions
        for point in self.points3:
            phaseOutline = Polygon(point[0])
            try:
                outline = outline - phaseOutline
            except:
                continue

        xs = []
        ys = []
        subres = int(np.ceil(np.sqrt(res)))
        try:
            oxlo, otlo, oxhi, othi = outline.bounds
        except:
            return
        xindices = np.linspace(oxlo, oxhi, subres)
        yindices = np.linspace(otlo, othi, subres)
        horizontal_splitters = [LineString([(x, yindices[0]), (x, yindices[-1])]) for x in xindices]
        vertical_splitters = [LineString([(xindices[0], y), (xindices[-1], y)]) for y in yindices]
        for splitter in vertical_splitters:
            try:
                outline = MultiPolygon(split(outline, splitter))
            except:
                continue
        for splitter in horizontal_splitters:
            try:
                outline = MultiPolygon(split(outline, splitter))
            except:
                continue
        for tempOutline in list(outline):
            if tempOutline.area < (1 / (10*res**2)):
                continue
            maxArea = max(tempOutline.area,maxArea)
            pxlo, ptlo, pxhi, pthi = tempOutline.bounds
            xstep = (pxhi - pxlo) / subres / 10
            ystep = (pthi - ptlo) / subres / 10
            xs.extend(np.linspace(pxlo + xstep, pxhi - xstep, subres))
            xs.extend(np.linspace(pxhi - xstep, pxlo + xstep, subres))
            ys.extend(np.linspace(pthi - ystep, ptlo + ystep, subres))
            ys.extend(np.linspace(pthi - ystep, ptlo + ystep, subres))

        if len(xs) > 0:
            with open(self.inputFileName, 'w') as inputFile:
                inputFile.write('! Python-generated input file for Thermochimica\n')
                inputFile.write(f'data file         = {self.datafile}\n')
                inputFile.write(f'temperature unit  = {self.tunit}\n')
                inputFile.write(f'pressure unit     = {self.punit}\n')
                inputFile.write(f'mass unit         = \'{self.munit}\'\n')
                inputFile.write( 'nEl               = 3 \n')
                inputFile.write(f'iEl               = {atomic_number_map.index(self.el1)+1} {atomic_number_map.index(self.el2)+1} {atomic_number_map.index(self.el3)+1}\n')
                inputFile.write(f'nCalc             = {len(xs)}\n')
                for i in range(len(xs)):
                    inputFile.write(f'{self.temperature} {self.pressure} {xs[i]} {ys[i]} {1-xs[i]-ys[i]}\n')
            print('Thermochimica calculation initiated.')
            subprocess.run(['./bin/RunCalculationList',self.inputFileName])
            print('Thermochimica calculation finished.')
            self.processPhaseDiagramData()
    def autoRefine2Phase(self,res):
        # Run iteratively
        nIt = 0
        while nIt < 3:
            nIt = nIt + 1
            maxGap = 0
            # Create arrays again with new data
            boundaries = []
            b = []
            for i in range(len(self.p1)):
                # If a miscibility gap label has been used unnecessarily, remove it
                if self.p1[i].find('#') > 0:
                    if not(self.p1[i][0:self.p1[i].find('#')] == self.p2[i]):
                        self.p1[i] = self.p1[i][0:self.p1[i].find('#')]
                if self.p2[i].find('#') > 0:
                    if not(self.p2[i][0:self.p2[i].find('#')] == self.p1[i]):
                        self.p2[i] = self.p2[i][0:self.p2[i].find('#')]
                repeat = False
                for j in range(len(boundaries)):
                    if (boundaries[j][0] == self.p1[i]) and (boundaries[j][1] == self.p2[i]):
                        b.append(j)
                        repeat = True
                if not(repeat):
                    boundaries.append([self.p1[i],self.p2[i]])
                    b.append(len(boundaries)-1)

            # Refine two-phase region density
            xs = []
            ys = []
            for j in range(len(boundaries)):
                inds = [i for i, k in enumerate(b) if k == j]
                if len(inds) < 2:
                    continue
                v1 = np.array([self.x1[inds[1],0],self.x1[inds[1],1],0])
                v2 = np.array([self.x2[inds[1],0],self.x2[inds[1],1],0])
                normal = np.cross(v1-v2,np.array([0,0,1]))
                order = [inds[i] for i, k in sorted(enumerate(self.x1[inds].tolist()), key=lambda coord: coord[1][0]*normal[0]+coord[1][1]*normal[1])]
                x1s = self.x1[order]
                x2s = self.x2[order]
                for i in range(len(x1s)-1):
                    gap = np.linalg.norm(x1s[i]-x1s[i+1])+np.linalg.norm(x2s[i]-x2s[i+1])
                    maxGap = max(gap,maxGap)
                    if gap > 1/res:
                        start = np.average([x1s[i],  x2s[i]],axis=0)
                        end   = np.average([x1s[i+1],x2s[i+1]],axis=0)
                        self.writeInputFile(start[0],end[0],start[1],end[1],np.ceil(gap * res))
                        self.runCalc()

            # Test the minimum difference between points to see if converged
            if maxGap <= 1/res:
                break
    def autoLabel(self):
        self.makeBackup()
        self.sgw.Element('Undo').Update(disabled = False)
        # Make list of boundaries and points belonging to them
        boundaries = []
        b = []
        for i in range(len(self.p1)):
            # If a miscibility gap label has been used unnecessarily, remove it
            if self.p1[i].find('#') > 0:
                if not(self.p1[i][0:self.p1[i].find('#')] == self.p2[i]):
                    self.p1[i] = self.p1[i][0:self.p1[i].find('#')]
            if self.p2[i].find('#') > 0:
                if not(self.p2[i][0:self.p2[i].find('#')] == self.p1[i]):
                    self.p2[i] = self.p2[i][0:self.p2[i].find('#')]
            repeat = False
            for j in range(len(boundaries)):
                if (boundaries[j][0] == self.p1[i]) and (boundaries[j][1] == self.p2[i]):
                    b.append(j)
                    repeat = True
            if not(repeat):
                boundaries.append([self.p1[i],self.p2[i]])
                b.append(len(boundaries)-1)

        # Make list of phases
        phases = []
        for i in range(len(boundaries)):
            repeat1 = False
            repeat2 = False
            for j in range(len(phases)):
                if (boundaries[i][0] == phases[j]):
                    repeat1 = True
                if (boundaries[i][1] == phases[j]):
                    repeat2 = True
            if not(repeat1 or boundaries[i][0].find('#') > 0):
                phases.append(boundaries[i][0])
            if not(repeat2 or boundaries[i][1].find('#') > 0):
                phases.append(boundaries[i][1])

        # label 1-phase regions
        for phase in phases:
            inds1 = [i for i, k in enumerate(self.p1) if k == phase]
            inds2 = [i for i, k in enumerate(self.p2) if k == phase]
            points = []
            if len(inds1) > 0:
                points.extend(self.x1[inds1].tolist())
            if len(inds2) > 0:
                points.extend(self.x2[inds2].tolist())
            for point in self.points1:
                if point[1] == phase:
                    points.append(point[0])
            points = np.array(points)
            average = np.average(points,axis=0)
            self.labels.append([[average[0],average[1]],phase])

        # label 2-phase regions
        for j in range(len(boundaries)):
            inds = [i for i, k in enumerate(b) if k == j]
            if len(inds) < 2:
                continue
            average = (np.average(self.x1[inds],axis=0) + np.average(self.x2[inds],axis=0)) / 2
            self.labels.append([[average[0],average[1]],'+'.join(boundaries[j])])

        # label 3-phase regions
        for point in self.points3:
            average = np.average(point[0],axis=0)
            self.labels.append([[average[0],average[1]],'+'.join(point[1])])
    def makeBackup(self):
        self.backup = CalculationWindow(self.parent, self.datafile, self.nElements, self.elements, False)
        self.backup.datafile = self.datafile
        self.backup.nElements = self.nElements
        self.backup.elements = copy.deepcopy(self.elements)
        self.backup.x1 = copy.deepcopy(self.x1)
        self.backup.x2 = copy.deepcopy(self.x2)
        self.backup.p1 = copy.deepcopy(self.p1)
        self.backup.p2 = copy.deepcopy(self.p2)
        self.backup.points3 = copy.deepcopy(self.points3)
        self.backup.points1 = copy.deepcopy(self.points1)
        self.backup.labels = copy.deepcopy(self.labels)
        self.backup.temperature = self.temperature
        self.backup.pressure = self.pressure
        self.backup.inputFileName = self.inputFileName
        self.backup.outputFileName = self.outputFileName
        self.backup.plotMarker = self.plotMarker
        self.backup.plotColor = self.plotColor
        self.backup.el1 = self.el1
        self.backup.el2 = self.el2
        self.backup.el3 = self.el3
        self.backup.tunit = self.tunit
        self.backup.punit = self.punit
        self.backup.munit = self.munit
        self.backup.exportFormat = self.exportFormat
        self.backup.exportFileName = self.exportFileName
        self.backup.exportDPI = self.exportDPI
        self.backup.resRef = self.resRef
        self.backup.resSmooth = self.resSmooth
        self.backup.tielines = self.tielines
        self.backup.tiegap = self.tiegap
    def activate(self):
        if not self.active:
            self.makeLayout()
            self.sgw = sg.Window(f'Phase Diagram Setup: {os.path.basename(self.datafile)}', self.layout, location = [400,0], finalize=True)
            windowList.append(self)
            self.active = True
            self.parent.children.append(self)
            self.sgw.Element('Refine').Update(disabled = False)
            self.sgw.Element('Auto Refine').Update(disabled = False)
            self.sgw.Element('Auto Smoothen').Update(disabled = False)
            self.sgw.Element('Add Label').Update(disabled = False)
            self.sgw.Element('Auto Label').Update(disabled = False)
            self.sgw.Element('Plot').Update(disabled = False)
            if len(self.labels) > 0:
                self.sgw.Element('Remove Label').Update(disabled = False)
    def makeLayout(self):
        elSelectLayout = [sg.Column([[sg.Text('Element 1')],[sg.Combo(self.elements[:self.nElements],default_value=self.elements[0],key='-el1-')]],vertical_alignment='t'),
                          sg.Column([[sg.Text('Element 2')],[sg.Combo(self.elements[:self.nElements],default_value=self.elements[1],key='-el2-')]],vertical_alignment='t'),
                          sg.Column([[sg.Text('Element 3')],[sg.Combo(self.elements[:self.nElements],default_value=self.elements[2],key='-el3-')]],vertical_alignment='t')]
        xLayout    = [sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstep-',size=(8,1))]],vertical_alignment='t')]
        tempLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperature-',size=(inputSize,1))],
                      [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')]],vertical_alignment='t')
                      ]
        presLayout = [sg.Column([[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(inputSize,1))],
                      [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')]],vertical_alignment='t')
                      ]
        self.layout = [elSelectLayout,xLayout,tempLayout,presLayout,[
            sg.Column([[sg.Button('Run', size = buttonSize)],
                       [sg.Button('Undo', disabled = True, size = buttonSize)],
                       [sg.Exit(size = buttonSize)]],vertical_alignment='t'),
            sg.Column([[sg.Button('Refine', disabled = True, size = buttonSize)],
                       [sg.Button('Auto Refine', disabled = True, size = buttonSize)],
                       [sg.Button('Auto Smoothen', disabled = True, size = buttonSize)]],vertical_alignment='t'),
            sg.Column([[sg.Button('Add Label', disabled = True, size = buttonSize)],
                       [sg.Button('Auto Label', disabled = True, size = buttonSize)],
                       [sg.Button('Remove Label', disabled = True, size = buttonSize)]],vertical_alignment='t'),
            sg.Column([[sg.Button('Plot', disabled = True, size = buttonSize)],
                       [sg.Button('Export Plot', disabled = True, size = buttonSize)],
                       [sg.Button('Plot Settings', size = buttonSize)]],vertical_alignment='t')
            ]]
    def exportPlot(self):
        try:
            self.currentPlot.savefig(f'{self.exportFileName}.{self.exportFormat}', format=self.exportFormat, dpi=self.exportDPI)
        except:
            errorLayout = [[sg.Text('The export failed, try changing plot settings.')],[sg.Button('Continue'), sg.Button('Cancel')]]
            errorWindow = sg.Window('Plot export failed', errorLayout, location = [400,0], finalize=True, keep_on_top = True)
            while True:
                event, values = errorWindow.read(timeout=timeout)
                if event == sg.WIN_CLOSED or event == 'Continue':
                    break
            errorWindow.close()

class RefineWindow:
    def __init__(self, parent, windowLayout):
        self.parent = parent
        windowList.append(self)
        self.sgw = sg.Window('Phase diagram refinement', windowLayout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            self.close()
        elif event =='Refine':
            cancelRun = False
            nxstep = 10
            try:
                tempstep = int(values['-nxstepr-'])
                if tempstep >= 0:
                    nxstep = tempstep
            except:
                pass
            if (float(nxstep)) > 500:
                cancelRun = True
                confirmLayout = [[sg.Text('The selected calculation is large and may take some time.')],[sg.Button('Continue'), sg.Button('Cancel')]]
                confirmWindow = sg.Window('Large calculation confirmation', confirmLayout, location = [400,0], finalize=True, keep_on_top = True)
                while True:
                    event, values = confirmWindow.read(timeout=timeout)
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Continue':
                        cancelRun = False
                        break
                confirmWindow.close()
            xlo1 = 0
            try:
                templo = float(values['-xlor1-'])
                if 0 <= templo <= 1:
                    xlo1 = templo
            except:
                pass
            xhi1 = 1
            try:
                temphi = float(values['-xhir1-'])
                if 0 <= temphi <= 1:
                    xhi1 = temphi
            except:
                pass
            xlo2 = 0
            try:
                templo = float(values['-xlor2-'])
                if 0 <= templo <= 1:
                    xlo2 = templo
            except:
                pass
            xhi2 = 1
            try:
                temphi = float(values['-xhir2-'])
                if 0 <= temphi <= 1:
                    xhi2 = temphi
            except:
                pass
            if not cancelRun:
                self.parent.makeBackup()
                self.parent.sgw.Element('Undo').Update(disabled = False)
                self.parent.writeInputFile(xlo1,xhi1,xlo2,xhi2,nxstep)
                self.parent.runCalc()
                self.parent.makePlot()

class LabelWindow:
    def __init__(self, parent, windowLayout):
        self.parent = parent
        windowList.append(self)
        self.sgw = sg.Window('Add phase label', windowLayout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            self.close()
        elif event =='Add Label':
            try:
                try:
                    x1lab = float(values['-x1lab-'])
                except ValueError:
                    num, den = values['-x1lab-'].split('/')
                    x1lab = float(num)/float(den)
                try:
                    x2lab = float(values['-x2lab-'])
                except ValueError:
                    num, den = values['-x2lab-'].split('/')
                    x2lab = float(num)/float(den)
                if (0 <= x1lab <= 1) and (0 <= x2lab <= (1-x1lab)):
                    self.parent.makeBackup()
                    self.parent.sgw.Element('Undo').Update(disabled = False)
                    self.parent.addLabel(x1lab,x2lab)
                    self.parent.makePlot()
                    self.parent.sgw.Element('Remove Label').Update(disabled = False)
            except:
                pass

class RemoveWindow:
    def __init__(self, parent, windowLayout):
        self.parent = parent
        windowList.append(self)
        self.sgw = sg.Window('Remove phase label', windowLayout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            self.close()
        if event == 'Remove Label(s)':
            self.parent.makeBackup()
            self.parent.sgw.Element('Undo').Update(disabled = False)
            tempLength = len(self.parent.labels)
            for i in reversed(range(tempLength)):
                try:
                    if values['-removeLabel'+str(i)+'-']:
                        del self.parent.labels[i]
                except:
                    continue
            if len(self.parent.labels) == 0:
                self.parent.sgw.Element('Remove Label').Update(disabled = True)
            self.parent.makePlot()
            self.close()

class SettingsWindow:
    def __init__(self, parent, windowLayout):
        self.parent = parent
        windowList.append(self)
        self.sgw = sg.Window('Plot Settings', windowLayout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED:
            self.close()
        elif event == '-mline-':
            self.parent.plotMarker = '-'
        elif event =='-mpoint-':
            self.parent.plotMarker = '.'
        elif event =='-mboth-':
            self.parent.plotMarker = '.-'
        elif event =='-mcolorful-':
            self.parent.plotColor = 'colorful'
        elif event =='-mbland-':
            self.parent.plotColor = 'bland'
        elif event =='Accept':
            self.parent.tielines = values['-tielines-']
            try:
                tempdensity = float(values['-tiedensity-'])
                if 0 < tempdensity < 1000:
                    self.parent.tiegap = 1/tempdensity
            except:
                pass
            try:
                self.parent.exportFileName = str(values['-filename-'])
            except:
                pass
            self.parent.exportFormat = values['-format-']
            try:
                tempDPI = int(values['-dpi-'])
                if tempDPI > 0 > 10000:
                    self.parent.exportDPI = int(values['-dpi-'])
            except:
                pass
            self.parent.makePlot()
            self.close()

if not(os.path.isfile('bin/InputScriptMode')):
    errorLayout = [[sg.Text('No Thermochimica executable available.')],
                   [sg.Text('Either Thermochimica has not been built (run make),')],
                   [sg.Text('or this script was not executed from Thermochimica root directory.')],
                   [sg.Button('Exit')]]
    errorWindow = sg.Window('Thermochimica Error Message', errorLayout, location = [0,0], finalize=True, keep_on_top = True)
    while True:
        event, values = errorWindow.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            break
    errorWindow.close()
    sys.exit()

windowList = []
dataWindow = DataWindow()
while len(windowList) > 0:
    for window in windowList:
        window.read()
