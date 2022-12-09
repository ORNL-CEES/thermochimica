import PySimpleGUI as sg
import math
import numpy as np
import thermoToolsGUI
import sys
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.geometry import LineString
from shapely.geometry import GeometryCollection
from shapely.prepared import prep
from shapely.ops import split
from functools import reduce
import operator
import thermoTools

# Classes

class pdPoint:
    def __init__(self,elements,temperature,concentration,phases,phaseConcentrations,energy,iterations):
        self.t = temperature
        self.runConcentration = concentration
        self.phaseConcentrations = phaseConcentrations
        self.phases = phases
        self.details = f'Temperature = {temperature:6.2f}\n'
        for elem,conc in zip(elements,self.runConcentration):
            self.details = self.details + f'Moles of {elem} = {conc:9.8f}\n'
        for phase,conc in zip(self.phases,self.phaseConcentrations):
            self.details = self.details + f'{phase} at {conc:5.4f}\n'
        self.details = self.details + f'Integral Gibbs Energy = {energy:.2f}\n'
        self.details = self.details + f'Number of GEM iterations = {iterations}'
        self.suppressed = False

class InspectWindow:
    def __init__(self,parent,endMember2,phases,windowList):
        self.parent = parent
        windowList.append(self)
        self.windowList = windowList
        dataColumn = [
            [sg.Text('Data Points')],
            [sg.Listbox(values=[], enable_events=True, size=(30, 50), key='-dataList-')]
        ]
        outputColumn = [
            [sg.Text('Calculation Details')],
            [sg.Multiline(key='-details-', size=(50,10), no_scrollbar=True)],
            [sg.Text(key = '-status-')],
            [sg.Button('Toggle Active/Suppressed Status', disabled = True)],
            [sg.Text('Filter points', font='underline')],
            [sg.Text('Temperature Range:')],
            [sg.Input(key='-tfilterlow-',size=(thermoToolsGUI.inputSize,1)),sg.Input(key='-tfilterhi-',size=(thermoToolsGUI.inputSize,1))],
            [sg.Text(f'{endMember2} Concentration Range:')],
            [sg.Input(key='-xfilterlow-',size=(thermoToolsGUI.inputSize,1)),sg.Input(key='-xfilterhi-',size=(thermoToolsGUI.inputSize,1))],
            [sg.Text('Contains Phases:')],
            [sg.Combo(['']+phases, key = '-pfilter1-'),sg.Combo(['']+phases, key = '-pfilter2-')],
            [sg.Text('Active/Suppressed Status:')],
            [sg.Combo(['','Active','Suppressed'], key = '-activefilter-')],
            [sg.Button('Apply Filter')]
        ]
        self.data = [[i, f'{point.t:6.2f} K {point.phaseConcentrations[0]:4.3f} {point.phaseConcentrations[1]:4.3f}'] for i,point in enumerate(self.parent.calculation.pdPoints)]
        self.sgw = sg.Window('Data inspection',
            [[sg.Pane([
                sg.Column(dataColumn, element_justification='l', expand_x=True, expand_y=True),
                sg.Column(outputColumn, element_justification='c', expand_x=True, expand_y=True)
            ], orientation='h', k='-PANE-')]],
            location = [0,0], finalize=True)
        self.sgw['-dataList-'].update(self.data)
        self.children = []
        self.index = -1
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in self.windowList:
            self.windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event == '-dataList-':
            self.index = values['-dataList-'][0][0]
            self.point = self.parent.calculation.pdPoints[self.index]
            self.sgw['-details-'].update(self.point.details)
            self.sgw['Toggle Active/Suppressed Status'].update(disabled = False)
            self.sgw['-status-'].update(f'{"Suppressed" if self.point.suppressed else "Active"}')
        elif event == 'Toggle Active/Suppressed Status':
            if self.index >= 0:
                self.point.suppressed = not(self.point.suppressed)
                self.parent.macro.append(f'macroPD.suppressed[{self.index}] = not(macroPD.suppressed[{self.index}])')
                self.sgw['-status-'].update(f'{"Suppressed" if self.point.suppressed else "Active"}')
        elif event == 'Apply Filter':
            tlo = -np.Inf
            thi  = np.Inf
            xlo = -np.Inf
            xhi  = np.Inf
            try:
                tlo = float(values['-tfilterlow-'])
            except:
                pass
            try:
                thi = float(values['-tfilterhi-'])
            except:
                pass
            try:
                xlo = float(values['-xfilterlow-'])
            except:
                pass
            try:
                xhi = float(values['-xfilterhi-'])
            except:
                pass
            self.data = []
            for i,p in enumerate(self.parent.calculation.pdPoints):
                # Check temperature
                tfilt = tlo <= p.t and thi >= p.t
                # Check concentration
                xfilt = (xlo <= p.phaseConcentrations[0] and xhi >= p.phaseConcentrations[0]) or (xlo <= p.phaseConcentrations[1] and xhi >= p.phaseConcentrations[1])
                # Check phases present
                pfilt = (values['-pfilter1-'] == '' or values['-pfilter1-'] == p.phases[0] or values['-pfilter1-'] == p.phases[1]) and (values['-pfilter2-'] == '' or values['-pfilter2-'] == p.phases[0] or values['-pfilter2-'] == p.phases[1])
                # Check active/suppressed status
                afilt = (values['-activefilter-'] == '') or ((values['-activefilter-'] == 'Suppressed') == p.suppressed)
                # If all filters pass, add to display list
                if tfilt and xfilt and pfilt and afilt:
                    self.data.append([i, f'{p.t:6.2f} K {p.phaseConcentrations[0]:4.3f} {p.phaseConcentrations[1]:4.3f}'])
            self.sgw['-dataList-'].update(self.data)

class RefineWindow:
    def __init__(self, parent,windowList):
        self.parent = parent
        windowList.append(self)
        self.windowList = windowList
        xRefLayout    = [sg.Column([[sg.Text('Start Concentration')],[sg.Input(key='-xlor-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t'),
                         sg.Column([[sg.Text('End Concentration')],[sg.Input(key='-xhir-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t'),
                         sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstepr-',size=(8,1))]],vertical_alignment='t')]
        tempRefLayout = [sg.Column([[sg.Text('Minimum Temperature')],[sg.Input(key='-temperaturer-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t'),
                         sg.Column([[sg.Text('Maximum Temperature')],[sg.Input(key='-endtemperaturer-',size=(thermoToolsGUI.inputSize,1))]],vertical_alignment='t'),
                         sg.Column([[sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstepr-',size=(8,1))]],vertical_alignment='t')]
        refineLayout = [xRefLayout,tempRefLayout,[sg.Button('Refine'), sg.Button('Cancel')]]
        self.sgw = sg.Window('Phase diagram refinement', refineLayout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in self.windowList:
            self.windowList.remove(self)
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
                confirmWindow = sg.Window('Large calculation confirmation', confirmLayout, location = [400,0], finalize=True, keep_on_top = True)
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
            # refine x-coords are going to come in scaled to axis (for compound endmembers)
            if self.parent.calculation.compoundScale:
                if xlo > 0:
                    xlo = 1/(1+((1-xlo)/xlo)*(self.parent.calculation.sum1/self.parent.calculation.sum2))
                if xhi > 0:
                    xhi = 1/(1+((1-xhi)/xhi)*(self.parent.calculation.sum1/self.parent.calculation.sum2))
            if not cancelRun:
                self.parent.calculation.makeBackup()
                self.parent.sgw.Element('Undo').Update(disabled = False)
                self.parent.calculation.writeInputFile(xlo,xhi,nxstep,tlo,thi,ntstep)
                self.parent.calculation.runCalc()
                self.parent.calculation.makePlot()
                self.parent.macro.append('macroPD.makeBackup()')
                self.parent.macro.append(f'macroPD.writeInputFile({xlo},{xhi},{nxstep},{tlo},{thi},{ntstep})')
                self.parent.macro.append('macroPD.runCalc()')

class LabelWindow:
    def __init__(self, parent, windowList):
        self.parent = parent
        xLabLayout  = [[sg.Text(f'{self.parent.calculation.massLabels[1].translate({ord(i):None for i in "{}_$"})} Concentration')],[sg.Input(key='-xlab-',size=(thermoToolsGUI.inputSize,1))]]
        tLabLayout  = [[sg.Text('Temperature')],[sg.Input(key='-tlab-',size=(thermoToolsGUI.inputSize,1))]]
        labelLayout = [xLabLayout,tLabLayout,[sg.Button('Add Label'), sg.Button('Cancel')]]
        self.sgw = sg.Window('Add phase label', labelLayout, location = [400,0], finalize=True)
        windowList.append(self)
        self.windowList = windowList
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in self.windowList:
            self.windowList.remove(self)
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
    def __init__(self, parent, windowList):
        self.parent = parent
        windowList.append(self)
        self.windowList = windowList
        headingsLayout = [[sg.Text('Label Text',   size = [55,1],justification='left'),
                           sg.Text('Concentration',size = [15,1],justification='center'),
                           sg.Text('Temperature',  size = [15,1],justification='center'),
                           sg.Text('Remove Label?',size = [15,1])]]
        labelListLayout = []
        for i in range(len(self.parent.calculation.labels)):
            labelListLayout.append([[sg.Text(self.parent.calculation.labels[i][1],size = [55,1],justification='left'),
                                     sg.Text("{:.3f}".format(float(self.parent.calculation.labels[i][0][0])),size = [15,1],justification='center'),
                                     sg.Text("{:.0f}".format(float(self.parent.calculation.labels[i][0][1])),size = [15,1],justification='center'),
                                     sg.Checkbox('',key='-removeLabel'+str(i)+'-',pad=[[40,0],[0,0]])]])
        removeLayout = [headingsLayout,labelListLayout,[sg.Button('Remove Label(s)'), sg.Button('Cancel')]]
        self.sgw = sg.Window('Remove phase label', removeLayout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in self.windowList:
            self.windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            self.close()
        if event == 'Remove Label(s)':
            self.parent.calculation.makeBackup()
            self.parent.macro.append('macroPD.makeBackup()')
            self.parent.sgw.Element('Undo').Update(disabled = False)
            tempLength = len(self.parent.calculation.labels)
            for i in reversed(range(tempLength)):
                try:
                    if values['-removeLabel'+str(i)+'-']:
                        del self.parent.calculation.labels[i]
                        self.parent.macro.append(f'del macroPD.labels[{i}]')
                except KeyError:
                    # If a new label was created since this window was opened, this will occur
                    continue
            if len(self.parent.calculation.labels) == 0:
                self.parent.sgw.Element('Remove Label').Update(disabled = True)
            self.parent.calculation.makePlot()
            self.close()

class SettingsWindow:
    def __init__(self, parent, windowList, onePhaseLabelsDisabled = False):
        self.parent = parent
        windowList.append(self)
        self.windowList = windowList
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
                           sg.Checkbox('Loaded Diagram', default=self.parent.calculation.showLoaded, key='-showLoaded-')],
                          [sg.Text('Auto-Label Settings:')],
                          [sg.Checkbox('1-Phase Regions', default=self.parent.calculation.label1phase, key='-label1phase-', disabled=onePhaseLabelsDisabled),
                           sg.Checkbox('2-Phase Regions', default=self.parent.calculation.label2phase, key='-label2phase-')],
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
        self.parent.macro.append(f'macroPD.showLoaded = {self.parent.calculation.showLoaded}')
        self.parent.macro.append(f'macroPD.label1phase = {self.parent.calculation.label1phase}')
        self.parent.macro.append(f'macroPD.label2phase = {self.parent.calculation.label2phase}')
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in self.windowList:
            self.windowList.remove(self)
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
            self.parent.calculation.showLoaded = values['-showLoaded-']
            self.parent.calculation.label1phase = values['-label1phase-']
            self.parent.calculation.label2phase = values['-label2phase-']
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

# Functions

def runMacro(calc):
    if 'macroPhaseDiagram' in sys.modules:
        del sys.modules['macroPhaseDiagram']
    import macroPhaseDiagram
    calc.calculation = macroPhaseDiagram.macroPD
    calc.calculation.active = True
    calc.calculation.interactivePlot = True
    enableButtons(calc)
    
def enableButtons(calc):
    calc.sgw.Element('Refine').Update(disabled = False)
    calc.sgw.Element('Auto Refine').Update(disabled = False)
    calc.sgw.Element('Auto Smoothen').Update(disabled = False)
    calc.sgw.Element('Add Label').Update(disabled = False)
    calc.sgw.Element('Auto Label').Update(disabled = False)
    calc.sgw.Element('Plot').Update(disabled = False)
    calc.sgw.Element('Undo').Update(disabled = False)
    calc.sgw.Element('Inspect').Update(disabled = False)
    calc.sgw.Element('Export Diagram Data').Update(disabled = False)
    calc.sgw.Element('Export Plot').Update(disabled = False)

def phaseBoundaries(calc):
    calc.boundaries = []
    calc.phases = []
    calc.b = []
    for point in calc.pdPoints:
        # If a miscibility gap label has been used unnecessarily, remove it
        if point.phases[0].find('#') > 0:
            if not(point.phases[0][0:point.phases[0].find('#')] == point.phases[1]):
                point.phases[0] = point.phases[0][0:point.phases[0].find('#')]
        if point.phases[1].find('#') > 0:
            if not(point.phases[1][0:point.phases[1].find('#')] == point.phases[0]):
                point.phases[1] = point.phases[1][0:point.phases[1].find('#')]
        repeat = False
        if point.suppressed:
            calc.b.append(-1)
        else:
            for j in range(len(calc.boundaries)):
                if (calc.boundaries[j][0] == point.phases[0]) and (calc.boundaries[j][1] == point.phases[1]):
                    calc.b.append(j)
                    repeat = True
            if not(repeat):
                calc.boundaries.append([point.phases[0],point.phases[1]])
                calc.b.append(len(calc.boundaries)-1)

    for i in range(len(calc.boundaries)):
        repeat1 = False
        repeat2 = False
        for j in range(len(calc.phases)):
            if (calc.boundaries[i][0] == calc.phases[j]):
                repeat1 = True
            if (calc.boundaries[i][1] == calc.phases[j]):
                repeat2 = True
        if not(repeat1 or calc.boundaries[i][0].find('#') > 0):
            calc.phases.append(calc.boundaries[i][0])
        if not(repeat2 or calc.boundaries[i][1].find('#') > 0):
            calc.phases.append(calc.boundaries[i][1])

    calc.congruentFound = [False for _ in range(len(calc.phases))]
    for j in range(len(calc.boundaries)):
        inds = [i for i, k in enumerate(calc.b) if k == j]
        if len(inds) < 2:
            continue
        ttt = [calc.pdPoints[i].t for i in inds]
        x1t = [calc.pdPoints[i].phaseConcentrations[0] for i in inds]
        x2t = [calc.pdPoints[i].phaseConcentrations[1] for i in inds]
        if x1t[0] > x2t[0]:
            dir = True
        else:
            dir = False
        extraBound = []
        for i in range(len(ttt)):
            if (x1t[i] > x2t[i]) != dir:
                # for miscibility gap, just flip them
                if calc.boundaries[j][0].find('#') > 0 or calc.boundaries[j][1].find('#') > 0:
                    temp = calc.pdPoints[inds[i]].phaseConcentrations[0]
                    calc.pdPoints[inds[i]].phaseConcentrations[0] = calc.pdPoints[inds[i]].phaseConcentrations[1]
                    calc.pdPoints[inds[i]].phaseConcentrations[1] = temp
                else:
                    extraBound.append(i)
        if len(extraBound):
            calc.congruentFound[calc.phases.index(calc.boundaries[j][0])] = True
            calc.congruentFound[calc.phases.index(calc.boundaries[j][1])] = True
            calc.boundaries.append(calc.boundaries[j])
            for k in extraBound:
                calc.b[inds[k]] = len(calc.boundaries)-1

    for j in range(len(calc.boundaries)):
        inds = [i for i, k in enumerate(calc.b) if k == j]
        if len(inds) < 2:
            continue
        ttt = [calc.pdPoints[i].t for i in inds]
        x1t = [calc.pdPoints[i].phaseConcentrations[0] for i in inds]
        x2t = [calc.pdPoints[i].phaseConcentrations[1] for i in inds]
        loc = False
        firstLoc = True
        for i in range(1,len(ttt)):
            if np.sqrt((ttt[i] - ttt[i-1])**2 + ((calc.maxt - calc.mint)*(x1t[i] - x1t[i-1]))**2 + ((calc.maxt - calc.mint)*(x2t[i] - x2t[i-1]))**2) > calc.gapLimit:
                loc = not(loc)
                if firstLoc:
                    calc.boundaries.append(calc.boundaries[j])
                    firstLoc = False
            if loc:
                calc.b[inds[i]] = len(calc.boundaries)-1

def autoRefine(calc,res,endpoints,useDiagramEdges=True,maxIts=4):
    nIt = 0
    while nIt < maxIts:
        nIt = nIt + 1
        maxArea = 0
        phaseBoundaries(calc)

        phasePolyPoints = [[] for _ in range(len(calc.phases))]

        if useDiagramEdges:
            for j in range(len(calc.x0data[1])):
                try:
                    i = calc.phases.index(calc.x0data[0][j])
                    phasePolyPoints[i].append([[0,calc.x0data[1][j]]])
                    phasePolyPoints[i].append([[0,calc.x0data[2][j]]])
                except:
                    continue
            for j in range(len(calc.x1data[1])):
                try:
                    i = calc.phases.index(calc.x1data[0][j])
                    phasePolyPoints[i].append([[1,calc.x1data[1][j]]])
                    phasePolyPoints[i].append([[1,calc.x1data[2][j]]])
                except:
                    continue

        # plot 2-phase region boundaries
        for j in range(len(calc.boundaries)):
            polygonPoints = []
            inds = [i for i, k in enumerate(calc.b) if k == j]
            if len(inds) < 2:
                continue
            ttt = [calc.pdPoints[i].t for i in inds]
            x1t = [calc.pdPoints[i].phaseConcentrations[0] for i in inds]
            x2t = [calc.pdPoints[i].phaseConcentrations[1] for i in inds]
            for i in range(len(inds)):
                polygonPoints.append([x1t[i],ttt[i]])
            for i in reversed(range(len(inds))):
                polygonPoints.append([x2t[i],ttt[i]])
            phaseOutline = Polygon(polygonPoints).buffer(0)
            calc.outline = calc.outline.buffer(0) - phaseOutline
            for i in range(len(calc.phases)):
                if calc.boundaries[j][0] == calc.phases[i]:
                    phasePolyPoints[i].append(polygonPoints[:len(inds)])
                if calc.boundaries[j][1] == calc.phases[i]:
                    phasePolyPoints[i].append(list(reversed(polygonPoints))[:len(inds)])

        for i in range(len(calc.phases)):
            if calc.congruentFound[i]:
                print(f'Warning: congruent phase transformation found, auto refine will skip {calc.phases[i]}')
                continue
            segcenters = []
            if len(phasePolyPoints[i]) < 2:
                continue
            for j in range(len(phasePolyPoints[i])):
                segcenters.append(tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), phasePolyPoints[i][j]), [len(phasePolyPoints[i][j])] * 2)))
            center = tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), segcenters), [len(segcenters)] * 2))
            sortcenters = sorted(segcenters, key=lambda coord: (-135 - math.degrees(math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360)
            sortedPolyPoints = []
            for j in range(len(phasePolyPoints[i])):
                k = segcenters.index(sortcenters[j])
                if sortcenters[j][1] > sortcenters[j-1][1]:
                    for l in range(len(phasePolyPoints[i][k])):
                        sortedPolyPoints.append(phasePolyPoints[i][k][l])
                else:
                    for l in reversed(range(len(phasePolyPoints[i][k]))):
                        sortedPolyPoints.append(phasePolyPoints[i][k][l])
            if len(sortedPolyPoints) > 2:
                phaseOutline = Polygon(sortedPolyPoints).buffer(0)
                try:
                    calc.outline = calc.outline - phaseOutline
                except:
                    pass

        xs = []
        ys = []
        subres = int(np.ceil(np.sqrt(res)))
        try:
            oxlo, otlo, oxhi, othi = calc.outline.bounds
        except:
            continue
        xindices = np.linspace(oxlo, oxhi, subres)
        yindices = np.linspace(otlo, othi, subres)
        horizontal_splitters = [LineString([(x, yindices[0]), (x, yindices[-1])]) for x in xindices]
        vertical_splitters = [LineString([(xindices[0], y), (xindices[-1], y)]) for y in yindices]
        # If the outline contains non-polygon shapes (like lines) it will be a GeometryCollection instead
        # and we need to remove those non-polygon shapes so it can be a MultiPolygon again
        if isinstance(calc.outline,GeometryCollection):
            calc.outline = MultiPolygon([shape for shape in list(calc.outline.geoms) if isinstance(shape,Polygon)])
        for splitter in vertical_splitters:
            try:
                calc.outline = MultiPolygon(split(calc.outline, splitter))
            except:
                continue
        for splitter in horizontal_splitters:
            try:
                calc.outline = MultiPolygon(split(calc.outline, splitter))
            except:
                continue
        for tempOutline in list(calc.outline.geoms):
            if (tempOutline.area / (calc.maxt-calc.mint)) < (1 / (10*res**2)):
                continue
            maxArea = max(tempOutline.area / (calc.maxt-calc.mint),maxArea)
            pxlo, ptlo, pxhi, pthi = tempOutline.bounds
            xstep = (pxhi - pxlo) / subres / 10
            ystep = (pthi - ptlo) / subres / 10
            xs.extend(np.linspace(pxlo + xstep, pxhi - xstep, subres))
            xs.extend(np.linspace(pxhi - xstep, pxlo + xstep, subres))
            ys.extend(np.linspace(pthi - ystep, ptlo + ystep, subres))
            ys.extend(np.linspace(pthi - ystep, ptlo + ystep, subres))

        if len(xs) > 0:
            calcList = []
            for i in range(len(xs)):
                concentration = endpoints[0]*(1-xs[i]) + endpoints[1]*xs[i]
                calcItem = [ys[i],calc.pressure]
                calcItem.extend(concentration)
                calcList.append(calcItem)
            thermoTools.WriteRunCalculationList(calc.inputFileName,calc.datafile,calc.elementsUsed,calcList,tunit=calc.tunit,punit=calc.punit,munit=calc.munit,printMode=0,fuzzyStoichiometry=calc.fuzzy,gibbsMinCheck=calc.fuzzy)
            print('Thermochimica calculation initiated.')
            thermoTools.RunRunCalculationList(calc.inputFileName)
            print('Thermochimica calculation finished.')
            calc.processPhaseDiagramData()

        # Test the minimum subgrid region area to see if converged
        if maxArea < 1 / (10*res**2):
            break
        elif any(calc.congruentFound):
            break

def autoRefine2Phase(calc,res,endpoints,maxIts=4):
    phaseBoundaries(calc)
    # Expand two-phase regions
    tres = (calc.maxt-calc.mint)/res
    xs = []
    ys = []
    for j in range(len(calc.boundaries)):
        inds = [i for i, k in enumerate(calc.b) if k == j]
        if len(inds) < 2:
            continue
        ttt = [calc.pdPoints[i].t for i in inds]
        x1t = [calc.pdPoints[i].phaseConcentrations[0] for i in inds]
        x2t = [calc.pdPoints[i].phaseConcentrations[1] for i in inds]
        tbound = max(calc.mint,ttt[0]-tres*3)
        for k in np.arange(tbound,max(calc.mint,ttt[0]-tres/3),tres/3):
            ys.append(k)
            xs.append((x1t[0] + x2t[0])/2)
            ys.append(k)
            xs.append((0.99*x1t[0] + 0.01*x2t[0]))
            ys.append(k)
            xs.append((0.01*x1t[0] + 0.99*x2t[0]))
        tbound = min(calc.maxt,ttt[-1]+tres*3)
        for k in np.arange(min(calc.maxt,ttt[-1]+tres/3),tbound,tres/3):
            ys.append(k)
            xs.append((x1t[-1] + x2t[-1])/2)
            ys.append(k)
            xs.append((0.99*x1t[-1] + 0.01*x2t[-1]))
            ys.append(k)
            xs.append((0.01*x1t[-1] + 0.99*x2t[-1]))

    if len(xs) > 0:
        calcList = []
        for i in range(len(xs)):
            concentration = endpoints[0]*(1-xs[i]) + endpoints[1]*xs[i]
            calcItem = [ys[i],calc.pressure]
            calcItem.extend(concentration)
            calcList.append(calcItem)
        thermoTools.WriteRunCalculationList(calc.inputFileName,calc.datafile,calc.elementsUsed,calcList,tunit=calc.tunit,punit=calc.punit,munit=calc.munit,printMode=0,fuzzyStoichiometry=calc.fuzzy,gibbsMinCheck=calc.fuzzy)
        print('Thermochimica calculation initiated.')
        thermoTools.RunRunCalculationList(calc.inputFileName)
        print('Thermochimica calculation finished.')
        calc.processPhaseDiagramData()

    nIt = 0
    while nIt < maxIts:
        nIt = nIt + 1
        maxGap = 0
        phaseBoundaries(calc)
        # Refine two-phase region density
        xs = []
        ys = []
        for j in range(len(calc.boundaries)):
            inds = [i for i, k in enumerate(calc.b) if k == j]
            if len(inds) < 2:
                continue
            ttt = [calc.pdPoints[i].t for i in inds]
            x1t = [calc.pdPoints[i].phaseConcentrations[0] for i in inds]
            x2t = [calc.pdPoints[i].phaseConcentrations[1] for i in inds]
            for i in range(len(ttt)-1):
                gap = np.sqrt(((ttt[i]-ttt[i+1])/(calc.maxt-calc.mint))**2+(x1t[i]-x1t[i+1])**2+(x2t[i]-x2t[i+1])**2)
                maxGap = max(gap,maxGap)
                if gap > 1/res:
                    step = tres*((ttt[i+1] - ttt[i])/(calc.maxt - calc.mint))/gap
                    try:
                        for k in np.arange(ttt[i] + step,ttt[i+1]-step,step):
                            ys.append(k)
                            progk = (k - ttt[i]) / (ttt[i+1] - ttt[i])
                            xs.append(progk * (x1t[i+1] + x2t[i+1]) / 2 + (1 - progk) * (x1t[i] +  x2t[i]) / 2)
                    except:
                        continue

        if len(xs) > 0:
            calcList = []
            for i in range(len(xs)):
                concentration = endpoints[0]*(1-xs[i]) + endpoints[1]*xs[i]
                calcItem = [ys[i],calc.pressure]
                calcItem.extend(concentration)
                calcList.append(calcItem)
            thermoTools.WriteRunCalculationList(calc.inputFileName,calc.datafile,calc.elementsUsed,calcList,tunit=calc.tunit,punit=calc.punit,munit=calc.munit,printMode=0,fuzzyStoichiometry=calc.fuzzy,gibbsMinCheck=calc.fuzzy)
            print('Thermochimica calculation initiated.')
            thermoTools.RunRunCalculationList(calc.inputFileName)
            print('Thermochimica calculation finished.')
            calc.processPhaseDiagramData()

        # Test the minimum difference between points to see if converged
        if maxGap <= 1/res:
            break
    calc.gapLimit = 3*tres

def autoLabel(calc):
    phaseBoundaries(calc)

    phasePolyPoints = [[] for _ in range(len(calc.phases))]

    if calc.label1phase:
        for j in range(len(calc.x0data[1])):
            try:
                i = calc.phases.index(calc.x0data[0][j])
                phasePolyPoints[i].append([[0,calc.x0data[1][j]]])
                phasePolyPoints[i].append([[0,calc.x0data[2][j]]])
            except:
                continue
        for j in range(len(calc.x1data[1])):
            try:
                i = calc.phases.index(calc.x1data[0][j])
                phasePolyPoints[i].append([[1,calc.x1data[1][j]]])
                phasePolyPoints[i].append([[1,calc.x1data[2][j]]])
            except:
                continue

    # plot 2-phase region boundaries
    for j in range(len(calc.boundaries)):
        polygonPoints = []
        inds = [i for i, k in enumerate(calc.b) if k == j]
        if len(inds) < 2:
            continue
        ttt = np.array([calc.pdPoints[i].t for i in inds])
        x1t = np.array([calc.pdPoints[i].phaseConcentrations[0] for i in inds])
        x2t = np.array([calc.pdPoints[i].phaseConcentrations[1] for i in inds])
        for i in range(len(inds)):
            polygonPoints.append([x1t[i],ttt[i]])
        for i in reversed(range(len(inds))):
            polygonPoints.append([x2t[i],ttt[i]])
        phaseOutline = Polygon(polygonPoints)#.buffer(0)
        center = list(phaseOutline.centroid.coords)[0]
        if calc.label2phase:
            calc.labels.append([[center[0],center[1]-calc.tshift],'+'.join(calc.boundaries[j])])
        for i in range(len(calc.phases)):
            if calc.boundaries[j][0] == calc.phases[i]:
                phasePolyPoints[i].append(polygonPoints[:len(inds)])
            if calc.boundaries[j][1] == calc.phases[i]:
                phasePolyPoints[i].append(list(reversed(polygonPoints))[:len(inds)])
    if calc.label1phase:
        for i in range(len(calc.phases)):
            if calc.congruentFound[i]:
                print(f'Warning: congruent phase transformation found, auto label will skip {calc.phases[i]}')
                continue
            segcenters = []
            if len(phasePolyPoints[i]) < 2:
                continue
            for j in range(len(phasePolyPoints[i])):
                segcenters.append(tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), phasePolyPoints[i][j]), [len(phasePolyPoints[i][j])] * 2)))
            center = tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), segcenters), [len(segcenters)] * 2))
            calc.labels.append([[center[0],center[1]-calc.tshift],calc.phases[i]])
