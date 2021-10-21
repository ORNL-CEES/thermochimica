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

timeout = 50
inputSize = 20
buttonSize = 12

# For boundaries of phase regions where both sides have (# phases) < (# elements), only plot points within phaseFractionTol of the boundary
phaseFractionTol = 1e-2
# Below this tolerance, set phase fraction = 0
phaseIncludeTol = 1e-8

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
                    elemLayout = [sg.Column(elem1Layout),sg.Column(elem2Layout)]
                else:
                    elemLayout = [sg.Column(elem1Layout,vertical_alignment='t', scrollable = True, vertical_scroll_only = True, expand_y = True),
                                  sg.Column(elem2Layout,vertical_alignment='t', scrollable = True, vertical_scroll_only = True, expand_y = True)]
                calcLayout = [tempLayout,
                              presLayout,
                              elemLayout,
                              [sg.Text('# of steps')],[sg.Input(key='-nxstep-',size=(8,1))],
                              [sg.Text('Mass unit')],
                              [sg.Combo(['moles'],default_value='moles',key='-munit-')],
                              [sg.Button('Run', size = buttonSize),
                               sg.Button('Refine', disabled = True, size = buttonSize),
                               sg.Button('Add Label', disabled = True, size = buttonSize),
                               sg.Button('Remove Label', disabled = True, size = buttonSize),
                               sg.Exit()]]
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
        self.sgw = sg.Window(f'Phase Diagram Setup: {os.path.basename(datafile)}', windowLayout, location = [400,0], finalize=True)
        self.children = []
        self.inputFileName = 'inputs/pythonPhaseDiagramInput.ti'
        self.outputFileName = 'thermoout.json'
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
                tlo = 300.0
                if values['-temperature-'] != '':
                    tlo = float(values['-temperature-'])
                thi = 1000.0
                if values['-endtemperature-'] != '':
                    thi = float(values['-endtemperature-'])
                ntstep = 10
                if values['-ntstep-'] != '':
                    ntstep = int(values['-ntstep-'])
                nxstep = 10
                self.pressure = 1
                if values['-pressure-'] != '':
                    self.pressure = values['-pressure-']
                self.mint = tlo
                self.maxt = thi
                self.points = []
                self.labels = []
                masses1 = [0.0]*self.nElements
                masses2 = [0.0]*self.nElements
                self.elementsUsed = []
                for i in range(self.nElements):
                    if values['-'+self.elements[i]+'1-'] != '':
                        masses1[i] = float(values['-'+self.elements[i]+'1-'])
                    if values['-'+self.elements[i]+'2-'] != '':
                        masses2[i] = float(values['-'+self.elements[i]+'2-'])
                    if (masses1[i] > 0.0) or (masses2[i] > 0.0):
                        self.elementsUsed.append(self.elements[i])
                for i in reversed(range(self.nElements)):
                    if not self.elements[i] in self.elementsUsed:
                        del masses1[i]
                        del masses2[i]
                self.nElementsUsed = len(self.elementsUsed)
                self.massLabels = ['','']
                for i in range(self.nElementsUsed):
                    if masses1[i] > 0:
                        self.massLabels[0] += self.elementsUsed[i]
                        if masses1[i] != 1:
                            if int(masses1[i]) == masses1[i]:
                                self.massLabels[0] += f'$_{ {int(masses1[i])} }$'
                            else:
                                self.massLabels[0] += f'$_{ {masses1[i]} }$'
                    if masses2[i] > 0:
                        self.massLabels[1] += self.elementsUsed[i]
                        if masses2[i] != 1:
                            if int(masses2[i]) == masses2[i]:
                                self.massLabels[1] += f'$_{ {int(masses2[i])} }$'
                            else:
                                self.massLabels[1] += f'$_{ {masses2[i]} }$'
                sum1 = sum(masses1)
                sum2 = sum(masses2)
                if (sum1 == 0) or (sum2 == 0) or (min(masses1) < 0) or (min(masses2) < 0):
                    alertLayout = [[sg.Text('One or more input masses are invalid.')],[sg.Button('Continue')]]
                    alertWindow = sg.Window('Invalid Mass Alert', alertLayout, location = [400,0], finalize=True)
                    while True:
                        event, values = alertWindow.read(timeout=timeout)
                        if event == sg.WIN_CLOSED or event == 'Continue':
                            break
                    alertWindow.close()
                    return
                masses1 = [mass / sum1 for mass in masses1]
                masses2 = [mass / sum2 for mass in masses2]
                self.plane = np.array([masses1,masses2])
                self.tunit = values['-tunit-']
                self.punit = values['-punit-']
                self.munit = values['-munit-']
                if values['-nxstep-'] != '':
                    nxstep = int(values['-nxstep-'])
                self.runCalc(0,1,nxstep,tlo,thi,ntstep)
                self.processPhaseDiagramData()
                self.makePlot()
                self.sgw.Element('Refine').Update(disabled = False)
                self.sgw.Element('Add Label').Update(disabled = False)
        elif event =='Refine':
            xRefLayout    = [sg.Column([[sg.Text('Start Concentration')],[sg.Input(key='-xlor-',size=(inputSize,1))]],vertical_alignment='t'),
                          sg.Column([[sg.Text('End Concentration')],[sg.Input(key='-xhir-',size=(inputSize,1))]],vertical_alignment='t'),
                          sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstepr-',size=(8,1))]],vertical_alignment='t')]
            tempRefLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperaturer-',size=(inputSize,1))]],vertical_alignment='t'),
                          sg.Column([[sg.Text('End Temperature')],[sg.Input(key='-endtemperaturer-',size=(inputSize,1))]],vertical_alignment='t'),
                          sg.Column([[sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstepr-',size=(8,1))]],vertical_alignment='t')]
            refineLayout = [xRefLayout,tempRefLayout,[sg.Button('Refine'), sg.Button('Cancel')]]
            refineWindow = RefineWindow(self, refineLayout)
            self.children.append(refineWindow)
        elif event =='Add Label':
            xLabLayout  = [[sg.Text(f'{self.massLabels[1].translate({ord(i):None for i in "{}_$"})} Concentration')],[sg.Input(key='-xlab-',size=(inputSize,1))]]
            tLabLayout  = [[sg.Text('Temperature')],[sg.Input(key='-tlab-',size=(inputSize,1))]]
            labelLayout = [xLabLayout,tLabLayout,[sg.Button('Add Label'), sg.Button('Cancel')]]
            labelWindow = LabelWindow(self,labelLayout)
            self.children.append(labelWindow)
        elif event =='Remove Label':
            headingsLayout = [[sg.Text('Label Text',   size = [55,1],justification='left'),
                               sg.Text('Concentration',size = [15,1],justification='center'),
                               sg.Text('Temperature',  size = [15,1],justification='center'),
                               sg.Text('Remove Label?',size = [15,1])]]
            labelListLayout = []
            for i in range(len(self.labels)):
                labelListLayout.append([[sg.Text(self.labels[i][1],size = [55,1],justification='left'),
                                         sg.Text("{:.3f}".format(self.labels[i][0][0]),size = [15,1],justification='center'),
                                         sg.Text("{:.0f}".format(self.labels[i][0][1]),size = [15,1],justification='center'),
                                         sg.Checkbox('',key='-removeLabel'+str(i)+'-',pad=[[40,0],[0,0]])]])
            removeLayout = [headingsLayout,labelListLayout,[sg.Button('Remove Label(s)'), sg.Button('Cancel')]]
            removeWindow = RemoveWindow(self, removeLayout)
            self.children.append(removeWindow)
    def runCalc(self,xlo,xhi,nxstep,tlo,thi,ntstep):
        xs = np.array([np.linspace((1-xlo)*self.plane[0,i] + xlo*self.plane[1,i],(1-xhi)*self.plane[0,i] + xhi*self.plane[1,i],nxstep) for i in range(self.nElementsUsed)]).T
        temps = np.linspace(tlo,thi,ntstep)
        with open(self.inputFileName, 'w') as inputFile:
            inputFile.write('! Python-generated input file for Thermochimica\n')
            inputFile.write(f'data file         = {self.datafile}\n')
            inputFile.write(f'temperature unit         = {self.tunit}\n')
            inputFile.write(f'pressure unit          = {self.punit}\n')
            inputFile.write(f'mass unit          = \'{self.munit}\'\n')
            inputFile.write(f'nEl         = {self.nElementsUsed} \n')
            inputFile.write(f'iEl         = {" ".join([str(atomic_number_map.index(elem)+1) for elem in self.elementsUsed])}\n')
            inputFile.write(f'nCalc       = {len(xs)*len(temps)}\n')
            for t in temps:
                for x in xs:
                    inputFile.write(f'{str(t)} {self.pressure} {" ".join([str(x[i]) for i in range(self.nElementsUsed)])}\n')
        print('Thermochimica calculation initiated.')
        subprocess.run(['./bin/RunCalculationList',self.inputFileName])
        print('Thermochimica calculation finished.')
    def processPhaseDiagramData(self):
        f = open(self.outputFileName,)
        data = json.load(f)
        f.close()
        if list(data.keys())[0] != '1':
            print('Output does not contain data series')
            exit()
        for i in list(data.keys()):
            self.mint = min(self.mint,data[i]['temperature'])
            self.maxt = max(self.maxt,data[i]['temperature'])
            nPhases = 0
            for phaseType in ['solution phases','pure condensed phases']:
                for phaseName in list(data[i][phaseType].keys()):
                    if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                        nPhases += 1
            if nPhases == 2:
                boundPhases = []
                phaseCompositions = np.zeros([nPhases,self.nElementsUsed])
                iPhase = 0
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                            boundPhases.append(phaseName)
                            for k in range(self.nElementsUsed):
                                if self.elementsUsed[k] in data[i][phaseType][phaseName]['elements'].keys():
                                    phaseCompositions[iPhase,k] = data[i][phaseType][phaseName]['elements'][self.elementsUsed[k]]["mole fraction of phase by element"]
                            iPhase += 1
                crossNorms = [np.linalg.norm(np.cross(phaseCompositions[k] - self.plane[0],self.plane[1] - phaseCompositions[k])) for k in range(nPhases)]
                if max(crossNorms) < phaseIncludeTol:
                    boundComps = [np.linalg.norm(phaseCompositions[k] - self.plane[0])/np.linalg.norm(self.plane[1] - self.plane[0]) for k in range(nPhases)]
                    self.points.append([data[i]['temperature'],boundComps[0],boundPhases])
                    self.points.append([data[i]['temperature'],boundComps[1],boundPhases])
                    continue
            if nPhases == self.nElementsUsed:
                allPhases = []
                phaseComps = []
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                            allPhases.append(phaseName)
                            tempComp = []
                            for element in self.elementsUsed:
                                tempComp.append(data[i][phaseType][phaseName]['elements'][element]['mole fraction of phase by element'])
                            phaseComps.append(tempComp)
                # Loop over possible phase zone intersections with plane of interest
                for j in range(self.nElementsUsed):
                    # Make list of phases on a (nElements-1) dimensional face through omission
                    omitComps = phaseComps.copy()
                    omitComps.remove(phaseComps[j])
                    omitPhase = allPhases.copy()
                    omitPhase.remove(allPhases[j])
                    lineComps = []
                    for k in range(self.nElementsUsed - 2):
                        lineComps.append([omitComps[0],omitComps[k+1]])
                    # Calculate intersection of this face with our plane of interest
                    intersect = self.line_intersection(lineComps)
                    # Check that the intersection is within the valid bounds
                    intSum = sum(intersect[1:])
                    intTest = intSum <= 1
                    for test in range(self.nElementsUsed - 1):
                        intTest = intTest and (0 <= intersect[test]) and (intersect[test] <= 1)
                    if intTest:
                        if intSum == 1:
                            # If we are on the far boundary, the first phase is not included
                            omitPhase.remove(omitPhase[0])
                        for k in range(self.nElementsUsed - 2):
                            # Check all the other boundaries
                            if intersect[k+1] == 0:
                                # If none of this is used, it is not included
                                omitPhase.remove(omitPhase[k+1])
                        self.points.append([data[i]['temperature'],intersect[0],omitPhase])
            elif nPhases > 1:
                boundPhases = []
                skipPoint = False
                phaseMoleSum = 0
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        phaseMoleSum += data[i][phaseType][phaseName]['moles']
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if data[i][phaseType][phaseName]['moles'] > phaseIncludeTol:
                            boundPhases.append(phaseName)
                            if phaseFractionTol < data[i][phaseType][phaseName]['moles']/phaseMoleSum < (1-phaseFractionTol):
                                skipPoint = True
                                break
                    if skipPoint:
                        break
                if skipPoint:
                    boundPhases = []
                    continue
                tempComp = np.zeros(self.nElementsUsed)
                for e in range(len(self.elementsUsed)):
                    if self.elementsUsed[e] in data[i]['elements'].keys():
                        tempComp[e] = data[i]['elements'][self.elementsUsed[e]]['moles']
                boundComps = np.linalg.norm(tempComp-self.plane[0])/np.linalg.norm(self.plane[1]-self.plane[0])
                self.points.append([data[i]['temperature'],boundComps,boundPhases])
    def makePlot(self):
        boundaries = []
        b = []
        for point in self.points:
            repeat = False
            for j in range(len(boundaries)):
                thisMatch = True
                if not (len(point[2]) == len(boundaries[j])):
                    continue
                for phase in point[2]:
                    if not (phase in boundaries[j]):
                        thisMatch = False
                        break
                if thisMatch:
                    b.append(j)
                    repeat = True
            if not(repeat):
                b.append(len(boundaries))
                boundaries.append(point[2])

        # Start figure
        fig = plt.figure()
        plt.ion()
        ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])

        for j in range(len(boundaries)):
            inds = [i for i, k in enumerate(b) if k == j]
            if len(inds) < 2:
                continue
            plotPoints = np.array([[self.points[i][1],self.points[i][0]] for i, k in enumerate(b) if k == j])
            ax.plot(plotPoints[:,0],plotPoints[:,1],'.')

        ax.set_xlim(0,1)
        ax.set_ylim(self.mint,self.maxt)
        title = " $-$ ".join(self.massLabels)
        ax.set_title(r'{0} phase diagram'.format(title))
        ax.set_xlabel(r'Mole fraction {0}'.format(self.massLabels[1]))
        ax.set_ylabel(r'Temperature [K]')
        for lab in self.labels:
            plt.text(float(lab[0][0]),float(lab[0][1]),lab[1], ha="center")
        plt.show()
        plt.pause(0.001)
    def addLabel(self,xlab,tlab):
        self.runCalc(xlab,xlab,1,tlab,tlab,1)
        f = open(self.outputFileName,)
        data = json.load(f)
        f.close()
        if list(data.keys())[0] != '1':
            print('Output does not contain data series')
            exit()
        labelName = []
        for phaseName in list(data['1']['solution phases'].keys()):
            if (data['1']['solution phases'][phaseName]['moles'] > phaseIncludeTol):
                labelName.append(phaseName)
        for phaseName in list(data['1']['pure condensed phases'].keys()):
            if (data['1']['pure condensed phases'][phaseName]['moles'] > phaseIncludeTol):
                labelName.append(phaseName)
        self.labels.append([[xlab,tlab],'+'.join(labelName)])
    def line_intersection(self, lines):
        l1 = np.array(self.plane)
        ls = np.array(lines)
        def diff(mu):
            sum = l1[0] + mu[0]*(l1[1] - l1[0])
            for i in range(self.nElementsUsed - 2):
                sum -= ls[i][0] + mu[i+1]*(ls[i][1] - ls[i][0])
            return sum
        return scipy.optimize.least_squares(diff, [0.5 for i in range(self.nElementsUsed-1)]).x

class RefineWindow():
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
            ntstep = 10
            if values['-ntstepr-'] != '':
                ntstep = int(values['-ntstepr-'])
            nxstep = 10
            if values['-nxstepr-'] != '':
                nxstep = int(values['-nxstepr-'])
            if (float(ntstep) * float(nxstep)) > 50000:
                cancelRun = True
                confirmLayout = [[sg.Text('The selected calculation is large and may take some time.')],[sg.Button('Continue'), sg.Button('Cancel')]]
                confirmWindow = sg.Window('Large calculation confirmation', confirmLayout, location = [400,0], finalize=True)
                while True:
                    event, values = confirmWindow.read(timeout=timeout)
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Continue':
                        cancelRun = False
                        break
                confirmWindow.close()
            tlo = float(values['-temperaturer-'])
            thi = float(values['-endtemperaturer-'])
            xhi = float(values['-xhir-'])
            xlo = float(values['-xlor-'])
            if not cancelRun:
                self.parent.runCalc(xlo,xhi,nxstep,tlo,thi,ntstep)
                self.parent.processPhaseDiagramData()
                self.parent.makePlot()

class LabelWindow():
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
            xlab = float(values['-xlab-'])
            tlab = float(values['-tlab-'])
            self.parent.addLabel(xlab,tlab)
            self.parent.processPhaseDiagramData()
            self.parent.makePlot()
            self.parent.sgw.Element('Remove Label').Update(disabled = False)

class RemoveWindow():
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
        if event == 'Remove Label(s)':
            tempLength = len(self.parent.labels)
            for i in reversed(range(tempLength)):
                if values['-removeLabel'+str(i)+'-']:
                    del self.parent.labels[i]
            if len(self.parent.labels) == 0:
                self.parent.sgw.Element('Remove Label').Update(disabled = True)
            self.parent.makePlot()
            self.close()

windowList = []
dataWindow = DataWindow()
while len(windowList) > 0:
    for window in windowList:
        window.read()
