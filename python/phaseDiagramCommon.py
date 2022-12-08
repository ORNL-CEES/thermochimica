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
