import json
import matplotlib.pyplot as plt
import numpy as np
import math
import subprocess
import copy
from itertools import cycle
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.geometry import MultiPoint
from shapely.geometry import LineString
from shapely.geometry import GeometryCollection
from shapely.prepared import prep
from shapely.ops import split
from functools import reduce
import operator
import csv
import thermoTools

phaseIncludeTol = 1e-8

class diagram:
    def __init__(self, datafile, active, interactivePlot):
        self.datafile = datafile
        self.active = active
        self.interactivePlot = interactivePlot
        self.mint = 1e5
        self.maxt = 0
        self.ts = np.empty([0])
        self.x1 = np.empty([0])
        self.x2 = np.empty([0])
        self.p1 = []
        self.p2 = []
        self.el1 = ''
        self.el2 = ''
        self.tunit = 'K'
        self.punit = 'atm'
        self.munit = 'moles'
        self.tshift = 0
        self.x0data = [[],[],[]]
        self.x1data = [[],[],[]]
        self.labels = []
        self.outline = MultiPolygon([])
        self.pressure = 1
        self.inputFileName = 'inputs/pythonPhaseDiagramInput.ti'
        self.outputFileName = 'outputs/thermoout.json'
        self.plotMarker = '-'
        self.plotColor = 'colorful'
        if self.active:
            self.backup = diagram(self.datafile, False, self.interactivePlot)
        else:
            self.backup = []
        self.currentPlot = []
        self.exportFormat = 'png'
        self.exportFileName = 'thermochimicaPhaseDiagram'
        self.exportDPI = 300
        self.resRef = 7
        self.resSmooth = 7
        self.gapLimit = np.Inf
        self.figureList = []
        self.boundaries = []
        self.phases = []
        self.b = []
        self.congruentFound = [False for i in range(len(self.phases))]
        self.label1phase = True
        self.label2phase = True
        self.experimentalData = []
        self.experimentNames = []
        self.experimentColor = 'bland'
        self.showExperiment = True
        self.pointDetails = []
        self.pointIndex = np.empty([0])
        self.suppressed = []
        self.loadedDiagram = []
        self.loaded = False
        self.showLoaded = True
        self.saveDataName = 'savedDiagram'
    def run(self,ntstep,nxstep,pressure,tunit,punit,xlo,xhi,tlo,thi,el1,el2,munit):
        self.pressure = pressure
        self.tunit = tunit
        self.punit = punit
        self.munit = munit
        self.el1 = el1
        self.el2 = el2
        self.writeInputFile(xlo,xhi,nxstep,tlo,thi,ntstep)
        self.ts = np.empty([0])
        self.x1 = np.empty([0])
        self.x2 = np.empty([0])
        self.p1 = []
        self.p2 = []
        self.x0data = [[],[],[]]
        self.x1data = [[],[],[]]
        # Check temperature unit for shift
        if self.tunit == 'K':
            self.tshift = 0
        elif self.tunit == 'C':
            self.tshift = 273.15
        # Initially reverse shift (opposite at plot)
        self.mint = tlo + self.tshift
        self.maxt = thi + self.tshift
        self.labels = []
        self.resRef = 7
        self.resSmooth = 7
        self.gapLimit = (self.maxt - self.mint) / 2
        self.experimentalData = []
        self.experimentNames = []
        self.pointDetails = []
        self.pointIndex = np.empty([0])
        self.suppressed = []
        self.loadedDiagram = []
        self.loaded = False
        self.saveDataName = 'savedDiagram'
        self.backup = diagram(self.datafile, False, self.interactivePlot)
        for fig in self.figureList:
            plt.close(fig=fig)
        self.runCalc()
        self.outline = MultiPolygon([Polygon([[0,self.mint], [0, self.maxt], [1, self.maxt], [1, self.mint]])])
    def refinery(self):
        self.refineLimit(0,(self.maxt-self.mint)/(self.resRef**2)/10)
        self.refineLimit(1,(self.maxt-self.mint)/(self.resRef**2)/10)
        self.autoRefine(self.resRef**2)
        self.resRef += 1
    def autoSmooth(self):
        self.autoRefine2Phase(self.resSmooth**2)
        self.resSmooth += 1
    def processPhaseDiagramData(self):
        f = open(self.outputFileName,)
        try:
            data = json.load(f)
            f.close()
        except:
            f.close()
            print('Data load failed, aborting phase diagram update')
            return
        if list(data.keys())[0] != '1':
            print('Output does not contain data series')
            exit()
        ts = self.ts.tolist()
        x1 = self.x1.tolist()
        x2 = self.x2.tolist()
        pointIndex = self.pointIndex.tolist()
        for i in list(data.keys()):
            try:
                self.mint = min(self.mint,data[i]['temperature'])
                self.maxt = max(self.maxt,data[i]['temperature'])
            except:
                continue
            nPhases = 0
            for phaseType in ['solution phases','pure condensed phases']:
                for phaseName in list(data[i][phaseType].keys()):
                    if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                        nPhases += 1
            if nPhases == 2:
                ts.append(data[i]['temperature'])
                boundPhases = []
                boundComps = []
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                            boundPhases.append(phaseName)
                            boundComps.append(data[i][phaseType][phaseName]['elements'][self.el2]['mole fraction of phase by element'])
                x1.append(boundComps[0])
                x2.append(boundComps[1])
                pointIndex.append(len(pointIndex))
                self.p1.append(boundPhases[0])
                self.p2.append(boundPhases[1])
                self.pointDetails.append(f'Temperature = {data[i]["temperature"]:6.2f}\nMoles of {self.el1} = {data[i]["elements"][self.el1]["moles"]:9.8f}\nMoles of {self.el2} = {data[i]["elements"][self.el2]["moles"]:9.8f}\nPhase 1 = {boundPhases[0]} at {boundComps[0]:5.4f} moles {self.el2}\nPhase 2 = {boundPhases[1]} at {boundComps[1]:5.4f} moles {self.el2}\nIntegral Gibbs Energy = {data[i]["integral Gibbs energy"]:.2f}\nNumber of GEM iterations = {data[i]["GEM iterations"]}')
                self.suppressed.append(False)
            elif nPhases == 1:
                if not(self.el2 in list(data[i]['elements'].keys())):
                    for phaseType in ['solution phases','pure condensed phases']:
                        for phaseName in list(data[i][phaseType].keys()):
                            if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                                pname = phaseName
                    if not(pname in self.x0data[0]):
                        self.x0data[0].append(pname)
                        self.x0data[1].append(data[i]['temperature'])
                        self.x0data[2].append(data[i]['temperature'])
                    pindex = self.x0data[0].index(pname)
                    self.x0data[1][pindex] = min(self.x0data[1][pindex],data[i]['temperature'])
                    self.x0data[2][pindex] = max(self.x0data[2][pindex],data[i]['temperature'])
                elif float(data[i]['elements'][self.el2]['moles']) == 1:
                    for phaseType in ['solution phases','pure condensed phases']:
                        for phaseName in list(data[i][phaseType].keys()):
                            if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                                pname = phaseName
                    if not(pname in self.x1data[0]):
                        self.x1data[0].append(pname)
                        self.x1data[1].append(data[i]['temperature'])
                        self.x1data[2].append(data[i]['temperature'])
                    pindex = self.x1data[0].index(pname)
                    self.x1data[1][pindex] = min(self.x1data[1][pindex],data[i]['temperature'])
                    self.x1data[2][pindex] = max(self.x1data[2][pindex],data[i]['temperature'])

        # Sort data here instead of repeatedly later
        self.ts = np.array(ts)
        self.x1 = np.array(x1)
        self.x2 = np.array(x2)
        self.pointIndex = np.array(pointIndex)
        sindex  = np.argsort(self.ts)
        self.ts = self.ts[sindex]
        self.x1 = self.x1[sindex]
        self.x2 = self.x2[sindex]
        self.pointIndex = self.pointIndex[sindex]
        self.p1 = [self.p1[i] for i in sindex]
        self.p2 = [self.p2[i] for i in sindex]

        if len(self.x0data[1]) > 1:
            x0sort = [i[0] for i in sorted(enumerate(self.x0data[1]), key=lambda x:x[1])]
            phaseOrder = []
            for i in x0sort:
                phaseOrder.append(self.x0data[0][i])
            xtemp = [[],[],[]]
            xtemp[0] = phaseOrder
            xtemp[1] = sorted(self.x0data[1])
            xtemp[2] = sorted(self.x0data[2])
            self.x0data = xtemp
        if len(self.x1data[1]) > 1:
            x1sort = [i[0] for i in sorted(enumerate(self.x1data[1]), key=lambda x:x[1])]
            phaseOrder = []
            for i in x1sort:
                phaseOrder.append(self.x1data[0][i])
            xtemp = [[],[],[]]
            xtemp[0] = phaseOrder
            xtemp[1] = sorted(self.x1data[1])
            xtemp[2] = sorted(self.x1data[2])
            self.x1data = xtemp
    def runCalc(self):
        print('Thermochimica calculation initiated.')
        subprocess.run(['./bin/PhaseDiagramDataGen',self.inputFileName])
        print('Thermochimica calculation finished.')
        self.processPhaseDiagramData()
    def phaseBoundaries(self):
        self.boundaries = []
        self.phases = []
        self.b = []
        for i in range(len(self.p1)):
            # If a miscibility gap label has been used unnecessarily, remove it
            if self.p1[i].find('#') > 0:
                if not(self.p1[i][0:self.p1[i].find('#')] == self.p2[i]):
                    self.p1[i] = self.p1[i][0:self.p1[i].find('#')]
            if self.p2[i].find('#') > 0:
                if not(self.p2[i][0:self.p2[i].find('#')] == self.p1[i]):
                    self.p2[i] = self.p2[i][0:self.p2[i].find('#')]
            repeat = False
            if self.suppressed[self.pointIndex[i]]:
                self.b.append(-1)
            else:
                for j in range(len(self.boundaries)):
                    if (self.boundaries[j][0] == self.p1[i]) and (self.boundaries[j][1] == self.p2[i]):
                        self.b.append(j)
                        repeat = True
                if not(repeat):
                    self.boundaries.append([self.p1[i],self.p2[i]])
                    self.b.append(len(self.boundaries)-1)

        for i in range(len(self.boundaries)):
            repeat1 = False
            repeat2 = False
            for j in range(len(self.phases)):
                if (self.boundaries[i][0] == self.phases[j]):
                    repeat1 = True
                if (self.boundaries[i][1] == self.phases[j]):
                    repeat2 = True
            if not(repeat1 or self.boundaries[i][0].find('#') > 0):
                self.phases.append(self.boundaries[i][0])
            if not(repeat2 or self.boundaries[i][1].find('#') > 0):
                self.phases.append(self.boundaries[i][1])

        self.congruentFound = [False for i in range(len(self.phases))]
        for j in range(len(self.boundaries)):
            inds = [i for i, k in enumerate(self.b) if k == j]
            if len(inds) < 2:
                continue
            ttt = self.ts[inds]
            x1t = self.x1[inds]
            x2t = self.x2[inds]
            if x1t[0] > x2t[0]:
                dir = True
            else:
                dir = False
            extraBound = []
            for i in range(len(ttt)):
                if (x1t[i] > x2t[i]) != dir:
                    # for miscibility gap, just flip them
                    if self.boundaries[j][0].find('#') > 0 or self.boundaries[j][1].find('#') > 0:
                        temp = self.x1[inds[i]]
                        self.x1[inds[i]] = self.x2[inds[i]]
                        self.x2[inds[i]] = temp
                    else:
                        extraBound.append(i)
            if len(extraBound):
                self.congruentFound[self.phases.index(self.boundaries[j][0])] = True
                self.congruentFound[self.phases.index(self.boundaries[j][1])] = True
                self.boundaries.append(self.boundaries[j])
                for k in extraBound:
                    self.b[inds[k]] = len(self.boundaries)-1

        for j in range(len(self.boundaries)):
            inds = [i for i, k in enumerate(self.b) if k == j]
            if len(inds) < 2:
                continue
            ttt = self.ts[inds]
            x1t = self.x1[inds]
            x2t = self.x2[inds]
            loc = False
            firstLoc = True
            for i in range(1,len(ttt)):
                if np.sqrt((ttt[i] - ttt[i-1])**2 + ((self.maxt - self.mint)*(x1t[i] - x1t[i-1]))**2 + ((self.maxt - self.mint)*(x2t[i] - x2t[i-1]))**2) > self.gapLimit:
                    loc = not(loc)
                    if firstLoc:
                        self.boundaries.append(self.boundaries[j])
                        firstLoc = False
                if loc:
                    self.b[inds[i]] = len(self.boundaries)-1
    def makePlot(self):
        self.phaseBoundaries()
        # Start figure
        fig = plt.figure()
        if self.interactivePlot:
            plt.ion()
        ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])

        bEdgeLine = [[False,False] for i in range(len(self.boundaries))]
        # Plot along x=0 and x=1 self.boundaries (this is the worst code I've ever written)
        for j in range(len(self.x0data[1])):
            if not self.x0data[0][j] in self.phases:
                continue
            i = self.phases.index(self.x0data[0][j])
            if j > 0:
                # ax.plot(0,self.x0data[1][j],'kv')
                match = []
                for k in range(len(self.boundaries)):
                    if (self.x0data[0][j] in self.boundaries[k]) and (self.x0data[0][j-1] in self.boundaries[k]):
                        inds = [i for i, l in enumerate(self.b) if l == k]
                        if len(inds) < 2:
                            continue
                        bind = self.boundaries[k].index(self.x0data[0][j])
                        if bind == 0:
                            minj = np.argmin(np.array(self.x1)[inds])
                            length = (0 - np.array(self.x1)[inds][minj])**2 + (self.x0data[1][j] - np.array(self.ts)[inds][minj])**2
                            match.append([length,k,np.array(self.x1)[inds][minj],np.array(self.ts)[inds][minj]])
                        elif bind == 1:
                            minj = np.argmin(np.array(self.x2)[inds])
                            length = (0 - np.array(self.x2)[inds][minj])**2 + (self.x0data[1][j] - np.array(self.ts)[inds][minj])**2
                            match.append([length,k,np.array(self.x2)[inds][minj],np.array(self.ts)[inds][minj]])
                if len(match) > 0:
                    match = np.array(match)
                    matchind = np.argmin(match[:,0])
                    k = int(match[matchind,1])
                    inds = [i for i, l in enumerate(self.b) if l == k]
                    ax.plot([0,match[matchind,2]],[self.x0data[1][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                    if match[matchind,3] == np.min(np.array(self.ts)[inds]):
                        bEdgeLine[k][0] = True
                    if match[matchind,3] == np.max(np.array(self.ts)[inds]):
                        bEdgeLine[k][1] = True
            if j < len(self.x0data[1]) - 1:
                # ax.plot(0,self.x0data[2][j],'k^')
                match = []
                for k in range(len(self.boundaries)):
                    if (self.x0data[0][j] in self.boundaries[k]) and (self.x0data[0][j+1] in self.boundaries[k]):
                        inds = [i for i, l in enumerate(self.b) if l == k]
                        if len(inds) < 2:
                            continue
                        bind = self.boundaries[k].index(self.x0data[0][j])
                        if bind == 0:
                            minj = np.argmin(np.array(self.x1)[inds])
                            length = (0 - np.array(self.x1)[inds][minj])**2 + (self.x0data[2][j] - np.array(self.ts)[inds][minj])**2
                            match.append([length,k,np.array(self.x1)[inds][minj],np.array(self.ts)[inds][minj]])
                        elif bind == 1:
                            minj = np.argmin(np.array(self.x2)[inds])
                            length = (0 - np.array(self.x2)[inds][minj])**2 + (self.x0data[2][j] - np.array(self.ts)[inds][minj])**2
                            match.append([length,k,np.array(self.x2)[inds][minj],np.array(self.ts)[inds][minj]])
                if len(match) > 0:
                    match = np.array(match)
                    matchind = np.argmin(match[:,0])
                    k = int(match[matchind,1])
                    inds = [i for i, l in enumerate(self.b) if l == k]
                    ax.plot([0,match[matchind,2]],[self.x0data[2][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                    if match[matchind,3] == np.min(np.array(self.ts)[inds]):
                        bEdgeLine[k][0] = True
                    if match[matchind,3] == np.max(np.array(self.ts)[inds]):
                        bEdgeLine[k][1] = True
        for j in range(len(self.x1data[1])):
            if not self.x1data[0][j] in self.phases:
                continue
            i = self.phases.index(self.x1data[0][j])
            if j > 0:
                # ax.plot(1,self.x1data[1][j],'kv')
                match = []
                for k in range(len(self.boundaries)):
                    if (self.x1data[0][j] in self.boundaries[k]) and (self.x1data[0][j-1] in self.boundaries[k]):
                        inds = [i for i, l in enumerate(self.b) if l == k]
                        if len(inds) < 2:
                            continue
                        bind = self.boundaries[k].index(self.x1data[0][j])
                        if bind == 0:
                            maxj = np.argmax(np.array(self.x1)[inds])
                            length = (1 - np.array(self.x1)[inds][maxj])**2 + (self.x1data[1][j] - np.array(self.ts)[inds][maxj])**2
                            match.append([length,k,np.array(self.x1)[inds][maxj],np.array(self.ts)[inds][maxj]])
                        elif bind == 1:
                            maxj = np.argmax(np.array(self.x2)[inds])
                            length = (1 - np.array(self.x2)[inds][maxj])**2 + (self.x1data[1][j] - np.array(self.ts)[inds][maxj])**2
                            match.append([length,k,np.array(self.x2)[inds][maxj],np.array(self.ts)[inds][maxj]])
                if len(match) > 0:
                    match = np.array(match)
                    matchind = np.argmin(match[:,0])
                    k = int(match[matchind,1])
                    inds = [i for i, l in enumerate(self.b) if l == k]
                    ax.plot([1,match[matchind,2]],[self.x1data[1][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                    if match[matchind,3] == np.min(np.array(self.ts)[inds]):
                        bEdgeLine[k][0] = True
                    if match[matchind,3] == np.max(np.array(self.ts)[inds]):
                        bEdgeLine[k][1] = True
            if j < len(self.x1data[1]) - 1:
                # ax.plot(1,self.x1data[2][j],'k^')
                match = []
                for k in range(len(self.boundaries)):
                    if (self.x1data[0][j] in self.boundaries[k]) and (self.x1data[0][j+1] in self.boundaries[k]):
                        inds = [i for i, l in enumerate(self.b) if l == k]
                        if len(inds) < 2:
                            continue
                        bind = self.boundaries[k].index(self.x1data[0][j])
                        if bind == 0:
                            maxj = np.argmax(np.array(self.x1)[inds])
                            length = (1 - np.array(self.x1)[inds][maxj])**2 + (self.x1data[2][j] - np.array(self.ts)[inds][maxj])**2
                            match.append([length,k,np.array(self.x1)[inds][maxj],np.array(self.ts)[inds][maxj]])
                        elif bind == 1:
                            maxj = np.argmax(np.array(self.x2)[inds])
                            length = (1 - np.array(self.x2)[inds][maxj])**2 + (self.x1data[2][j] - np.array(self.ts)[inds][maxj])**2
                            match.append([length,k,np.array(self.x2)[inds][maxj],np.array(self.ts)[inds][maxj]])
                if len(match) > 0:
                    match = np.array(match)
                    matchind = np.argmin(match[:,0])
                    k = int(match[matchind,1])
                    inds = [i for i, l in enumerate(self.b) if l == k]
                    ax.plot([1,match[matchind,2]],[self.x1data[2][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                    if match[matchind,3] == np.min(np.array(self.ts)[inds]):
                        bEdgeLine[k][0] = True
                    if match[matchind,3] == np.max(np.array(self.ts)[inds]):
                        bEdgeLine[k][1] = True

        # plot 2-phase region boundaries
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(self.boundaries))))
        for j in range(len(self.boundaries)):
            if self.plotColor == 'colorful':
                c = next(color)
            else:
                c = 'k'
            inds = [i for i, k in enumerate(self.b) if k == j]
            if len(inds) < 2:
                continue
            ttt = self.ts[inds]
            x1t = self.x1[inds]
            x2t = self.x2[inds]
            ax.plot(x1t,ttt-self.tshift,self.plotMarker,c=c)
            ax.plot(x2t[::-1],ttt[::-1]-self.tshift,self.plotMarker,c=c)
            minj = np.argmin(ttt)
            maxj = np.argmax(ttt)
            # plot invariant temperatures
            if (ttt[minj] > self.mint) and not(bEdgeLine[j][0]):
                ax.plot([x1t[minj],x2t[minj]],[ttt[minj]-self.tshift,ttt[minj]-self.tshift],self.plotMarker,c=c)
            if (ttt[maxj] < self.maxt) and not(bEdgeLine[j][1]):
                ax.plot([x1t[maxj],x2t[maxj]],[ttt[maxj]-self.tshift,ttt[maxj]-self.tshift],self.plotMarker,c=c)

        # Plot experimental data
        if self.showExperiment:
            markerList = cycle(['o','v','<','^','s'])
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(self.experimentalData))))
            for e in range(len(self.experimentalData)):
                if self.experimentColor == 'colorful':
                    c = next(color)
                else:
                    c = 'k'
                m = next(markerList)
                ax.plot(self.experimentalData[e][:,0],self.experimentalData[e][:,1],m,c=c,label=self.experimentNames[e])

        ax.set_xlim(0,1)
        ax.set_ylim(self.mint-self.tshift,self.maxt-self.tshift)
        ax.set_title(str(self.el1) + ' + ' + str(self.el2) + ' binary phase diagram')
        ax.set_xlabel('Mole fraction ' + str(self.el2))
        tunit_display = r"$^\circ$" if self.tunit == "C" else ""
        ax.set_ylabel(f'Temperature [{tunit_display}{self.tunit}]')
        if len(self.experimentalData) > 0 and self.showExperiment:
            ax.legend(loc=0)
        for lab in self.labels:
            plt.text(float(lab[0][0]),float(lab[0][1]),lab[1], ha='center')

        # Plot loaded phase diagram
        if self.loaded and self.showLoaded:
            bEdgeLine = [[False,False] for i in range(len(self.loadedDiagram.boundaries))]
            # Plot along x=0 and x=1 self.boundaries (this is the worst code I've ever written)
            for j in range(len(self.loadedDiagram.x0data[1])):
                if not self.loadedDiagram.x0data[0][j] in self.loadedDiagram.phases:
                    continue
                i = self.loadedDiagram.phases.index(self.loadedDiagram.x0data[0][j])
                if j > 0:
                    # ax.plot(0,self.loadedDiagram.x0data[1][j],'kv')
                    match = []
                    for k in range(len(self.loadedDiagram.boundaries)):
                        if (self.loadedDiagram.x0data[0][j] in self.loadedDiagram.boundaries[k]) and (self.loadedDiagram.x0data[0][j-1] in self.loadedDiagram.boundaries[k]):
                            inds = [i for i, l in enumerate(self.loadedDiagram.b) if l == k]
                            if len(inds) < 2:
                                continue
                            bind = self.loadedDiagram.boundaries[k].index(self.loadedDiagram.x0data[0][j])
                            if bind == 0:
                                minj = np.argmin(np.array(self.loadedDiagram.x1)[inds])
                                length = (0 - np.array(self.loadedDiagram.x1)[inds][minj])**2 + (self.loadedDiagram.x0data[1][j] - np.array(self.loadedDiagram.ts)[inds][minj])**2
                                match.append([length,k,np.array(self.loadedDiagram.x1)[inds][minj],np.array(self.loadedDiagram.ts)[inds][minj]])
                            elif bind == 1:
                                minj = np.argmin(np.array(self.loadedDiagram.x2)[inds])
                                length = (0 - np.array(self.loadedDiagram.x2)[inds][minj])**2 + (self.loadedDiagram.x0data[1][j] - np.array(self.loadedDiagram.ts)[inds][minj])**2
                                match.append([length,k,np.array(self.loadedDiagram.x2)[inds][minj],np.array(self.loadedDiagram.ts)[inds][minj]])
                    if len(match) > 0:
                        match = np.array(match)
                        matchind = np.argmin(match[:,0])
                        k = int(match[matchind,1])
                        inds = [i for i, l in enumerate(self.loadedDiagram.b) if l == k]
                        ax.plot([0,match[matchind,2]],[self.loadedDiagram.x0data[1][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                        if match[matchind,3] == np.min(np.array(self.loadedDiagram.ts)[inds]):
                            bEdgeLine[k][0] = True
                        if match[matchind,3] == np.max(np.array(self.loadedDiagram.ts)[inds]):
                            bEdgeLine[k][1] = True
                if j < len(self.loadedDiagram.x0data[1]) - 1:
                    # ax.plot(0,self.loadedDiagram.x0data[2][j],'k^')
                    match = []
                    for k in range(len(self.loadedDiagram.boundaries)):
                        if (self.loadedDiagram.x0data[0][j] in self.loadedDiagram.boundaries[k]) and (self.loadedDiagram.x0data[0][j+1] in self.loadedDiagram.boundaries[k]):
                            inds = [i for i, l in enumerate(self.loadedDiagram.b) if l == k]
                            if len(inds) < 2:
                                continue
                            bind = self.loadedDiagram.boundaries[k].index(self.loadedDiagram.x0data[0][j])
                            if bind == 0:
                                minj = np.argmin(np.array(self.loadedDiagram.x1)[inds])
                                length = (0 - np.array(self.loadedDiagram.x1)[inds][minj])**2 + (self.loadedDiagram.x0data[2][j] - np.array(self.loadedDiagram.ts)[inds][minj])**2
                                match.append([length,k,np.array(self.loadedDiagram.x1)[inds][minj],np.array(self.loadedDiagram.ts)[inds][minj]])
                            elif bind == 1:
                                minj = np.argmin(np.array(self.loadedDiagram.x2)[inds])
                                length = (0 - np.array(self.loadedDiagram.x2)[inds][minj])**2 + (self.loadedDiagram.x0data[2][j] - np.array(self.loadedDiagram.ts)[inds][minj])**2
                                match.append([length,k,np.array(self.loadedDiagram.x2)[inds][minj],np.array(self.loadedDiagram.ts)[inds][minj]])
                    if len(match) > 0:
                        match = np.array(match)
                        matchind = np.argmin(match[:,0])
                        k = int(match[matchind,1])
                        inds = [i for i, l in enumerate(self.loadedDiagram.b) if l == k]
                        ax.plot([0,match[matchind,2]],[self.loadedDiagram.x0data[2][j]-self.tshift,match[matchind,3]]-self.tshift,'k-')
                        if match[matchind,3] == np.min(np.array(self.loadedDiagram.ts)[inds]):
                            bEdgeLine[k][0] = True
                        if match[matchind,3] == np.max(np.array(self.loadedDiagram.ts)[inds]):
                            bEdgeLine[k][1] = True
            for j in range(len(self.loadedDiagram.x1data[1])):
                if not self.loadedDiagram.x1data[0][j] in self.loadedDiagram.phases:
                    continue
                i = self.loadedDiagram.phases.index(self.loadedDiagram.x1data[0][j])
                if j > 0:
                    # ax.plot(1,self.loadedDiagram.x1data[1][j],'kv')
                    match = []
                    for k in range(len(self.loadedDiagram.boundaries)):
                        if (self.loadedDiagram.x1data[0][j] in self.loadedDiagram.boundaries[k]) and (self.loadedDiagram.x1data[0][j-1] in self.loadedDiagram.boundaries[k]):
                            inds = [i for i, l in enumerate(self.loadedDiagram.b) if l == k]
                            if len(inds) < 2:
                                continue
                            bind = self.loadedDiagram.boundaries[k].index(self.loadedDiagram.x1data[0][j])
                            if bind == 0:
                                maxj = np.argmax(np.array(self.loadedDiagram.x1)[inds])
                                length = (1 - np.array(self.loadedDiagram.x1)[inds][maxj])**2 + (self.loadedDiagram.x1data[1][j] - np.array(self.loadedDiagram.ts)[inds][maxj])**2
                                match.append([length,k,np.array(self.loadedDiagram.x1)[inds][maxj],np.array(self.loadedDiagram.ts)[inds][maxj]])
                            elif bind == 1:
                                maxj = np.argmax(np.array(self.loadedDiagram.x2)[inds])
                                length = (1 - np.array(self.loadedDiagram.x2)[inds][maxj])**2 + (self.loadedDiagram.x1data[1][j] - np.array(self.loadedDiagram.ts)[inds][maxj])**2
                                match.append([length,k,np.array(self.loadedDiagram.x2)[inds][maxj],np.array(self.loadedDiagram.ts)[inds][maxj]])
                    if len(match) > 0:
                        match = np.array(match)
                        matchind = np.argmin(match[:,0])
                        k = int(match[matchind,1])
                        inds = [i for i, l in enumerate(self.loadedDiagram.b) if l == k]
                        ax.plot([1,match[matchind,2]],[self.loadedDiagram.x1data[1][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                        if match[matchind,3] == np.min(np.array(self.loadedDiagram.ts)[inds]):
                            bEdgeLine[k][0] = True
                        if match[matchind,3] == np.max(np.array(self.loadedDiagram.ts)[inds]):
                            bEdgeLine[k][1] = True
                if j < len(self.loadedDiagram.x1data[1]) - 1:
                    # ax.plot(1,self.loadedDiagram.x1data[2][j],'k^')
                    match = []
                    for k in range(len(self.loadedDiagram.boundaries)):
                        if (self.loadedDiagram.x1data[0][j] in self.loadedDiagram.boundaries[k]) and (self.loadedDiagram.x1data[0][j+1] in self.loadedDiagram.boundaries[k]):
                            inds = [i for i, l in enumerate(self.loadedDiagram.b) if l == k]
                            if len(inds) < 2:
                                continue
                            bind = self.loadedDiagram.boundaries[k].index(self.loadedDiagram.x1data[0][j])
                            if bind == 0:
                                maxj = np.argmax(np.array(self.loadedDiagram.x1)[inds])
                                length = (1 - np.array(self.loadedDiagram.x1)[inds][maxj])**2 + (self.loadedDiagram.x1data[2][j] - np.array(self.loadedDiagram.ts)[inds][maxj])**2
                                match.append([length,k,np.array(self.loadedDiagram.x1)[inds][maxj],np.array(self.loadedDiagram.ts)[inds][maxj]])
                            elif bind == 1:
                                maxj = np.argmax(np.array(self.loadedDiagram.x2)[inds])
                                length = (1 - np.array(self.loadedDiagram.x2)[inds][maxj])**2 + (self.loadedDiagram.x1data[2][j] - np.array(self.loadedDiagram.ts)[inds][maxj])**2
                                match.append([length,k,np.array(self.loadedDiagram.x2)[inds][maxj],np.array(self.loadedDiagram.ts)[inds][maxj]])
                    if len(match) > 0:
                        match = np.array(match)
                        matchind = np.argmin(match[:,0])
                        k = int(match[matchind,1])
                        inds = [i for i, l in enumerate(self.loadedDiagram.b) if l == k]
                        ax.plot([1,match[matchind,2]],[self.loadedDiagram.x1data[2][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                        if match[matchind,3] == np.min(np.array(self.loadedDiagram.ts)[inds]):
                            bEdgeLine[k][0] = True
                        if match[matchind,3] == np.max(np.array(self.loadedDiagram.ts)[inds]):
                            bEdgeLine[k][1] = True

            # plot 2-phase region boundaries
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(self.loadedDiagram.boundaries))))
            for j in range(len(self.loadedDiagram.boundaries)):
                c = 'k'
                inds = [i for i, k in enumerate(self.loadedDiagram.b) if k == j]
                if len(inds) < 2:
                    continue
                ttt = self.loadedDiagram.ts[inds]
                x1t = self.loadedDiagram.x1[inds]
                x2t = self.loadedDiagram.x2[inds]
                ax.plot(x1t,ttt-self.tshift,'--',c=c)
                ax.plot(x2t[::-1],ttt[::-1]-self.tshift,'--',c=c)
                minj = np.argmin(ttt)
                maxj = np.argmax(ttt)
                # plot invariant temperatures
                if (ttt[minj] > self.loadedDiagram.mint) and not(bEdgeLine[j][0]):
                    ax.plot([x1t[minj],x2t[minj]],[ttt[minj]-self.tshift,ttt[minj]-self.tshift],'--',c=c)
                if (ttt[maxj] < self.loadedDiagram.maxt) and not(bEdgeLine[j][1]):
                    ax.plot([x1t[maxj],x2t[maxj]],[ttt[maxj]-self.tshift,ttt[maxj]-self.tshift],'--',c=c)

        plt.show()
        if self.interactivePlot:
            plt.pause(0.001)
        self.currentPlot = fig
        self.figureList.append(fig)
    def writeInputFile(self,xlo,xhi,nxstep,tlo,thi,ntstep):
        with open(self.inputFileName, 'w') as inputFile:
            inputFile.write('! Python-generated input file for Thermochimica\n')
            if float(nxstep) > 0:
                xstep = (float(xhi)-float(xlo))/float(nxstep)
            else:
                xstep = 0
            inputFile.write('x          = ' + str(xlo) + ':' + str(xhi) + ':' + str(xstep) + '\n')
            if float(ntstep) > 0:
                tstep = (float(thi)-float(tlo))/float(ntstep)
            else:
                tstep = 0
            inputFile.write('temperature          = ' + str(tlo) + ':' + str(thi) + ':' + str(tstep) + '\n')
            inputFile.write('pressure          = ' + str(self.pressure) + '\n')
            inputFile.write('temperature unit         = ' + self.tunit + '\n')
            inputFile.write('pressure unit          = ' + self.punit + '\n')
            inputFile.write('mass unit          = \'' + self.munit + '\'\n')
            inputFile.write('iEl         = ' + str(thermoTools.atomic_number_map.index(self.el1)+1) + ' ' + str(thermoTools.atomic_number_map.index(self.el2)+1) + '\n')
            inputFile.write('data file         = ' + self.datafile + '\n')
    def addLabel(self,xlab,tlab):
        self.writeInputFile(xlab,xlab,0,tlab,tlab,0)
        subprocess.run(['./bin/PhaseDiagramDataGen',self.inputFileName])
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
        self.labels.append([[xlab,tlab],'+'.join(labelName)])
        self.processPhaseDiagramData()
    def refineLimit(self,x,res):
        maxit = 10
        if x == 0:
            for i in range(len(self.x0data[1])-1):
                nit = 0
                while ((self.x0data[1][i+1] - self.x0data[2][i]) > res) and (nit < maxit):
                    nit += 1
                    self.writeInputFile(0,0.001,2,self.x0data[2][i]-self.tshift,self.x0data[1][i+1]-self.tshift,4)
                    self.runCalc()
        if x == 1:
            for i in range(len(self.x1data[1])-1):
                nit = 0
                while ((self.x1data[1][i+1] - self.x1data[2][i]) > res) and (nit < maxit):
                    nit += 1
                    self.writeInputFile(0.999,1,2,self.x1data[2][i]-self.tshift,self.x1data[1][i+1]-self.tshift,4)
                    self.runCalc()
    def autoRefine(self,res):
        nIt = 0
        while nIt < 4:
            nIt = nIt + 1
            maxArea = 0
            self.phaseBoundaries()

            phasePolyPoints = [[] for i in range(len(self.phases))]

            for j in range(len(self.x0data[1])):
                try:
                    i = self.phases.index(self.x0data[0][j])
                    phasePolyPoints[i].append([[0,self.x0data[1][j]]])
                    phasePolyPoints[i].append([[0,self.x0data[2][j]]])
                except:
                    continue
            for j in range(len(self.x1data[1])):
                try:
                    i = self.phases.index(self.x1data[0][j])
                    phasePolyPoints[i].append([[1,self.x1data[1][j]]])
                    phasePolyPoints[i].append([[1,self.x1data[2][j]]])
                except:
                    continue

            # plot 2-phase region boundaries
            for j in range(len(self.boundaries)):
                polygonPoints = []
                inds = [i for i, k in enumerate(self.b) if k == j]
                if len(inds) < 2:
                    continue
                ttt = self.ts[inds]
                x1t = self.x1[inds]
                x2t = self.x2[inds]
                for i in range(len(inds)):
                    polygonPoints.append([x1t[i],ttt[i]])
                for i in reversed(range(len(inds))):
                    polygonPoints.append([x2t[i],ttt[i]])
                phaseOutline = Polygon(polygonPoints).buffer(0)
                self.outline = self.outline.buffer(0) - phaseOutline
                minj = np.argmin(ttt)
                maxj = np.argmax(ttt)
                for i in range(len(self.phases)):
                    if self.boundaries[j][0] == self.phases[i]:
                        phasePolyPoints[i].append(polygonPoints[:len(inds)])
                    if self.boundaries[j][1] == self.phases[i]:
                        phasePolyPoints[i].append(list(reversed(polygonPoints))[:len(inds)])

            for i in range(len(self.phases)):
                if self.congruentFound[i]:
                    print(f'Warning: congruent phase transformation found, auto refine will skip {self.phases[i]}')
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
                        self.outline = self.outline - phaseOutline
                    except:
                        pass

            xs = []
            ys = []
            subres = int(np.ceil(np.sqrt(res)))
            try:
                oxlo, otlo, oxhi, othi = self.outline.bounds
            except:
                continue
            xindices = np.linspace(oxlo, oxhi, subres)
            yindices = np.linspace(otlo, othi, subres)
            horizontal_splitters = [LineString([(x, yindices[0]), (x, yindices[-1])]) for x in xindices]
            vertical_splitters = [LineString([(xindices[0], y), (xindices[-1], y)]) for y in yindices]
            # If the outline contains non-polygon shapes (like lines) it will be a GeometryCollection instead
            # and we need to remove those non-polygon shapes so it can be a MultiPolygon again
            if isinstance(self.outline,GeometryCollection):
                self.outline = MultiPolygon([shape for shape in list(self.outline.geoms) if isinstance(shape,Polygon)])
            for splitter in vertical_splitters:
                try:
                    self.outline = MultiPolygon(split(self.outline, splitter))
                except:
                    continue
            for splitter in horizontal_splitters:
                try:
                    self.outline = MultiPolygon(split(self.outline, splitter))
                except:
                    continue
            for tempOutline in list(self.outline.geoms):
                if (tempOutline.area / (self.maxt-self.mint)) < (1 / (10*res**2)):
                    continue
                maxArea = max(tempOutline.area / (self.maxt-self.mint),maxArea)
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
                    calc = [ys[i],self.pressure,1-xs[i],xs[i]]
                    calcList.append(calc)
                thermoTools.WriteRunCalculationList(self.inputFileName,self.datafile,[self.el1,self.el2],calcList,tunit=self.tunit,punit=self.punit,munit=self.munit,printMode=0)
                print('Thermochimica calculation initiated.')
                thermoTools.RunRunCalculationList(self.inputFileName)
                print('Thermochimica calculation finished.')
                self.processPhaseDiagramData()

            # Test the minimum subgrid region area to see if converged
            if maxArea < 1 / (10*res**2):
                break
            elif any(self.congruentFound):
                break
    def autoRefine2Phase(self,res):
        self.phaseBoundaries()
        # Expand two-phase regions
        tres = (self.maxt-self.mint)/res
        xs = []
        ys = []
        for j in range(len(self.boundaries)):
            inds = [i for i, k in enumerate(self.b) if k == j]
            if len(inds) < 2:
                continue
            ttt = self.ts[inds]
            x1t = self.x1[inds]
            x2t = self.x2[inds]
            tbound = max(self.mint,ttt[0]-tres*3)
            for k in np.arange(tbound,max(self.mint,ttt[0]-tres/3),tres/3):
                ys.append(k)
                xs.append((x1t[0] + x2t[0])/2)
                ys.append(k)
                xs.append((0.99*x1t[0] + 0.01*x2t[0]))
                ys.append(k)
                xs.append((0.01*x1t[0] + 0.99*x2t[0]))
            tbound = min(self.maxt,ttt[-1]+tres*3)
            for k in np.arange(min(self.maxt,ttt[-1]+tres/3),tbound,tres/3):
                ys.append(k)
                xs.append((x1t[-1] + x2t[-1])/2)
                ys.append(k)
                xs.append((0.99*x1t[-1] + 0.01*x2t[-1]))
                ys.append(k)
                xs.append((0.01*x1t[-1] + 0.99*x2t[-1]))

        if len(xs) > 0:
            calcList = []
            for i in range(len(xs)):
                calc = [ys[i],self.pressure,1-xs[i],xs[i]]
                calcList.append(calc)
            thermoTools.WriteRunCalculationList(self.inputFileName,self.datafile,[self.el1,self.el2],calcList,tunit=self.tunit,punit=self.punit,munit=self.munit,printMode=0)
            print('Thermochimica calculation initiated.')
            thermoTools.RunRunCalculationList(self.inputFileName)
            print('Thermochimica calculation finished.')
            self.processPhaseDiagramData()

        nIt = 0
        while nIt < 4:
            nIt = nIt + 1
            maxGap = 0
            self.phaseBoundaries()
            # Refine two-phase region density
            xs = []
            ys = []
            for j in range(len(self.boundaries)):
                inds = [i for i, k in enumerate(self.b) if k == j]
                if len(inds) < 2:
                    continue
                ttt = self.ts[inds]
                x1t = self.x1[inds]
                x2t = self.x2[inds]
                for i in range(len(ttt)-1):
                    gap = np.sqrt(((ttt[i]-ttt[i+1])/(self.maxt-self.mint))**2+(x1t[i]-x1t[i+1])**2+(x2t[i]-x2t[i+1])**2)
                    maxGap = max(gap,maxGap)
                    if gap > 1/res:
                        step = tres*((ttt[i+1] - ttt[i])/(self.maxt - self.mint))/gap
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
                    calc = [ys[i],self.pressure,1-xs[i],xs[i]]
                    calcList.append(calc)
                thermoTools.WriteRunCalculationList(self.inputFileName,self.datafile,[self.el1,self.el2],calcList,tunit=self.tunit,punit=self.punit,munit=self.munit,printMode=0)
                print('Thermochimica calculation initiated.')
                thermoTools.RunRunCalculationList(self.inputFileName)
                print('Thermochimica calculation finished.')
                self.processPhaseDiagramData()

            # Test the minimum difference between points to see if converged
            if maxGap <= 1/res:
                break
        self.gapLimit = 3*tres
    def autoLabel(self):
        self.phaseBoundaries()

        phasePolyPoints = [[] for i in range(len(self.phases))]

        for j in range(len(self.x0data[1])):
            try:
                i = self.phases.index(self.x0data[0][j])
                phasePolyPoints[i].append([[0,self.x0data[1][j]]])
                phasePolyPoints[i].append([[0,self.x0data[2][j]]])
            except:
                continue
        for j in range(len(self.x1data[1])):
            try:
                i = self.phases.index(self.x1data[0][j])
                phasePolyPoints[i].append([[1,self.x1data[1][j]]])
                phasePolyPoints[i].append([[1,self.x1data[2][j]]])
            except:
                continue

        # plot 2-phase region boundaries
        for j in range(len(self.boundaries)):
            polygonPoints = []
            inds = [i for i, k in enumerate(self.b) if k == j]
            if len(inds) < 2:
                continue
            ttt = self.ts[inds]
            x1t = self.x1[inds]
            x2t = self.x2[inds]
            for i in range(len(inds)):
                polygonPoints.append([x1t[i],ttt[i]])
            for i in reversed(range(len(inds))):
                polygonPoints.append([x2t[i],ttt[i]])
            phaseOutline = Polygon(polygonPoints)#.buffer(0)
            center = list(phaseOutline.centroid.coords)[0]
            if self.label2phase:
                self.labels.append([[center[0],center[1]-self.tshift],'+'.join(self.boundaries[j])])
            for i in range(len(self.phases)):
                if self.boundaries[j][0] == self.phases[i]:
                    phasePolyPoints[i].append(polygonPoints[:len(inds)])
                if self.boundaries[j][1] == self.phases[i]:
                    phasePolyPoints[i].append(list(reversed(polygonPoints))[:len(inds)])

        if self.label1phase:
            for i in range(len(self.phases)):
                if self.congruentFound[i]:
                    print(f'Warning: congruent phase transformation found, auto label will skip {self.phases[i]}')
                    continue
                segcenters = []
                if len(phasePolyPoints[i]) < 2:
                    continue
                for j in range(len(phasePolyPoints[i])):
                    segcenters.append(tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), phasePolyPoints[i][j]), [len(phasePolyPoints[i][j])] * 2)))
                center = tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), segcenters), [len(segcenters)] * 2))
                self.labels.append([[center[0],center[1]-self.tshift],self.phases[i]])
    def makeBackup(self):
        self.backup = diagram(self.datafile, False, self.interactivePlot)
        self.backup.mint = self.mint
        self.backup.maxt = self.maxt
        self.backup.ts = copy.deepcopy(self.ts)
        self.backup.x1 = copy.deepcopy(self.x1)
        self.backup.x2 = copy.deepcopy(self.x2)
        self.backup.p1 = copy.deepcopy(self.p1)
        self.backup.p2 = copy.deepcopy(self.p2)
        self.backup.x0data = copy.deepcopy(self.x0data)
        self.backup.x1data = copy.deepcopy(self.x1data)
        self.backup.labels = copy.deepcopy(self.labels)
        self.backup.outline = copy.deepcopy(self.outline)
        self.backup.pressure = self.pressure
        self.backup.inputFileName = self.inputFileName
        self.backup.outputFileName = self.outputFileName
        self.backup.plotMarker = self.plotMarker
        self.backup.plotColor = self.plotColor
        self.backup.el1 = self.el1
        self.backup.el2 = self.el2
        self.backup.tunit = self.tunit
        self.backup.punit = self.punit
        self.backup.munit = self.munit
        self.backup.tshift = self.tshift
        self.backup.exportFormat = self.exportFormat
        self.backup.exportFileName = self.exportFileName
        self.backup.exportDPI = self.exportDPI
        self.backup.resRef = self.resRef
        self.backup.resSmooth = self.resSmooth
        self.backup.gapLimit = self.gapLimit
        self.backup.label1phase = self.label1phase
        self.backup.label2phase = self.label2phase
        self.backup.experimentalData = self.experimentalData
        self.backup.experimentNames = self.experimentNames
        self.backup.experimentColor = self.experimentColor
        self.backup.pointDetails = self.pointDetails
        self.backup.pointIndex = self.pointIndex
        self.backup.suppressed = self.suppressed
        self.backup.loadedDiagram = self.loadedDiagram
        self.backup.loaded = self.loaded
        self.backup.saveDataName = self.saveDataName
    def exportPlot(self):
        try:
            self.currentPlot.savefig(f'{self.exportFileName}.{self.exportFormat}', format=self.exportFormat, dpi=self.exportDPI)
            return 0
        except:
            return 1
    def addData(self,datafile,expName):
        newData = []
        with open(datafile) as f:
            data = csv.reader(f)
            # next(data, None)  # skip the header WHY SHOULD THERE BE A HEADER
            for row in data:
                newrow = []
                for number in row:
                    newrow.append(float(number))
                newData.append(newrow)
        self.experimentalData.append(np.array(newData))
        self.experimentNames.append(expName)
