import json
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import copy
from itertools import cycle
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from functools import reduce
import operator
import csv
import thermoTools
from phaseDiagramCommon import *

phaseIncludeTol = 1e-8

class diagram:
    def __init__(self, datafile, active, interactivePlot):
        self.datafile = datafile
        self.active = active
        self.interactivePlot = interactivePlot
        self.mint = 1e5
        self.maxt = 0
        self.pdPoints = []
        self.el1 = ''
        self.el2 = ''
        self.elementsUsed = []
        self.massLabels = ['','']
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
        self.compoundScale = False
        self.gapLimit = np.Inf
        self.figureList = []
        self.boundaries = []
        self.phases = []
        self.b = []
        self.congruentFound = [False for _ in range(len(self.phases))]
        self.label1phase = True
        self.label2phase = True
        self.experimentalData = []
        self.experimentNames = []
        self.experimentColor = 'bland'
        self.showExperiment = True
        self.loadedDiagram = []
        self.loaded = False
        self.showLoaded = True
        self.saveDataName = 'savedDiagram'
        self.fuzzy = False
    def run(self,ntstep,nxstep,pressure,tunit,punit,xlo,xhi,tlo,thi,el1,el2,munit,fuzzy=False):
        self.pressure = pressure
        self.tunit = tunit
        self.punit = punit
        self.munit = munit
        self.el1 = el1
        self.el2 = el2
        self.elementsUsed = [el1,el2]
        self.massLabels   = [el1,el2]
        self.writeInputFile(xlo,xhi,nxstep,tlo,thi,ntstep)
        self.pdPoints = []
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
        # Get fuzzy stoichiometry setting
        self.fuzzy = fuzzy
        self.labels = []
        self.resRef = 7
        self.resSmooth = 7
        self.gapLimit = (self.maxt - self.mint) / 2
        self.experimentalData = []
        self.experimentNames = []
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
        autoRefine(self,self.resRef**2,np.array([[1,0],[0,1]]))
        self.resRef += 1
    def autoSmooth(self):
        autoRefine2Phase(self,self.resSmooth**2,np.array([[1,0],[0,1]]))
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
        elem = [self.el1,self.el2]
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
                t = data[i]['temperature']
                boundPhases = []
                boundComps = []
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > phaseIncludeTol):
                            boundPhases.append(phaseName)
                            boundComps.append(data[i][phaseType][phaseName]['elements'][self.el2]['mole fraction of phase by element'])
                x = [boundComps[0],boundComps[1]]
                p = [boundPhases[0],boundPhases[1]]
                conc = [data[i]["elements"][self.el1]["moles"],data[i]["elements"][self.el2]["moles"]]
                en = data[i]["integral Gibbs energy"]
                it = data[i]["GEM iterations"]
                self.pdPoints.append(pdPoint(elem,t,conc,p,x,en,it))
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
        self.pdPoints.sort(key=lambda x: x.t)

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
    def makePlot(self):
        phaseBoundaries(self)
        # Start figure
        fig = plt.figure()
        plt.ioff()
        if self.interactivePlot:
            plt.ion()
        ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])

        self.plotLines(ax,diagram=self,linestyle=self.plotMarker,plotColor=self.plotColor)

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
            self.plotLines(ax,diagram=self.loadedDiagram,linestyle='--',plotColor='black')

        plt.show()
        if self.interactivePlot:
            plt.pause(0.001)
        self.currentPlot = fig
        self.figureList.append(fig)
    def plotLines(self,ax,diagram,linestyle,plotColor):
        bEdgeLine = [[False,False] for _ in range(len(diagram.boundaries))]
        # Plot along x=0 and x=1 self.boundaries (this is the worst code I've ever written)
        for j in range(len(diagram.x0data[1])):
            if not diagram.x0data[0][j] in diagram.phases:
                continue
            i = diagram.phases.index(diagram.x0data[0][j])
            if j > 0:
                # ax.plot(0,diagram.x0data[1][j],'kv')
                match = []
                for k in range(len(diagram.boundaries)):
                    if (diagram.x0data[0][j] in diagram.boundaries[k]) and (diagram.x0data[0][j-1] in diagram.boundaries[k]):
                        inds = [i for i, l in enumerate(diagram.b) if l == k]
                        if len(inds) < 2:
                            continue
                        bind = diagram.boundaries[k].index(diagram.x0data[0][j])
                        if bind == 0:
                            t = np.array([diagram.pdPoints[i].t for i in inds])
                            x = np.array([diagram.pdPoints[i].phaseConcentrations[0] for i in inds])
                            minj = np.argmin(x)
                            length = (0 - x[minj])**2 + (diagram.x0data[1][j] - t[minj])**2
                            match.append([length,k,x[minj],t[minj]])
                        elif bind == 1:
                            t = np.array([diagram.pdPoints[i].t for i in inds])
                            x = np.array([diagram.pdPoints[i].phaseConcentrations[1] for i in inds])
                            minj = np.argmin(x)
                            length = (0 - x[minj])**2 + (diagram.x0data[1][j] - t[minj])**2
                            match.append([length,k,x[minj],t[minj]])
                if len(match) > 0:
                    match = np.array(match)
                    matchind = np.argmin(match[:,0])
                    k = int(match[matchind,1])
                    inds = [i for i, l in enumerate(diagram.b) if l == k]
                    t = np.array([diagram.pdPoints[i].t for i in inds])
                    ax.plot([0,match[matchind,2]],[diagram.x0data[1][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                    if match[matchind,3] == np.min(t):
                        bEdgeLine[k][0] = True
                    if match[matchind,3] == np.max(t):
                        bEdgeLine[k][1] = True
            if j < len(diagram.x0data[1]) - 1:
                # ax.plot(0,diagram.x0data[2][j],'k^')
                match = []
                for k in range(len(diagram.boundaries)):
                    if (diagram.x0data[0][j] in diagram.boundaries[k]) and (diagram.x0data[0][j+1] in diagram.boundaries[k]):
                        inds = [i for i, l in enumerate(diagram.b) if l == k]
                        if len(inds) < 2:
                            continue
                        bind = diagram.boundaries[k].index(diagram.x0data[0][j])
                        if bind == 0:
                            t = np.array([diagram.pdPoints[i].t for i in inds])
                            x = np.array([diagram.pdPoints[i].phaseConcentrations[0] for i in inds])
                            minj = np.argmin(x)
                            length = (0 - x[minj])**2 + (diagram.x0data[2][j] - t[minj])**2
                            match.append([length,k,x[minj],t[minj]])
                        elif bind == 1:
                            t = np.array([diagram.pdPoints[i].t for i in inds])
                            x = np.array([diagram.pdPoints[i].phaseConcentrations[1] for i in inds])
                            minj = np.argmin(x)
                            length = (0 - x[minj])**2 + (diagram.x0data[2][j] - t[minj])**2
                            match.append([length,k,x[minj],t[minj]])
                if len(match) > 0:
                    match = np.array(match)
                    matchind = np.argmin(match[:,0])
                    k = int(match[matchind,1])
                    inds = [i for i, l in enumerate(diagram.b) if l == k]
                    t = np.array([diagram.pdPoints[i].t for i in inds])
                    ax.plot([0,match[matchind,2]],[diagram.x0data[2][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                    if match[matchind,3] == np.min(t):
                        bEdgeLine[k][0] = True
                    if match[matchind,3] == np.max(t):
                        bEdgeLine[k][1] = True
        for j in range(len(diagram.x1data[1])):
            if not diagram.x1data[0][j] in diagram.phases:
                continue
            i = diagram.phases.index(diagram.x1data[0][j])
            if j > 0:
                # ax.plot(1,diagram.x1data[1][j],'kv')
                match = []
                for k in range(len(diagram.boundaries)):
                    if (diagram.x1data[0][j] in diagram.boundaries[k]) and (diagram.x1data[0][j-1] in diagram.boundaries[k]):
                        inds = [i for i, l in enumerate(diagram.b) if l == k]
                        if len(inds) < 2:
                            continue
                        bind = diagram.boundaries[k].index(diagram.x1data[0][j])
                        if bind == 0:
                            t = np.array([diagram.pdPoints[i].t for i in inds])
                            x = np.array([diagram.pdPoints[i].phaseConcentrations[0] for i in inds])
                            maxj = np.argmax(x)
                            length = (1 - x[maxj])**2 + (diagram.x1data[1][j] - t[maxj])**2
                            match.append([length,k,x[maxj],t[maxj]])
                        elif bind == 1:
                            t = np.array([diagram.pdPoints[i].t for i in inds])
                            x = np.array([diagram.pdPoints[i].phaseConcentrations[1] for i in inds])
                            maxj = np.argmax(x)
                            length = (1 - x[maxj])**2 + (diagram.x1data[1][j] - t[maxj])**2
                            match.append([length,k,x[maxj],t[maxj]])
                if len(match) > 0:
                    match = np.array(match)
                    matchind = np.argmin(match[:,0])
                    k = int(match[matchind,1])
                    inds = [i for i, l in enumerate(diagram.b) if l == k]
                    t = np.array([diagram.pdPoints[i].t for i in inds])
                    ax.plot([1,match[matchind,2]],[diagram.x1data[1][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                    if match[matchind,3] == np.min(t):
                        bEdgeLine[k][0] = True
                    if match[matchind,3] == np.max(t):
                        bEdgeLine[k][1] = True
            if j < len(diagram.x1data[1]) - 1:
                # ax.plot(1,diagram.x1data[2][j],'k^')
                match = []
                for k in range(len(diagram.boundaries)):
                    if (diagram.x1data[0][j] in diagram.boundaries[k]) and (diagram.x1data[0][j+1] in diagram.boundaries[k]):
                        inds = [i for i, l in enumerate(diagram.b) if l == k]
                        if len(inds) < 2:
                            continue
                        bind = diagram.boundaries[k].index(diagram.x1data[0][j])
                        if bind == 0:
                            t = np.array([diagram.pdPoints[i].t for i in inds])
                            x = np.array([diagram.pdPoints[i].phaseConcentrations[0] for i in inds])
                            maxj = np.argmax(x)
                            length = (1 - x[maxj])**2 + (diagram.x1data[2][j] - x[maxj])**2
                            match.append([length,k,x[maxj],t[maxj]])
                        elif bind == 1:
                            t = np.array([diagram.pdPoints[i].t for i in inds])
                            x = np.array([diagram.pdPoints[i].phaseConcentrations[1] for i in inds])
                            maxj = np.argmax(x)
                            length = (1 - x[maxj])**2 + (diagram.x1data[2][j] - t[maxj])**2
                            match.append([length,k,x[maxj],t[maxj]])
                if len(match) > 0:
                    match = np.array(match)
                    matchind = np.argmin(match[:,0])
                    k = int(match[matchind,1])
                    inds = [i for i, l in enumerate(diagram.b) if l == k]
                    t = np.array([diagram.pdPoints[i].t for i in inds])
                    ax.plot([1,match[matchind,2]],[diagram.x1data[2][j]-self.tshift,match[matchind,3]-self.tshift],'k-')
                    if match[matchind,3] == np.min(t):
                        bEdgeLine[k][0] = True
                    if match[matchind,3] == np.max(t):
                        bEdgeLine[k][1] = True

        # plot 2-phase region boundaries
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(diagram.boundaries))))
        for j in range(len(diagram.boundaries)):
            if plotColor == 'colorful':
                c = next(color)
            else:
                c = 'k'
            inds = [i for i, k in enumerate(diagram.b) if k == j]
            if len(inds) < 2:
                continue
            ttt = np.array([diagram.pdPoints[i].t for i in inds])
            x1t = np.array([diagram.pdPoints[i].phaseConcentrations[0] for i in inds])
            x2t = np.array([diagram.pdPoints[i].phaseConcentrations[1] for i in inds])
            ax.plot(x1t,ttt-self.tshift,linestyle,c=c)
            ax.plot(x2t[::-1],ttt[::-1]-self.tshift,linestyle,c=c)
            minj = np.argmin(ttt)
            maxj = np.argmax(ttt)
            # plot invariant temperatures
            if (ttt[minj] > diagram.mint) and not(bEdgeLine[j][0]):
                ax.plot([x1t[minj],x2t[minj]],[ttt[minj]-self.tshift,ttt[minj]-self.tshift],linestyle,c=c)
            if (ttt[maxj] < diagram.maxt) and not(bEdgeLine[j][1]):
                ax.plot([x1t[maxj],x2t[maxj]],[ttt[maxj]-self.tshift,ttt[maxj]-self.tshift],linestyle,c=c)
    def writeInputFile(self,xlo,xhi,nxstep,tlo,thi,ntstep):
        with open(self.inputFileName, 'w') as inputFile:
            inputFile.write('! Python-generated input file for Thermochimica\n')
            if float(nxstep) > 0:
                xstep = (float(xhi)-float(xlo))/float(nxstep)
            else:
                xstep = 0
            inputFile.write(f'x                 = {xlo}:{xhi}:{xstep}\n')
            if float(ntstep) > 0:
                tstep = (float(thi)-float(tlo))/float(ntstep)
            else:
                tstep = 0
            inputFile.write(f'temperature       = {tlo}:{thi}:{tstep}\n')
            inputFile.write(f'pressure          = {self.pressure}\n')
            inputFile.write(f'temperature unit  = \'{self.tunit}\'\n')
            inputFile.write(f'pressure unit     = \'{self.punit}\'\n')
            inputFile.write(f'mass unit         = \'{self.munit}\'\n')
            inputFile.write(f'iEl               = {thermoTools.atomic_number_map.index(self.el1)+1} {thermoTools.atomic_number_map.index(self.el2)+1}\n')
            inputFile.write(f'data file         = {self.datafile}\n')
            # Fuzzy stoichiometry settings
            inputFile.write(f'fuzzy             = {".TRUE." if self.fuzzy else ".FALSE."}\n')
            inputFile.write(f'gibbs min         = {".TRUE." if self.fuzzy else ".FALSE."}\n')
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
    def autoLabel(self):
        autoLabel(self)
    def makeBackup(self):
        self.backup = diagram(self.datafile, False, self.interactivePlot)
        self.backup.mint = self.mint
        self.backup.maxt = self.maxt
        self.backup.pdPoints = copy.deepcopy(self.pdPoints)
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
        self.backup.elementsUsed = copy.deepcopy(self.elementsUsed)
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
        self.backup.loadedDiagram = self.loadedDiagram
        self.backup.loaded = self.loaded
        self.backup.saveDataName = self.saveDataName
    def exportPlot(self):
        # Make sure there is an open plot to save
        if not plt.fignum_exists(self.currentPlot.number):
            self.makePlot()
        try:
            self.currentPlot.savefig(f'outputs/{self.exportFileName}.{self.exportFormat}', format=self.exportFormat, dpi=self.exportDPI)
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
