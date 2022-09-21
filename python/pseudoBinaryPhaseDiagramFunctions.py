import json
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import copy
import scipy.optimize
from itertools import cycle
import csv
import thermoTools

# For boundaries of phase regions where both sides have (# phases) < (# elements), only plot points within phaseFractionTol of the boundary
phaseFractionTol = 1e-2
# Below this tolerance, set phase fraction = 0
phaseIncludeTol = 1e-8

class diagram:
    def __init__(self, datafile, active, interactivePlot):
        self.datafile = datafile
        self.active = active
        self.interactivePlot = interactivePlot
        self.children = []
        self.inputFileName = 'inputs/pythonPhaseDiagramInput.ti'
        self.outputFileName = 'thermoout.json'
        self.plotMarker = '-'
        self.plotColor = 'colorful'
        if self.active:
            self.backup = diagram(self.datafile, False, self.interactivePlot)
        else:
            self.backup = []
        self.currentPlot = []
        self.exportFormat = 'png'
        self.exportFileName = 'thermochimicaPseudoBinaryPhaseDiagram'
        self.exportDPI = 100
        self.mint = 1e5
        self.maxt = 0
        self.pressure = 1
        self.points = []
        self.labels = []
        self.elementsUsed = []
        self.nElementsUsed = 0
        self.massLabels = ['','']
        self.sum1 = 0
        self.sum2 = 0
        self.plane = np.array([0,0])
        self.tunit = 'K'
        self.punit = 'atm'
        self.munit = 'moles'
        self.tshift = 0
        # I don't think anyone is going to change this scale, so consider this a debug setting
        self.normalizeX = False
        self.figureList = []
        self.experimentalData = []
        self.experimentNames = []
        self.experimentColor = 'bland'
        self.showExperiment = True
        self.showExperimentLegend = True
    def initRun(self,pressure,tunit,punit,plane,sum1,sum2,mint,maxt,elementsUsed,massLabels,munit,tshift):
        self.mint = mint
        self.maxt = maxt
        self.pressure = pressure
        self.points = []
        self.labels = []
        self.elementsUsed = elementsUsed
        self.nElementsUsed = len(elementsUsed)
        self.massLabels = massLabels
        self.sum1 = sum1
        self.sum2 = sum2
        self.plane = np.array(plane)
        self.tunit = tunit
        self.punit = punit
        self.munit = munit
        self.tshift = tshift
    def runCalc(self,xlo,xhi,nxstep,tlo,thi,ntstep):
        xs = np.array([np.linspace((1-xlo)*self.plane[0,i] + xlo*self.plane[1,i],(1-xhi)*self.plane[0,i] + xhi*self.plane[1,i],nxstep) for i in range(self.nElementsUsed)]).T
        temps = np.linspace(tlo,thi,ntstep)
        calcList = []
        for t in temps:
            ioff = 0
            for x in xs:
                toff = 0
                if (t > tlo) and (t < thi - (thi-tlo)/ntstep):
                    toff = ioff * ((thi-tlo)/ntstep)/10
                    ioff += 1
                    ioff = ioff % 10
                calc = [t+toff,self.pressure]
                calc.extend([x[i] for i in range(self.nElementsUsed)])
                calcList.append(calc)
        thermoTools.WriteRunCalculationList(self.inputFileName,self.datafile,self.elementsUsed,calcList,tunit=self.tunit,punit=self.punit,munit=self.munit,printMode=0)
        print('Thermochimica calculation initiated.')
        subprocess.run(['./bin/RunCalculationList',self.inputFileName])
        print('Thermochimica calculation finished.')
    def processPhaseDiagramData(self):
        f = open(self.outputFileName,)
        data = json.load(f)
        f.close()
        if list(data.keys())[0] != '1':
            print('Output does not contain data series')
            return
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
                    self.points.append([[data[i]['temperature'],boundComps[0],boundPhases],[data[i]['temperature'],boundComps[1],boundPhases]])
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
                temppoints = []
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
                        # if intSum == 1:
                        #     # If we are on the far boundary, the first phase is not included
                        #     omitPhase.remove(omitPhase[0])
                        # for k in range(self.nElementsUsed - 2):
                        #     # Check all the other boundaries
                        #     if intersect[k+1] == 0:
                        #         # If none of this is used, it is not included
                        #         omitPhase.remove(omitPhase[k+1])
                        # self.points.append([data[i]['temperature'],intersect[0],omitPhase])
                        temppoints.append([data[i]['temperature'],intersect[0],allPhases])
                self.points.append(temppoints)
            elif nPhases > 1 and False:
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
                if not (len(point[0][2]) == len(boundaries[j])):
                    continue
                for phase in point[0][2]:
                    if not (phase in boundaries[j]):
                        thisMatch = False
                        break
                if thisMatch:
                    b.append(j)
                    repeat = True
            if not(repeat):
                b.append(len(boundaries))
                boundaries.append(point[0][2])

        for j in range(len(boundaries)):
            inds = [i for i, k in enumerate(b) if k == j]
            if len(inds) < 2:
                continue
            temppoints1 =[self.points[i][0][1] for i in inds]
            temppoints2 =[self.points[i][1][1] for i in inds]
            if temppoints1[0] > temppoints2[0]:
                dir = True
            else:
                dir = False
            extraBound = []
            for i in range(len(temppoints1)):
                if (temppoints1[i] > temppoints2[i]) != dir:
                    extraBound.append(i)
            if len(extraBound):
                boundaries.append(boundaries[j])
                for k in extraBound:
                    b[inds[k]] = len(boundaries)-1
        
        # Start figure
        fig = plt.figure()
        plt.ion()
        ax = fig.add_axes([0.12, 0.1, 0.85, 0.85])

        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(boundaries))))
        for j in range(len(boundaries)):
            if self.plotColor == 'colorful':
                c = next(color)
            else:
                c = 'k'
            inds = [i for i, k in enumerate(b) if k == j]
            if len(inds) < 2:
                continue
            plotPoints = np.empty([0,2])
            temppoints = np.array([[self.points[i][0][1],self.points[i][0][0]] for i in inds])
            plotPoints = np.append(plotPoints,temppoints[temppoints[:,1].argsort()], axis=0)
            temppoints = np.array([[self.points[i][1][1],self.points[i][1][0]] for i in inds])
            plotPoints = np.append(plotPoints,temppoints[temppoints[:,1].argsort()][::-1], axis=0)
            if self.normalizeX:
                plotX = plotPoints[:,0]
            else:
                plotX = self.unscaleX(plotPoints[:,0])
            ax.plot(plotX,plotPoints[:,1]-self.tshift,self.plotMarker,c=c, label='_nolegend_')

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
        title = " $-$ ".join(self.massLabels)
        ax.set_title(r'{0} phase diagram'.format(title))
        ax.set_xlabel(r'Mole fraction {0}'.format(self.massLabels[1]))
        tunit_display = r"$^\circ$" if self.tunit == "C" else ""
        ax.set_ylabel(f'Temperature [{tunit_display}{self.tunit}]')
        for lab in self.labels:
            plt.text(float(lab[0][0]),float(lab[0][1]),lab[1], ha="center")
        if self.showExperimentLegend and len(self.experimentalData):
            ax.legend(loc = 0)
        plt.show()
        plt.pause(0.001)
        self.currentPlot = fig
        self.figureList.append(fig)
    def addLabel(self,xlab,tlab):
        # label x-coords are going to come in scaled to axis
        xlabrun = xlab
        if not self.normalizeX:
            if xlab > 0:
                xlabrun = 1/(1+((1-xlab)/xlab)*(self.sum1/self.sum2))

        self.runCalc(xlabrun,xlabrun,1,tlab,tlab,1)
        f = open(self.outputFileName,)
        data = json.load(f)
        f.close()
        if len(data.keys()) == 0:
            print('Output does not contain data series')
            return
        if list(data.keys())[0] != '1':
            print('Output does not contain data series')
            return
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
    def makeBackup(self):
        self.backup = diagram(self.datafile, False, self.interactivePlot)
        self.backup.mint = self.mint
        self.backup.maxt = self.maxt
        self.backup.labels = copy.deepcopy(self.labels)
        self.backup.pressure = self.pressure
        self.backup.inputFileName = self.inputFileName
        self.backup.outputFileName = self.outputFileName
        self.backup.plotMarker = self.plotMarker
        self.backup.plotColor = self.plotColor
        self.backup.tunit = self.tunit
        self.backup.punit = self.punit
        self.backup.munit = self.munit
        self.backup.tshift = self.tshift
        self.backup.exportFormat = self.exportFormat
        self.backup.exportFileName = self.exportFileName
        self.backup.exportDPI = self.exportDPI
        self.backup.points = copy.deepcopy(self.points)
        self.backup.elementsUsed = copy.deepcopy(self.elementsUsed)
        self.backup.nElementsUsed = self.nElementsUsed
        self.backup.massLabels = copy.deepcopy(self.massLabels)
        self.backup.sum1 = self.sum1
        self.backup.sum2 = self.sum2
        self.backup.plane = copy.deepcopy(self.plane)
        self.backup.normalizeX = self.normalizeX 
        self.backup.experimentalData = copy.deepcopy(self.experimentalData)
        self.backup.experimentNames = copy.deepcopy(self.experimentNames)
        self.backup.experimentColor = self.experimentColor 
        self.backup.showExperiment = self.showExperiment 
    def exportPlot(self):
        try:
            self.currentPlot.savefig(f'{self.exportFileName}.{self.exportFormat}', format=self.exportFormat, dpi=self.exportDPI)
            return 0
        except:
            return 1
    def unscaleX(self,scaledX):
        unscaledX = copy.deepcopy(scaledX)
        for i in range(len(scaledX)):
            if scaledX[i] > 0:
                unscaledX[i] = 1/(1+(self.sum2/self.sum1)*(1-scaledX[i])/scaledX[i])
        return unscaledX
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
        self.makePlot()