import PySimpleGUI as sg
import json
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import subprocess
import scipy.optimize

fname = 'pdmoru700-1000example.json'
elements = ['Mo', 'Pd', 'Ru']
nElements = len(elements)

mint = 1e5
maxt = 0

# For boundaries of phase regions where both sides have (# phases) < (# elements), only plot points within phaseFractionTol of the boundary
phaseFractionTol = 1e-2

plane = np.array([[0,0.5,0.5],[0.5,0,0.5]])
# plane = np.array([[0.05,0.00,0.95],[0.05,0.95,0.00]])

def line_intersection(line1, lines):
    l1 = np.array(line1)
    ls = np.array(lines)
    def diff(mu):
        sum = l1[0] + mu[0]*(l1[1] - l1[0])
        for i in range(nElements - 2):
            sum -= ls[i][0] + mu[i+1]*(ls[i][1] - ls[i][0])
        return sum
    # try:
    result = scipy.optimize.least_squares(diff, [0.5 for i in range(nElements-1)]).x
    # except:
    #     result = [-1000 for i in range(nElements-1)]
    #     print()
    return result

points = []

f = open(fname,)
data = json.load(f)
f.close()
if list(data.keys())[0] != '1':
    print('Output does not contain data series')
    exit()
for i in list(data.keys()):
    mint = min(mint,data[i]['temperature'])
    maxt = max(maxt,data[i]['temperature'])
    if (data[i]['# solution phases'] + data[i]['# pure condensed phases']) == nElements:
        allPhases = []
        phaseComps = []
        for phaseType in ['solution phases','pure condensed phases']:
            for phaseName in list(data[i][phaseType].keys()):
                if (data[i][phaseType][phaseName]['moles'] > 0):
                    allPhases.append(phaseName)
                    tempComp = []
                    for element in elements:
                        tempComp.append(data[i][phaseType][phaseName]['elements'][element]['mole fraction of phase by element'])
                    phaseComps.append(tempComp)
        # Loop over possible phase zone intersections with plane of interest
        for j in range(nElements):
            # Make list of phases on a (nElements-1) dimensional face through omission
            omitComps = phaseComps.copy()
            omitComps.remove(phaseComps[j])
            omitPhase = allPhases.copy()
            omitPhase.remove(allPhases[j])
            lineComps = []
            for k in range(nElements - 2):
                lineComps.append([omitComps[0],omitComps[k+1]])
            # Calculate intersection of this face with our plane of interest
            intersect = line_intersection(plane,lineComps)
            # Check that the intersection is within the valid bounds
            intSum = sum(intersect[1:])
            intTest = intSum <= 1
            for test in range(nElements - 1):
                intTest = intTest and (0 <= intersect[test]) and (intersect[test] <= 1)
            if intTest:
                if intSum == 1:
                    # If we are on the far boundary, the first phase is not included
                    omitPhase.remove(omitPhase[0])
                for k in range(nElements - 2):
                    # Check all the other boundaries
                    if intersect[k+1] == 0:
                        # If none of this is used, it is not included
                        omitPhase.remove(omitPhase[k+1])
                points.append([data[i]['temperature'],intersect[0],omitPhase])
    elif (data[i]['# solution phases'] + data[i]['# pure condensed phases']) > 0:
        boundPhases = []
        skipPoint = False
        phaseMoleSum = 0
        for phaseType in ['solution phases','pure condensed phases']:
            for phaseName in list(data[i][phaseType].keys()):
                phaseMoleSum += data[i][phaseType][phaseName]['moles']
        for phaseType in ['solution phases','pure condensed phases']:
            for phaseName in list(data[i][phaseType].keys()):
                if data[i][phaseType][phaseName]['moles'] > 0:
                    boundPhases.append(phaseName)
                    if phaseFractionTol < data[i][phaseType][phaseName]['moles']/phaseMoleSum < (1-phaseFractionTol):
                        skipPoint = True
                        break
            if skipPoint:
                break
        if skipPoint:
            boundPhases = []
            continue
        tempComp = np.zeros(nElements)
        for e in range(len(elements)):
            if elements[e] in data[i]['elements'].keys():
                tempComp[e] = data[i]['elements'][elements[e]]['moles']
        boundComps = np.linalg.norm(tempComp-plane[0])/np.linalg.norm(plane[1]-plane[0])
        points.append([data[i]['temperature'],boundComps,boundPhases])

boundaries = []
b = []
for point in points:
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
# plt.ion()
ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])

for j in range(len(boundaries)):
    inds = [i for i, k in enumerate(b) if k == j]
    if len(inds) < 2:
        continue
    plotPoints = np.array([[points[i][1],points[i][0]] for i, k in enumerate(b) if k == j])
    print()
    ax.plot(plotPoints[:,0],plotPoints[:,1],'.')

ax.set_xlim(0,1)
# ax.set_ylim(mint,maxt)
# ax.set_title(str(el1) + ' + ' + str(el2) + ' binary phase diagram')
# ax.set_xlabel('Mole fraction ' + str(el2))
ax.set_ylabel('Temperature [K]')
# for lab in labels:
#     plt.text(float(lab[0][0]),float(lab[0][1]),lab[1], ha="center")
plt.show()
plt.pause(0.001)
