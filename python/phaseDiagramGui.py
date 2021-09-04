import PySimpleGUI as sg
import json
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import subprocess
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
from shapely.prepared import prep
from descartes import PolygonPatch
from functools import reduce
import operator

def processPhaseDiagramData(fname, elx, ts, x1, x2, p1, p2, mint, maxt, x0data, x1data):
    f = open(fname,)
    data = json.load(f)
    f.close()
    if list(data.keys())[0] != '1':
        print('Output does not contain data series')
        exit()
    for i in list(data.keys()):
        mint = min(mint,data[i]['temperature'])
        maxt = max(maxt,data[i]['temperature'])
        if (data[i]['# solution phases'] + data[i]['# pure condensed phases']) == 2:
            ts.append(data[i]['temperature'])
            boundPhases = []
            boundComps = []
            for phaseType in ['solution phases','pure condensed phases']:
                for phaseName in list(data[i][phaseType].keys()):
                    if (data[i][phaseType][phaseName]['moles'] > 0):
                        boundPhases.append(phaseName)
                        boundComps.append(data[i][phaseType][phaseName]['elements'][elx]['mole fraction of phase by element'])
            x1.append(boundComps[0])
            x2.append(boundComps[1])
            p1.append(boundPhases[0])
            p2.append(boundPhases[1])
        elif (data[i]['# solution phases'] + data[i]['# pure condensed phases']) == 1:
            if not(elx in list(data[i]['elements'].keys())):
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > 0):
                            pname = phaseName
                if not(pname in x0data[0]):
                    x0data[0].append(pname)
                    x0data[1].append(data[i]['temperature'])
                    x0data[2].append(data[i]['temperature'])
                pindex = x0data[0].index(pname)
                x0data[1][pindex] = min(x0data[1][pindex],data[i]['temperature'])
                x0data[2][pindex] = max(x0data[2][pindex],data[i]['temperature'])
            elif float(data[i]['elements'][elx]['moles']) == 1:
                for phaseType in ['solution phases','pure condensed phases']:
                    for phaseName in list(data[i][phaseType].keys()):
                        if (data[i][phaseType][phaseName]['moles'] > 0):
                            pname = phaseName
                if not(pname in x1data[0]):
                    x1data[0].append(pname)
                    x1data[1].append(data[i]['temperature'])
                    x1data[2].append(data[i]['temperature'])
                pindex = x1data[0].index(pname)
                x1data[1][pindex] = min(x1data[1][pindex],data[i]['temperature'])
                x1data[2][pindex] = max(x1data[2][pindex],data[i]['temperature'])
    if len(x0data[1]) > 1:
        x0sort = [i[0] for i in sorted(enumerate(x0data[1]), key=lambda x:x[1])]
        phaseOrder = []
        for i in x0sort:
            phaseOrder.append(x0data[0][i])
        xtemp = [[],[],[]]
        xtemp[0] = phaseOrder
        xtemp[1] = sorted(x0data[1])
        xtemp[2] = sorted(x0data[2])
        x0data = xtemp
    if len(x1data[1]) > 1:
        x1sort = [i[0] for i in sorted(enumerate(x1data[1]), key=lambda x:x[1])]
        phaseOrder = []
        for i in x1sort:
            phaseOrder.append(x1data[0][i])
        xtemp = [[],[],[]]
        xtemp[0] = phaseOrder
        xtemp[1] = sorted(x1data[1])
        xtemp[2] = sorted(x1data[2])
        x1data = xtemp
    return mint, maxt

def runCalc(el1, el2, ts, x1, x2, p1, p2, mint, maxt, labels, x0data, x1data):
    print('Thermochimica calculation initiated.')
    subprocess.run(['./bin/PhaseDiagramDataGen',filename])
    print('Thermochimica calculation finished.')

    fname = 'thermoout.json'

    mint, maxt = processPhaseDiagramData(fname, el2, ts, x1, x2, p1, p2, mint, maxt, x0data, x1data)
    makePlot(el1, el2, ts, x1, x2, p1, p2, mint, maxt, labels, x0data, x1data)
    return mint, maxt

def clockwiseangle_and_distance(point):
    # Vector between point and the origin: v = p - o
    vector = [point[0]-origin[0], point[1]-origin[1]]
    # Length of vector: ||v||
    lenvector = math.hypot(vector[0], vector[1])
    # If length is zero there is no angle
    if lenvector == 0:
        return -math.pi, 0
    # Normalize vector: v/||v||
    normalized = [vector[0]/lenvector, vector[1]/lenvector]
    dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
    diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
    angle = math.atan2(diffprod, dotprod)
    # Negative angles represent counter-clockwise angles so we need to subtract them
    # from 2*pi (360 degrees)
    if angle < 0:
        return 2*math.pi+angle, lenvector
    # I return first the angle because that's the primary sorting criterium
    # but if two vectors have the same angle then the shorter distance should come first.
    return angle, lenvector

def makePlot(el1, el2, ts, x1, x2, p1, p2, mint, maxt, labels, x0data, x1data):
    boundaries = []
    phases = []
    b = []
    for i in range(len(p1)):
        # If a miscibility gap label has been used unnecessarily, remove it
        if p1[i].find('#2') > 0:
            if not(p1[i][0:p1[i].find('#2')] == p2[i]):
                p1[i] = p1[i][0:p1[i].find('#2')]
        if p2[i].find('#2') > 0:
            if not(p2[i][0:p2[i].find('#2')] == p1[i]):
                p2[i] = p2[i][0:p2[i].find('#2')]

        repeat = False
        for j in range(len(boundaries)):
            if (boundaries[j][0] == p1[i]) and (boundaries[j][1] == p2[i]):
                b.append(j)
                repeat = True
        if not(repeat):
            boundaries.append([p1[i],p2[i]])
            b.append(len(boundaries)-1)

    for i in range(len(boundaries)):
        repeat1 = False
        repeat2 = False
        for j in range(len(phases)):
            if (boundaries[i][0] == phases[j]):
                repeat1 = True
            if (boundaries[i][1] == phases[j]):
                repeat2 = True
        if not(repeat1):
            phases.append(boundaries[i][0])
        if not(repeat2):
            phases.append(boundaries[i][1])

    phasePolyPoints = [[] for i in range(len(phases))]
    # Start figure
    fig = plt.figure()
    ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])

    bEdgeLine = [[False,False] for i in range(len(boundaries))]
    # Plot along x=0 and x=1 boundaries
    for j in range(len(x0data[1])):
        i = phases.index(x0data[0][j])
        phasePolyPoints[i].append([[0,x0data[1][j]]])
        phasePolyPoints[i].append([[0,x0data[2][j]]])
        if j > 0:
            ax.plot(0,x0data[1][j],'kv')
            for k in range(len(boundaries)):
                if (x0data[0][j] in boundaries[k]) and (x0data[0][j-1] in boundaries[k]):
                    inds = [i for i, l in enumerate(b) if l == k]
                    bind = boundaries[k].index(x0data[0][j])
                    if bind == 0:
                        minj = np.argmin(np.array(x1)[inds])
                        ax.plot([0,np.array(x1)[inds][minj]],[x0data[1][j],np.array(ts)[inds][minj]],'k-')
                    elif bind == 1:
                        minj = np.argmin(np.array(x2)[inds])
                        ax.plot([0,np.array(x2)[inds][minj]],[x0data[1][j],np.array(ts)[inds][minj]],'k-')
                    if np.array(ts)[inds][minj] == np.min(np.array(ts)[inds]):
                        bEdgeLine[k][0] = True
                    if np.array(ts)[inds][minj] == np.max(np.array(ts)[inds]):
                        bEdgeLine[k][1] = True
        if j < len(x0data[1]) - 1:
            ax.plot(0,x0data[2][j],'k^')
            for k in range(len(boundaries)):
                if (x0data[0][j] in boundaries[k]) and (x0data[0][j+1] in boundaries[k]):
                    inds = [i for i, l in enumerate(b) if l == k]
                    bind = boundaries[k].index(x0data[0][j])
                    if bind == 0:
                        minj = np.argmin(np.array(x1)[inds])
                        ax.plot([0,np.array(x1)[inds][minj]],[x0data[2][j],np.array(ts)[inds][minj]],'k-')
                    elif bind == 1:
                        minj = np.argmin(np.array(x2)[inds])
                        ax.plot([0,np.array(x2)[inds][minj]],[x0data[2][j],np.array(ts)[inds][minj]],'k-')
                    if np.array(ts)[inds][minj] == np.min(np.array(ts)[inds]):
                        bEdgeLine[k][0] = True
                    if np.array(ts)[inds][minj] == np.max(np.array(ts)[inds]):
                        bEdgeLine[k][1] = True
    for j in range(len(x1data[1])):
        i = phases.index(x1data[0][j])
        phasePolyPoints[i].append([[1,x1data[1][j]]])
        phasePolyPoints[i].append([[1,x1data[2][j]]])
        if j > 0:
            ax.plot(1,x1data[1][j],'kv')
            for k in range(len(boundaries)):
                if (x1data[0][j] in boundaries[k]) and (x1data[0][j-1] in boundaries[k]):
                    inds = [i for i, l in enumerate(b) if l == k]
                    bind = boundaries[k].index(x1data[0][j])
                    if bind == 0:
                        maxj = np.argmax(np.array(x1)[inds])
                        ax.plot([1,np.array(x1)[inds][maxj]],[x1data[1][j],np.array(ts)[inds][maxj]],'k-')
                    elif bind == 1:
                        maxj = np.argmax(np.array(x2)[inds])
                        ax.plot([1,np.array(x2)[inds][maxj]],[x1data[1][j],np.array(ts)[inds][maxj]],'k-')
                    if np.array(ts)[inds][maxj] == np.min(np.array(ts)[inds]):
                        bEdgeLine[k][0] = True
                    if np.array(ts)[inds][maxj] == np.max(np.array(ts)[inds]):
                        bEdgeLine[k][1] = True
        if j < len(x0data[1]) - 1:
            ax.plot(1,x1data[2][j],'k^')
            for k in range(len(boundaries)):
                if (x1data[0][j] in boundaries[k]) and (x1data[0][j+1] in boundaries[k]):
                    inds = [i for i, l in enumerate(b) if l == k]
                    bind = boundaries[k].index(x1data[0][j])
                    if bind == 0:
                        maxj = np.argmax(np.array(x1)[inds])
                        ax.plot([1,np.array(x1)[inds][maxj]],[x1data[2][j],np.array(ts)[inds][maxj]],'k-')
                    elif bind == 1:
                        maxj = np.argmax(np.array(x2)[inds])
                        ax.plot([1,np.array(x2)[inds][maxj]],[x1data[2][j],np.array(ts)[inds][maxj]],'k-')
                    if np.array(ts)[inds][maxj] == np.min(np.array(ts)[inds]):
                        bEdgeLine[k][0] = True
                    if np.array(ts)[inds][maxj] == np.max(np.array(ts)[inds]):
                        bEdgeLine[k][1] = True

    outline = Polygon([[0,mint], [0, maxt], [1, maxt], [1, mint]])
    # plot 2-phase region boundaries
    for j in range(len(boundaries)):
        polygonPoints = []
        inds = [i for i, k in enumerate(b) if k == j]
        ttt = np.array(ts)[inds]
        sindex = np.argsort(ttt)
        ttt = ttt[sindex]
        x1t = np.array(x1)[inds]
        x1t = x1t[sindex]
        x2t = np.array(x2)[inds]
        x2t = x2t[sindex]
        ax.plot(np.array(x1)[inds],np.array(ts)[inds],'.')
        ax.plot(np.array(x2)[inds],np.array(ts)[inds],'.')
        for i in range(len(inds)):
            polygonPoints.append([x1t[i],ttt[i]])
        for i in reversed(range(len(inds))):
            polygonPoints.append([x2t[i],ttt[i]])
        phaseOutline = Polygon(polygonPoints).buffer(0)
        outline = outline - phaseOutline
        # patch = PolygonPatch(phaseOutline.buffer(0))
        # ax.add_patch(patch)
        minj = np.argmin(np.array(ts)[inds])
        maxj = np.argmax(np.array(ts)[inds])
        # plot invariant temperatures
        if (np.array(ts)[inds][minj] > mint) and not(bEdgeLine[j][0]):
            ax.plot([np.array(x1)[inds][minj],np.array(x2)[inds][minj]],[np.array(ts)[inds][minj],np.array(ts)[inds][minj]],'k-')
        if (np.array(ts)[inds][maxj] < maxt) and not(bEdgeLine[j][1]):
            ax.plot([np.array(x1)[inds][maxj],np.array(x2)[inds][maxj]],[np.array(ts)[inds][maxj],np.array(ts)[inds][maxj]],'k-')
        for i in range(len(phases)):
            if boundaries[j][0] == phases[i]:
                phasePolyPoints[i].append(polygonPoints[:len(inds)])
            if boundaries[j][1] == phases[i]:
                phasePolyPoints[i].append(list(reversed(polygonPoints))[:len(inds)])
    for i in range(len(phases)):
        segcenters = []
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
        phaseOutline = Polygon(sortedPolyPoints).buffer(0)
        outline = outline - phaseOutline

    resolution = 100
    x, y = np.meshgrid(np.arange(0, 1, 1/resolution), np.arange(mint, maxt, (maxt-mint)/resolution))
    points = MultiPoint(list(zip(x.flatten(),y.flatten()+x.flatten()*(maxt-mint)/resolution)))
    prepOutline = prep(outline)
    valid_points = list(filter(prepOutline.contains,points))
    xs = [point.x for point in valid_points]
    ys = [point.y for point in valid_points]
    ax.plot(xs,ys,'k*')

    # patch = PolygonPatch(outline.buffer(0))
    # ax.add_patch(patch)
    ax.set_xlim(0,1)
    ax.set_ylim(mint,maxt)
    ax.set_title(str(el1) + ' + ' + str(el2) + ' binary phase diagram')
    ax.set_xlabel('Mole fraction ' + str(el2))
    ax.set_ylabel('Temperature [K]')
    for lab in labels:
        plt.text(float(lab[0][0]),float(lab[0][1]),lab[1])
    plt.show()

def writeInputFile(filename,xlo,xhi,nxstep,tlo,thi,ntstep,pressure,tunit,punit,munit,el1,el2,datafile):
    with open(filename, 'w') as inputFile:
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
        inputFile.write('pressure          = ' + str(pressure) + '\n')
        inputFile.write('temperature unit         = ' + tunit + '\n')
        inputFile.write('pressure unit          = ' + punit + '\n')
        inputFile.write('mass unit          = \'' + munit + '\'\n')
        inputFile.write('iEl         = ' + str(atomic_number_map.index(el1)+1) + ' ' + str(atomic_number_map.index(el2)+1) + '\n')
        inputFile.write('data file         = ' + datafile + '\n')

def addLabel(filename,xlab,tlab,pressure,tunit,punit,munit,el1,el2,datafile,mint,maxt,labels,x0data,x1data,ts,x1,x2,p1,p2):
    writeInputFile(filename,xlab,xlab,0,tlab,tlab,0,pressure,tunit,punit,munit,el1,el2,datafile)
    subprocess.run(['./bin/PhaseDiagramDataGen',filename])
    fname = 'thermoout.json'
    f = open(fname,)
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
    labels.append([[xlab,tlab],'+'.join(labelName)])
    mint, maxt = runCalc(el1, el2, ts, x1, x2, p1, p2, mint, maxt, labels, x0data, x1data)
    return mint, maxt

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

dataWindow = sg.Window('Thermochimica database selection', file_list_column, location = [0,0], finalize=True)
folder = os.getcwd()+'/data'
try:
    file_list = os.listdir(folder)
except:
    file_list = []

fnames = [
    f
    for f in file_list
    if os.path.isfile(os.path.join(folder, f))
    and f.lower().endswith((".dat", ".DAT"))
]
dataWindow["-FILE LIST-"].update(fnames)

timeout = 50
inputSize = 20

while True:
    event, values = dataWindow.read()
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    elif event == "-FOLDER-":
        folder = values["-FOLDER-"]
        try:
            file_list = os.listdir(folder)
        except:
            file_list = []
        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(folder, f))
            and f.lower().endswith((".dat", ".DAT"))
        ]
        dataWindow["-FILE LIST-"].update(fnames)
    elif event == "-FILE LIST-":  # A file was chosen from the listbox
        try:
            datafile = os.path.join(
                folder, values["-FILE LIST-"][0]
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
                    print(elements[i]+' not in list') # if the name is bogus (or e(phase)), discard
                    elements.remove(elements[i])
                    nElements = nElements - 1
            elSelectLayout = [sg.Column([[sg.Text('Element 1')],[sg.Combo(elements[:nElements],default_value=elements[0],key='-el1-')]],vertical_alignment='t'),
                              sg.Column([[sg.Text('Element 2')],[sg.Combo(elements[:nElements],default_value=elements[1],key='-el2-')]],vertical_alignment='t')]
            xLayout    = [sg.Column([[sg.Text('Start Element 2 Concentration')],[sg.Input(key='-xlo-',size=(inputSize,1))],
                          [sg.Text('Concentration unit')],[sg.Combo(['mole fraction'],default_value='mole fraction',key='-munit-')]],vertical_alignment='t'),
                          sg.Column([[sg.Text('End Element 2 Concentration')],[sg.Input(key='-xhi-',size=(inputSize,1))],
                          ],vertical_alignment='t'),
                          sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstep-',size=(8,1))]],vertical_alignment='t')]
            tempLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperature-',size=(inputSize,1))],
                          [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')]],vertical_alignment='t'),
                          sg.Column([[sg.Text('End Temperature')],[sg.Input(key='-endtemperature-',size=(inputSize,1))],
                          ],vertical_alignment='t'),
                          sg.Column([[sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstep-',size=(8,1))]],vertical_alignment='t')]
            presLayout = [sg.Column([[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(inputSize,1))],
                          [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')]],vertical_alignment='t')
                          ]
            setupLayout = [elSelectLayout,xLayout,tempLayout,presLayout,[sg.Button('Run'), sg.Button('Refine', disabled = True),
                            sg.Button('Add Label', disabled = True), sg.Button('Remove Label', disabled = True), sg.Exit()]]
            setupWindow = sg.Window('Phase diagram setup', setupLayout, location = [400,0], finalize=True)
            while True:
                event, values = setupWindow.read(timeout=timeout)
                eventd, valuesd = dataWindow.read(timeout=timeout)
                if event == sg.WIN_CLOSED or event == 'Exit' or eventd == sg.WIN_CLOSED or eventd == 'Exit':
                    break
                elif event =='Run':
                    cancelRun = False
                    ntstep = values['-ntstep-']
                    nxstep = values['-nxstep-']
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
                    tlo= values['-temperature-']
                    thi = values['-endtemperature-']
                    pressure = values['-pressure-']
                    filename = 'inputs/pythonPhaseDiagramInput.ti'
                    tunit = values['-tunit-']
                    punit = values['-punit-']
                    xhi = values['-xhi-']
                    xlo = values['-xlo-']
                    el1 = values['-el1-']
                    el2 = values['-el2-']
                    if (str(el1) == str(el2)) or (float(tlo) == float(thi)):
                        cancelRun = True
                        repeatLayout = [[sg.Text('Values cannot be equal.')],[sg.Button('Cancel')]]
                        repeatWindow = sg.Window('Repeat value notification', repeatLayout, location = [400,0], finalize=True)
                        while True:
                            event, values = repeatWindow.read(timeout=timeout)
                            if event == sg.WIN_CLOSED or event == 'Cancel':
                                break
                            elif event == 'Continue':
                                cancelRun = False
                                break
                        repeatWindow.close()
                    munit = values['-munit-']
                    if cancelRun:
                        continue
                    writeInputFile(filename,xlo,xhi,nxstep,tlo,thi,ntstep,pressure,tunit,punit,munit,el1,el2,datafile)
                    ts = []
                    x1 = []
                    x2 = []
                    p1 = []
                    p2 = []
                    x0data = [[],[],[]]
                    x1data = [[],[],[]]
                    mint = 1e6
                    maxt = 0
                    labels = []
                    mint, maxt = runCalc(el1, el2, ts, x1, x2, p1, p2, mint, maxt, labels, x0data, x1data)
                    setupWindow.Element('Refine').Update(disabled = False)
                    setupWindow.Element('Add Label').Update(disabled = False)
                elif event =='Refine':
                    xRefLayout    = [sg.Column([[sg.Text('Start Concentration')],[sg.Input(key='-xlor-',size=(inputSize,1))]],vertical_alignment='t'),
                                  sg.Column([[sg.Text('End Concentration')],[sg.Input(key='-xhir-',size=(inputSize,1))]],vertical_alignment='t'),
                                  sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstepr-',size=(8,1))]],vertical_alignment='t')]
                    tempRefLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperaturer-',size=(inputSize,1))]],vertical_alignment='t'),
                                  sg.Column([[sg.Text('End Temperature')],[sg.Input(key='-endtemperaturer-',size=(inputSize,1))]],vertical_alignment='t'),
                                  sg.Column([[sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstepr-',size=(8,1))]],vertical_alignment='t')]
                    refineLayout = [xRefLayout,tempRefLayout,[sg.Button('Refine'), sg.Button('Cancel')]]
                    refineWindow = sg.Window('Phase diagram refinement', refineLayout, location = [400,0], finalize=True)
                    while True:
                        event, values = refineWindow.read(timeout=timeout)
                        if event == sg.WIN_CLOSED or event == 'Cancel':
                            break
                        elif event =='Refine':
                            cancelRun = False
                            ntstep = values['-ntstepr-']
                            nxstep = values['-nxstepr-']
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
                            tlo = values['-temperaturer-']
                            thi = values['-endtemperaturer-']
                            xhi = values['-xhir-']
                            xlo = values['-xlor-']
                            if cancelRun:
                                continue
                            writeInputFile(filename,xlo,xhi,nxstep,tlo,thi,ntstep,pressure,tunit,punit,munit,el1,el2,datafile)
                            mint, maxt = runCalc(el1, el2, ts, x1, x2, p1, p2, mint, maxt, labels, x0data, x1data)
                    refineWindow.close()
                elif event =='Add Label':
                    xLabLayout    = [[sg.Text('Element 2 Concentration')],[sg.Input(key='-xlab-',size=(inputSize,1))]]
                    tLabLayout = [[sg.Text('Temperature')],[sg.Input(key='-tlab-',size=(inputSize,1))]]
                    labelLayout = [xLabLayout,tLabLayout,[sg.Button('Add Label'), sg.Button('Cancel')]]
                    labelWindow = sg.Window('Add phase label', labelLayout, location = [400,0], finalize=True)
                    while True:
                        event, values = labelWindow.read(timeout=timeout)
                        if event == sg.WIN_CLOSED or event == 'Cancel':
                            break
                        elif event =='Add Label':
                            xlab = values['-xlab-']
                            tlab = values['-tlab-']
                            mint, maxt = addLabel(filename,xlab,tlab,pressure,tunit,punit,munit,el1,el2,datafile,mint,maxt,labels,x0data,x1data,ts,x1,x2,p1,p2)
                    labelWindow.close()
                    setupWindow.Element('Remove Label').Update(disabled = False)
                elif event =='Remove Label':
                    headingsLayout = [[sg.Text('Label Text',   size = [55,1],justification='left'),
                                       sg.Text('Concentration',size = [15,1],justification='center'),
                                       sg.Text('Temperature',  size = [15,1],justification='center'),
                                       sg.Text('Remove Label?',size = [15,1])]]
                    labelListLayout = []
                    for i in range(len(labels)):
                        labelListLayout.append([[sg.Text(labels[i][1],size = [55,1],justification='left'),
                                                 sg.Text(labels[i][0][0],size = [15,1],justification='center'),
                                                 sg.Text(labels[i][0][1],size = [15,1],justification='center'),
                                                 sg.Checkbox('',key='-removeLabel'+str(i)+'-',pad=[[40,0],[0,0]])]])
                    removeLayout = [headingsLayout,labelListLayout,[sg.Button('Remove Label(s)'), sg.Button('Cancel')]]
                    removeWindow = sg.Window('Add phase label', removeLayout, location = [400,0], finalize=True)
                    while True:
                        event, values = removeWindow.read(timeout=timeout)
                        if event == sg.WIN_CLOSED or event == 'Cancel':
                            break
                        if event == 'Remove Label(s)':
                            tempLength = len(labels)
                            for i in reversed(range(tempLength)):
                                if values['-removeLabel'+str(i)+'-']:
                                    del labels[i]
                            if len(labels) == 0:
                                setupWindow.Element('Remove Label').Update(disabled = True)
                            break
                    removeWindow.close()
                    makePlot(el1, el2, ts, x1, x2, p1, p2, mint, maxt, labels, x0data, x1data)
            setupWindow.close()
        except:
            pass
