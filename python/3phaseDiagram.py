import PySimpleGUI as sg
import json
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import math
import os
import subprocess

def processPhaseDiagramData(fname, elx, ely, x1, x2, p1, p2, x0data, x1data):
    f = open(fname,)
    data = json.load(f)
    f.close()
    if list(data.keys())[0] != '1':
        print('Output does not contain data series')
        exit()
    for i in list(data.keys()):
        if (data[i]['# solution phases'] + data[i]['# pure condensed phases']) == 2:
            boundPhases = []
            boundComps = []
            for phaseType in ['solution phases','pure condensed phases']:
                for phaseName in list(data[i][phaseType].keys()):
                    if (data[i][phaseType][phaseName]['moles'] > 0):
                        boundPhases.append(phaseName)
                        tempComps = [0,0]
                        if elx in list(data[i][phaseType][phaseName]['elements'].keys()):
                            tempComps[0] = data[i][phaseType][phaseName]['elements'][elx]['mole fraction of phase by element']
                        if ely in list(data[i][phaseType][phaseName]['elements'].keys()):
                            tempComps[1] = data[i][phaseType][phaseName]['elements'][ely]['mole fraction of phase by element']
                        boundComps.append(tempComps)
            x1.append(boundComps[0])
            x2.append(boundComps[1])
            p1.append(boundPhases[0])
            p2.append(boundPhases[1])

def runCalc(el1, el2, el3, x1, x2, p1, p2, labels, x0data, x1data):
    print('Thermochimica calculation initiated.')
    subprocess.run(['./bin/Phase3DiagramDataGen',filename])
    print('Thermochimica calculation finished.')

    fname = 'thermoout.json'

    processPhaseDiagramData(fname, el1, el2, x1, x2, p1, p2, x0data, x1data)
    makePlot(el1, el2, el3, x1, x2, p1, p2, labels, x0data, x1data)

def fmt(x,pos=None):
    return '{:.2f}'.format(1-x)

def makePlot(el1, el2, el3, x1, x2, p1, p2, labels, x0data, x1data):
    boundaries = []
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

    # Start figure
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.85])

    # plot 2-phase region boundaries
    for j in range(len(boundaries)):
        inds = [i for i, k in enumerate(b) if k == j]
        ax.plot(1-(np.array(x1)[inds,0]+np.array(x1)[inds,1]/2),np.array(x1)[inds,1],'.')
        ax.plot(1-(np.array(x2)[inds,0]+np.array(x2)[inds,1]/2),np.array(x2)[inds,1],'.')
        minj = np.argmin(np.array(x1)[inds,0])
        maxj = np.argmax(np.array(x1)[inds,0])
        if minj == maxj:
            minj = np.argmin(np.array(x2)[inds,0])
            maxj = np.argmax(np.array(x2)[inds,0])
        # plot 2-phase to 3-phase boundaries
        ax.plot([1-(np.array(x1)[inds,0][minj]+np.array(x1)[inds,1][minj]/2),1-(np.array(x2)[inds,0][minj]+np.array(x2)[inds,1][minj]/2)],
        [np.array(x1)[inds,1][minj],np.array(x2)[inds,1][minj]],'k-')
        ax.plot([1-(np.array(x1)[inds,0][maxj]+np.array(x1)[inds,1][maxj]/2),1-(np.array(x2)[inds,0][maxj]+np.array(x2)[inds,1][maxj]/2)],
        [np.array(x1)[inds,1][maxj],np.array(x2)[inds,1][maxj]],'k-')

    ax.plot([0,0.5,1,0],[0,1,0,0],'k-')
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_title(str(el1) + ' + ' + str(el2) + ' + ' + str(el3) + ' ternary phase diagram')
    ax.set_xlabel('Mole fraction ' + str(el1))
    ax.set_ylabel('Mole fraction ' + str(el2))
    for lab in labels:
        plt.text(float(lab[0][0]),float(lab[0][1]),lab[1])
    # reverse x tick labels for ternary plot
    ax.xaxis.set_major_formatter(FuncFormatter(fmt))
    ax2 = ax.twinx()
    ax2.yaxis.set_major_formatter(FuncFormatter(fmt))
    ax2.set_ylabel('Mole fraction ' + str(el3))
    plt.sca(ax)
    plt.show()

def writeInputFile(filename,x1lo,x1hi,x2lo,x2hi,nx1step,nx2step,t,pressure,tunit,punit,munit,el1,el2,el3,datafile):
    with open(filename, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        if float(nx1step) > 0:
            x1step = (float(x1hi)-float(x1lo))/float(nx1step)
        else:
            x1step = 0
        if float(nx2step) > 0:
            x2step = (float(x2hi)-float(x2lo))/float(nx2step)
        else:
            x2step = 0
        inputFile.write('x1          = ' + str(x1lo) + ':' + str(x1hi) + ':' + str(x1step) + '\n')
        inputFile.write('x2          = ' + str(x2lo) + ':' + str(x2hi) + ':' + str(x2step) + '\n')
        inputFile.write('temperature          = ' + str(t) + '\n')
        inputFile.write('pressure          = ' + str(pressure) + '\n')
        inputFile.write('temperature unit         = ' + tunit + '\n')
        inputFile.write('pressure unit          = ' + punit + '\n')
        inputFile.write('mass unit          = \'' + munit + '\'\n')
        inputFile.write('iEl         = ' + str(atomic_number_map.index(el1)+1) + ' ' + str(atomic_number_map.index(el2)+1) + ' ' + str(atomic_number_map.index(el3)+1) + '\n')
        inputFile.write('data file         = ' + datafile + '\n')

def addLabel(filename,x1lab,x2lab,temperature,pressure,tunit,punit,munit,el1,el2,datafile,labels,x0data,x1data,x1,x2,p1,p2):
    writeInputFile(filename,x1lab,x1lab,x2lab,x2lab,0,0,temperature,pressure,tunit,punit,munit,el1,el2,el3,datafile)
    subprocess.run(['./bin/Phase3DiagramDataGen',filename])
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
    labels.append([[1-(x1lab+x2lab/2),x2lab],'+'.join(labelName)])
    runCalc(el1, el2, el3, x1, x2, p1, p2, labels, x0data, x1data)

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

# fname = 'thermoout.json'

filename = 'inputs/python3PhaseDiagramInput.ti'
datafile = 'data/Kaye_NobleMetals.dat'
el1 = 'Mo'
el2 = 'Pd'
el3 = 'Ru'
x1 = []
x2 = []
p1 = []
p2 = []
x0data = [[],[],[]]
x1data = [[],[],[]]
labels = []

x1lo = 0
x1hi = 1
x2lo = 0
x2hi = 1
nxstep = 30
t = 900
pressure = 1
tunit = 'K'
punit = 'atm'
munit = 'mole fraction'
writeInputFile(filename,x1lo,x1hi,x2lo,x2hi,nxstep,nxstep,t,pressure,tunit,punit,munit,el1,el2,el3,datafile)
runCalc(el1, el2, el3, x1, x2, p1, p2, labels, x0data, x1data)
x1lab = 0.7
x2lab = 0.2
addLabel(filename,x1lab,x2lab,t,pressure,tunit,punit,munit,el1,el2,datafile,labels,x0data,x1data,x1,x2,p1,p2)
