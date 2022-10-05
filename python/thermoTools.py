import subprocess
import matplotlib.pyplot as plt
import numpy as np
import json

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

def WriteRunCalculationList(filename,datafile,elements,calcList,tunit='K',punit='atm',munit='moles',printMode=2,heatCapacity=False,writeJson=True,debugMode=False,reinitialization=False,minSpecies=None):
    nElements = len(elements)
    with open(filename, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        inputFile.write(f'data file         = {datafile}\n')
        inputFile.write(f'temperature unit  = {tunit}\n')
        inputFile.write(f'pressure unit     = {punit}\n')
        inputFile.write(f'mass unit         = \'{munit}\'\n')
        inputFile.write(f'print mode        = {printMode}\n')
        inputFile.write(f'heat capacity     = {".TRUE." if heatCapacity else ".FALSE."}\n')
        inputFile.write(f'write json        = {".TRUE." if writeJson else ".FALSE."}\n')
        inputFile.write(f'debug mode        = {".TRUE." if debugMode else ".FALSE."}\n')
        inputFile.write(f'reinitialization  = {".TRUE." if reinitialization else ".FALSE."}\n')
        if minSpecies:
            # Preserve Thermochimica default unless set
            inputFile.write(f'min species       = {minSpecies}\n')
        inputFile.write(f'nEl               = {nElements} \n')
        inputFile.write(f'iEl               = {" ".join([str(atomic_number_map.index(elem)+1) for elem in elements])}\n')
        inputFile.write(f'nCalc             = {len(calcList)}\n')
        # Write calculations list
        for calc in calcList:
            inputFile.write(f'{calc[0]} {calc[1]} {" ".join([str(calc[i]) for i in range(2,len(calc))])}\n')

def WriteInputScript(filename,datafile,elements,tstart,tend,ntstep,pstart,pend,npstep,masses,tunit='K',punit='atm',munit='moles',printMode=2,heatCapacity=False,writeJson=True,debugMode=False,reinitialization=False,minSpecies=None,stepTogether=False):
    nElements = len(elements)
    with open(filename, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        if ntstep > 0:
            tstep = (float(tend)-float(tstart))/ntstep
            inputFile.write(f'temperature       = {tstart}:{tend}:{tstep}\n')
        else:
            inputFile.write(f'temperature       = {tstart}\n')
        if npstep > 0:
            pstep = (float(pend)-float(pstart))/npstep
            inputFile.write(f'pressure          = {pstart}:{pend}:{pstep}\n')
        else:
            inputFile.write(f'pressure          = {pstart}\n')
        for i in range(nElements):
            inputFile.write(f'mass(' + str(atomic_number_map.index(elements[i])+1) + ')           = ' + str(masses[i]) + '\n')
        inputFile.write(f'temperature unit  = {tunit}\n')
        inputFile.write(f'pressure unit     = {punit}\n')
        inputFile.write(f'step together     = {".TRUE." if stepTogether else ".FALSE."}\n')
        inputFile.write(f'mass unit         = {munit}\n')
        inputFile.write(f'data file         = {datafile}\n')
        inputFile.write(f'print mode        = {printMode}\n')
        inputFile.write(f'heat capacity     = {".TRUE." if heatCapacity else ".FALSE."}\n')
        inputFile.write(f'debug mode        = {".TRUE." if debugMode else ".FALSE."}\n')
        inputFile.write(f'write json        = {".TRUE." if writeJson else ".FALSE."}\n')
        inputFile.write(f'debug mode        = {".TRUE." if debugMode else ".FALSE."}\n')
        inputFile.write(f'reinitialization  = {".TRUE." if reinitialization else ".FALSE."}\n')
        if minSpecies:
            # Preserve Thermochimica default unless set
            inputFile.write(f'min species       = {minSpecies}\n')

def RunRunCalculationList(filename,checkOutput=False):
    if checkOutput:
        thermoOut = subprocess.check_output(['./bin/RunCalculationList',filename]).decode("utf-8")
        return thermoOut
    else:
        subprocess.run(['./bin/RunCalculationList',filename])

def RunInputScript(filename,checkOutput=False):
    if checkOutput:
        thermoOut = subprocess.check_output(['./bin/InputScriptMode',filename]).decode("utf-8")
        return thermoOut
    else:
        subprocess.run(['./bin/InputScriptMode',filename])

def makePlot(datafile,xkey,yen,ykey,leg,ylab,yen2=None,ykey2=None,leg2=None,ylab2=None,plotColor='colorful',plotColor2='colorful',plotMarker='-.',plotMarker2='--*',xlog=False,ylog=False,ylog2=False):
    # Do plot setup
    x,y,y2,xlab,legend,legend2 = plotDataSetup(datafile,xkey,yen,ykey,leg,yen2=yen2,ykey2=ykey2,leg2=leg2)

    # Start figure
    fig = plt.figure()
    plt.ion()
    lns=[]
    if True in yen2:
        ax = fig.add_axes([0.2, 0.1, 0.65, 0.85])
    else:
        ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(y))))
    for yi in range(len(y)):
        if plotColor == 'colorful':
            c = next(color)
        else:
            c = 'k'
        lns = lns + ax.plot(x,y[yi],plotMarker,c=c,label=legend[yi])
    ax.set_xlabel(xlab)
    if xlog:
        ax.set_xscale('log')
    ax.set_ylabel(ylab)
    if True in yen2:
        ax2 = ax.twinx()
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(y2))))
        for yi in range(len(y2)):
            if plotColor2 == 'colorful':
                c = next(color)
            else:
                c = 'k'
            lns = lns + ax2.plot(x,y2[yi],plotMarker2,c=c,label=legend2[yi])
        ax2.set_ylabel(ylab2)
        if ylog2:
            ax2.set_yscale('log')
    labs = [l.get_label() for l in lns]
    if ylog:
        ax.set_yscale('log')
    ax.legend(lns, labs, loc=0)
    plt.show()
    plt.pause(0.001)

    return x, y, y2, legend, legend2, xlab

def plotDataSetup(datafile,xkey,yen,ykey,leg,yen2=None,ykey2=None,leg2=None):
    # Init x-axis
    x = []
    
    # Init left-hand y-axis
    yused = []
    legend = []
    for i in range(len(yen)):
        if yen[i]:
            yused.append(ykey[i])
            legend.append(leg[i])
    y = [[] for _ in range(len(yused))]
    
    # Init right-hand y-axis
    yused2 = []
    legend2 = []
    for i in range(len(yen2)):
        if yen2[i]:
            yused2.append(ykey2[i])
            legend2.append(leg2[i])
    y2 = [[] for _ in range(len(yused2))]
    
    # Loop over all calculations and get requested values
    data = readDatabase(datafile)
    for j in data.keys():
        try:
            for yi in range(len(yused)):
                if len(yused[yi]) == 1:
                    y[yi].append(data[j][yused[yi][0]])
                elif len(yused[yi]) == 3:
                    y[yi].append(data[j][yused[yi][0]][yused[yi][1]][yused[yi][2]])
                elif len(yused[yi]) == 5:
                    if yused[yi][4] == 'vapor pressure':
                        y[yi].append(data[j][yused[yi][0]][yused[yi][1]][yused[yi][2]][yused[yi][3]]['mole fraction']*data[j]['pressure'])
                    else:
                        y[yi].append(data[j][yused[yi][0]][yused[yi][1]][yused[yi][2]][yused[yi][3]][yused[yi][4]])
            for yi in range(len(yused2)):
                if len(yused2[yi]) == 1:
                    y2[yi].append(data[j][yused2[yi][0]])
                elif len(yused2[yi]) == 3:
                    y2[yi].append(data[j][yused2[yi][0]][yused2[yi][1]][yused2[yi][2]])
                elif len(yused2[yi]) == 5:
                    if yused2[yi][4] == 'vapor pressure':
                        y2[yi].append(data[j][yused2[yi][0]][yused2[yi][1]][yused2[yi][2]][yused2[yi][3]]['mole fraction']*data[j]['pressure'])
                    else:
                        y2[yi].append(data[j][yused2[yi][0]][yused2[yi][1]][yused2[yi][2]][yused2[yi][3]][yused2[yi][4]])
            if xkey == 'iteration':
                x.append(int(j))
                xlab = 'Iteration'
            else:
                x.append(data[j][xkey])
                if xkey == 'temperature':
                    xlab = 'Temperature [K]'
                elif xkey == 'pressure':
                    xlab = 'Pressure [atm]'
        except:
            # do nothing
            continue
    
    return x,y,y2,xlab,legend,legend2

def exportPlotScript(filename,datafile,xkey,yen,ykey,leg,ylab,yen2=None,ykey2=None,leg2=None,ylab2=None,plotColor='colorful',plotColor2='colorful',plotMarker='-.',plotMarker2='--*',xlog=False,ylog=False,ylog2=False):
    # Don't want to call makePlot() from here because then it is too hard to reconfigure the plot
    # Call plotDataSetup() instead
    with open(filename, 'w') as f:
        f.write('# Thermochimica-generated plot script\n')
        # Imports
        f.write('import matplotlib.pyplot as plt\n')
        f.write('import numpy as np\n')
        f.write('import thermoTools\n')
        f.write('\n')
        # Copy data
        f.write(f'datafile = \'{datafile}\'\n')
        f.write(f'xkey     = \'{xkey}\'\n')
        f.write(f'yen      = {yen}\n')
        f.write(f'ykey     = {ykey}\n')
        f.write(f'leg      = {leg}\n')
        f.write(f'yen2     = {yen2}\n')
        f.write(f'ykey2    = {ykey2}\n')
        f.write(f'leg2     = {leg2}\n')
        f.write(f'ylab     = \'{ylab}\'\n')
        # Call plotDataSetup (enables use of datafile path)
        f.write('x,y,y2,xlab,legend,legend2 = thermoTools.plotDataSetup(datafile,xkey,yen,ykey,leg,yen2=yen2,ykey2=ykey2,leg2=leg2)\n')
        # Figure generation commands
        f.write('lns=[]\n')
        f.write('# Start figure\n')
        f.write('fig = plt.figure()\n')
        if True in yen2:
            f.write('ax  = fig.add_axes([0.2, 0.1, 0.65, 0.85])\n')
        else:
            f.write('ax  = fig.add_axes([0.2, 0.1, 0.75, 0.85])\n')
        f.write('color = iter(plt.cm.rainbow(np.linspace(0, 1, len(y))))\n')
        f.write('for yi in range(len(y)):\n')
        if plotColor == 'colorful':
            f.write('    c = next(color)\n')
        else:
            f.write('    c = \'k\'\n')
        f.write(f'    lns = lns + ax.plot(x,y[yi],\'{plotMarker}\',c=c,label=legend[yi])\n')
        if True in yen2:
            f.write(f'ylab2 = \'{ylab2}\'\n')
            f.write('ax2 = ax.twinx()\n')
            f.write('color = iter(plt.cm.rainbow(np.linspace(0, 1, len(y2))))\n')
            f.write('for yi in range(len(y2)):\n')
            if plotColor2 == 'colorful':
                f.write('    c = next(color)\n')
            else:
                f.write('    c = \'k\'\n')
            f.write(f'    lns = lns + ax2.plot(x,y2[yi],\'{plotMarker2}\',c=c,label=legend2[yi])\n')
            f.write('ax2.set_ylabel(ylab2)\n')
            if ylog2:
                f.write("ax2.set_yscale('log')\n")
        f.write('labs = [l.get_label() for l in lns]\n')
        f.write('ax.legend(lns, labs, loc=0)\n')
        f.write('ax.set_xlabel(xlab)\n')
        f.write('ax.set_ylabel(ylab)\n')
        if xlog:
            f.write("ax.set_xscale('log')\n")
        if ylog:
            f.write("ax.set_yscale('log')\n")
        f.write('plt.show()\n')

def readDatabase(datafile):
    f = open(datafile,)
    data = json.load(f)
    f.close()
    return data

