import subprocess
import matplotlib.pyplot as plt
import numpy as np
import json
import shutil

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

def WriteRunCalculationList(filename,datafile,elements,calcList,tunit='K',punit='atm',munit='moles',printMode=2,heatCapacity=False,writeJson=True,debugMode=False,reinitialization=False,minSpecies=None,excludePhases=None,excludePhasesExcept=None,fuzzyStoichiometry=False,fuzzyMagnitude=-1,gibbsMinCheck=False):
    nElements = len(elements)
    with open(filename, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        inputFile.write(f'data file         = {datafile}\n')
        inputFile.write(f'temperature unit  = \'{tunit}\'\n')
        inputFile.write(f'pressure unit     = \'{punit}\'\n')
        inputFile.write(f'mass unit         = \'{munit}\'\n')
        # 0: no output, 1: condensed output, 2: standard output
        inputFile.write(f'print mode        = {printMode}\n')
        # Toggle for running properties requiring derivatives: enthalpy, entropy, heat capacity
        inputFile.write(f'heat capacity     = {".TRUE." if heatCapacity else ".FALSE."}\n')
        # Toggle for writing JSON output database
        inputFile.write(f'write json        = {".TRUE." if writeJson else ".FALSE."}\n')
        # Outputs a huge amount of information to terminal
        inputFile.write(f'debug mode        = {".TRUE." if debugMode else ".FALSE."}\n')
        # Toggle reinitialization mode
        inputFile.write(f'reinitialization  = {".TRUE." if reinitialization else ".FALSE."}\n')
        # Sets minimum number of species in solution phase for it to be included in the calculation (traditionally 2, 1 also makes sense)
        # Preserve Thermochimica default unless set
        if minSpecies:
            inputFile.write(f'min species       = {minSpecies}\n')
        # Number of elements to be used
        inputFile.write(f'nEl               = {nElements} \n')
        # List of elements (by number)
        inputFile.write(f'iEl               = {" ".join([str(atomic_number_map.index(elem)+1) for elem in elements])}\n')
        # Discard phases by name from database
        if excludePhases:
            inputFile.write(f'number excluded   = {len(excludePhases)}\n')
            inputFile.write(f'phases excluded   = {" ".join(excludePhases)}\n')
        # Discard all phases except these
        if excludePhasesExcept:
            inputFile.write(f'number excluded except = {len(excludePhasesExcept)}\n')
            inputFile.write(f'phases excluded except = {" ".join(excludePhasesExcept)}\n')
        # Fuzzy stoichiometry settings
        inputFile.write(f'fuzzy             = {".TRUE." if fuzzyStoichiometry else ".FALSE."}\n')
        if (fuzzyStoichiometry and (fuzzyMagnitude >= 0)):
            inputFile.write(f'fuzzy magnitude   = {fuzzyMagnitude}\n')
        inputFile.write(f'gibbs min         = {".TRUE." if gibbsMinCheck else ".FALSE."}\n')

        # Number of calculations to be run in list
        inputFile.write(f'nCalc             = {len(calcList)}\n')
        # Write calculations list
        for calc in calcList:
            inputFile.write(f'{calc[0]} {calc[1]} {" ".join([str(calc[i]) for i in range(2,len(calc))])}\n')

def WriteInputScript(filename,datafile,elements,tstart,tend,ntstep,pstart,pend,npstep,masses,tunit='K',punit='atm',munit='moles',printMode=2,heatCapacity=False,writeJson=True,debugMode=False,reinitialization=False,minSpecies=None,stepTogether=False,excludePhases=None,excludePhasesExcept=None,fuzzyStoichiometry=False,fuzzyMagnitude=-1,gibbsMinCheck=False):
    nElements = len(elements)
    with open(filename, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        # Temperature stepping
        if ntstep > 0:
            tstep = (float(tend)-float(tstart))/ntstep
            inputFile.write(f'temperature       = {tstart}:{tend}:{tstep}\n')
        else:
            inputFile.write(f'temperature       = {tstart}\n')
        # Pressure stepping
        if npstep > 0:
            pstep = (float(pend)-float(pstart))/npstep
            inputFile.write(f'pressure          = {pstart}:{pend}:{pstep}\n')
        else:
            inputFile.write(f'pressure          = {pstart}\n')
        # List of elements (by number)
        for i in range(nElements):
            inputFile.write(f'mass(' + str(atomic_number_map.index(elements[i])+1) + ')           = ' + str(masses[i]) + '\n')
        inputFile.write(f'temperature unit  = \'{tunit}\'\n')
        inputFile.write(f'pressure unit     = \'{punit}\'\n')
        inputFile.write(f'mass unit         = \'{munit}\'\n')
        # Toggles whether loops over temperature and pressure are nested (false) or simultaneous (true)
        inputFile.write(f'step together     = {".TRUE." if stepTogether else ".FALSE."}\n')
        inputFile.write(f'data file         = {datafile}\n')
        # 0: no output, 1: condensed output, 2: standard output
        inputFile.write(f'print mode        = {printMode}\n')
        # Toggle for running properties requiring derivatives: enthalpy, entropy, heat capacity
        inputFile.write(f'heat capacity     = {".TRUE." if heatCapacity else ".FALSE."}\n')
        # Outputs a huge amount of information to terminal
        inputFile.write(f'debug mode        = {".TRUE." if debugMode else ".FALSE."}\n')
        # Toggle for writing JSON output database
        inputFile.write(f'write json        = {".TRUE." if writeJson else ".FALSE."}\n')
        # Toggle reinitialization mode
        inputFile.write(f'reinitialization  = {".TRUE." if reinitialization else ".FALSE."}\n')
        # Discard phases by name from database
        if excludePhases:
            inputFile.write(f'number excluded   = {len(excludePhases)}\n')
            inputFile.write(f'phases excluded   = {" ".join(excludePhases)}\n')
        # Discard all phases except these
        if excludePhasesExcept:
            inputFile.write(f'number excluded except = {len(excludePhasesExcept)}\n')
            inputFile.write(f'phases excluded except = {" ".join(excludePhasesExcept)}\n')
        # Sets minimum number of species in solution phase for it to be included in the calculation (traditionally 2, 1 also makes sense)
        # Preserve Thermochimica default unless set
        if minSpecies:
            inputFile.write(f'min species       = {minSpecies}\n')
        # Fuzzy stoichiometry settings
        inputFile.write(f'fuzzy             = {".TRUE." if fuzzyStoichiometry else ".FALSE."}\n')
        if (fuzzyStoichiometry and (fuzzyMagnitude >= 0)):
            inputFile.write(f'fuzzy magnitude   = {fuzzyMagnitude}\n')
        inputFile.write(f'gibbs min         = {".TRUE." if gibbsMinCheck else ".FALSE."}\n')

def RunRunCalculationList(filename,checkOutput=False,jsonName=None,thermochimica_path = '.'):
    thermoOut = None
    if checkOutput:
        thermoOut = subprocess.check_output([f'{thermochimica_path}/bin/RunCalculationList',filename]).decode("utf-8")
    else:
        subprocess.run([f'{thermochimica_path}/bin/RunCalculationList',filename])
    if jsonName:
        try:
            shutil.copy2(f'{thermochimica_path}/outputs/thermoout.json', f'{thermochimica_path}/outputs/{jsonName}')
        except:
            pass
    return thermoOut

def RunInputScript(filename,checkOutput=False,jsonName=None,thermochimica_path = '.'):
    thermoOut = None
    if checkOutput:
        thermoOut = subprocess.check_output([f'{thermochimica_path}/bin/InputScriptMode',filename]).decode("utf-8")
    else:
        subprocess.run([f'{thermochimica_path}/bin/InputScriptMode',filename])
    if jsonName:
        try:
            shutil.copy2(f'{thermochimica_path}/outputs/thermoout.json', f'{thermochimica_path}/outputs/{jsonName}')
        except:
            pass
    return thermoOut

def makePlot(datafile,xkey,yused,legend=None,yused2=None,legend2=None,plotColor='colorful',plotColor2='colorful',plotMarker='.-',plotMarker2='*--',xlog=False,ylog=False,ylog2=False,xinv=False,xinvScale=1,interactive=False,directory='outputs/',sort=True):
    # Prepend filename with directory
    datafile = f'{directory}{datafile}'
    # Do plot setup
    x,y,y2,xlab,ylab,ylab2 = plotDataSetup(datafile,xkey,yused,yused2=yused2,sort=sort)

    # If inverted x-axis requested, calculate new x values
    if xinv:
        x = [xinvScale/v for v in x]

    # Start figure
    fig = plt.figure()
    if interactive:
        plt.ion()
    lns=[]
    if yused2:
        ax = fig.add_axes([0.2, 0.1, 0.65, 0.85])
    else:
        ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(y))))
    for yi in range(len(y)):
        if plotColor == 'colorful':
            c = next(color)
        else:
            c = 'k'
        if legend:
            lns = lns + ax.plot(x,y[yi],plotMarker,c=c,label=legend[yi])
        else:
            ax.plot(x,y[yi],plotMarker,c=c)

    # Set x-axis label
    ax.set_xlabel(xlab)
    if xinv:
        ax.set_xlabel(f'{xinvScale}/{xlab}')
    
    if xlog:
        ax.set_xscale('log')
    ax.set_ylabel(ylab)
    if yused2:
        ax2 = ax.twinx()
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(y2))))
        for yi in range(len(y2)):
            if plotColor2 == 'colorful':
                c = next(color)
            else:
                c = 'k'
            if legend2:
                lns = lns + ax2.plot(x,y2[yi],plotMarker2,c=c,label=legend2[yi])
            else:
                ax2.plot(x,y2[yi],plotMarker2,c=c)
        ax2.set_ylabel(ylab2)
        if ylog2:
            ax2.set_yscale('log')
    labs = [l.get_label() for l in lns]
    if ylog:
        ax.set_yscale('log')
    if legend or legend2:
        ax.legend(lns, labs, loc=0)
    plt.show()
    plt.pause(0.001)

    return fig

def selectData(yen,ykey,leg,yen2=None,ykey2=None,leg2=None):
    # Select for left-hand y-axis
    yused = []
    legend = []
    for i in range(len(yen)):
        if yen[i]:
            yused.append(ykey[i])
            legend.append(leg[i])

    # Select for right-hand y-axis
    yused2 = []
    legend2 = []
    for i in range(len(yen2)):
        if yen2[i]:
            yused2.append(ykey2[i])
            legend2.append(leg2[i])

    return yused, legend, yused2, legend2

def plotDataSetup(datafile,xkey,yused,yused2=None,sort=True):
    # Init axes
    x = []
    y = [[] for _ in range(len(yused))]
    if yused2:
        y2 = [[] for _ in range(len(yused2))]
    else:
        y2 = []
    
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
            if yused2:
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
            else:
                x.append(data[j][xkey])
        except:
            # do nothing
            continue
    
    # Sort data unless asked not to
    if sort:
        sortOrder = [i for i,_ in sorted(enumerate(x),key=lambda s:s[1])]
        x = [x[i] for i in sortOrder]
        y = [[y[yi][i] for i in sortOrder] for yi in range(len(yused))]
        if yused2:
            y2 = [[y2[yi][i] for i in sortOrder] for yi in range(len(yused2))]


    # Set axis labels based on tags
    if xkey == 'iteration':
        xlab = 'Iteration'
    elif xkey == 'temperature':
        xlab = 'Temperature [K]'
    elif xkey == 'pressure':
        xlab = 'Pressure [atm]'
    else:
        xlab = ''

    # Use first key, last entry
    if yused:
        ykey = yused[0][-1]
        ylab = getYlab(ykey)
    else:
        ylab = None

    if yused2:
        ykey = yused2[0][-1]
        ylab2 = getYlab(ykey)
    else:
        ylab2 = None

    return x,y,y2,xlab,ylab,ylab2

def getYlab(ykey):
    if ykey == 'temperature':
        ylab = 'Temperature [K]'
    elif ykey == 'pressure':
        ylab = 'Pressure [atm]'
    elif ykey == 'integral Gibbs energy':
        ylab = 'Integral Gibbs Energy [J]'
    elif ykey == 'functional norm':
        ylab = 'Functional Norm'
    elif ykey == 'GEM iterations':
        ylab = 'GEM Iterations'
    elif ykey == 'heat capacity':
        ylab = 'Heat Capacity'
    elif ykey == 'enthalpy':
        ylab = 'Enthalpy'
    elif ykey == 'entropy':
        ylab = 'Entropy'
    elif ykey == '# phases':
        ylab = 'Number of Stable Phases'
    elif ykey == '# solution phases':
        ylab = 'Number of Stable Phases'
    elif ykey == '# pure condensed phases':
        ylab = 'Number of Stable Phases'
    elif ykey == 'moles':
        ylab = 'Moles'
    elif ykey == 'chemical potential':
        ylab = 'Chemical Potential [J]'
    elif ykey == 'driving force':
        ylab = 'Driving Force [J]'
    elif ykey == 'moles of element in phase':
        ylab = 'Moles of Element in Phase'
    elif ykey == 'mole fraction of phase by element':
        ylab = 'Mole Fraction of Phase by Element'
    elif ykey == 'mole fraction of element by phase':
        ylab = 'Mole Fraction of Element by Phase'
    elif ykey == 'mole fraction':
        ylab = 'Mole Fraction'
    elif ykey == 'mole fraction of endmembers':
        ylab = 'Mole Fraction'
    elif ykey == 'vapor pressure':
        ylab = 'Vapor Pressure [atm]'
    elif ykey == 'moles of elements':
        ylab = 'Moles'
    elif ykey == 'element potential':
        ylab = 'Element Potential [J]'
    else:
        ylab = ''
    
    return ylab    

def exportPlotScript(filename,datafile,xkey,yused,legend=None,yused2=None,legend2=None,plotColor='colorful',plotColor2='colorful',plotMarker='.-',plotMarker2='*--',xlog=False,ylog=False,ylog2=False):
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
        f.write(f'yused    = {yused}\n')
        f.write(f'legend   = {legend}\n')
        f.write(f'yused2   = {yused2}\n')
        f.write(f'legend2  = {legend2}\n')
        # Call plotDataSetup (enables use of datafile path)
        f.write('x,y,y2,xlab,ylab,ylab2 = thermoTools.plotDataSetup(datafile,xkey,yused,yused2=yused2)\n')
        # Figure generation commands
        f.write('lns=[]\n')
        f.write('# Start figure\n')
        f.write('fig = plt.figure()\n')
        if yused2:
            f.write('ax  = fig.add_axes([0.2, 0.1, 0.65, 0.85])\n')
        else:
            f.write('ax  = fig.add_axes([0.2, 0.1, 0.75, 0.85])\n')
        f.write('color = iter(plt.cm.rainbow(np.linspace(0, 1, len(y))))\n')
        f.write('for yi in range(len(y)):\n')
        if plotColor == 'colorful':
            f.write('    c = next(color)\n')
        else:
            f.write('    c = \'k\'\n')
        if legend:
            f.write(f'    lns = lns + ax.plot(x,y[yi],\'{plotMarker}\',c=c,label=legend[yi])\n')
        else:
            f.write(f'    ax.plot(x,y[yi],\'{plotMarker}\',c=c)\n')
        if yused2:
            f.write('ax2 = ax.twinx()\n')
            f.write('color = iter(plt.cm.rainbow(np.linspace(0, 1, len(y2))))\n')
            f.write('for yi in range(len(y2)):\n')
            if plotColor2 == 'colorful':
                f.write('    c = next(color)\n')
            else:
                f.write('    c = \'k\'\n')
            if legend2:
                f.write(f'    lns = lns + ax2.plot(x,y2[yi],\'{plotMarker2}\',c=c,label=legend2[yi])\n')
            else:
                f.write(f'    ax2.plot(x,y2[yi],\'{plotMarker2}\',c=c)\n')
            f.write('ax2.set_ylabel(ylab2)\n')
            if ylog2:
                f.write("ax2.set_yscale('log')\n")
        f.write('labs = [l.get_label() for l in lns]\n')
        if legend or legend2:
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

