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