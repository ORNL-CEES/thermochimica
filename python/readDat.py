import math
datafile = 'data/FPDB12-19-2019-mod.dat'
with open(datafile) as f:
    f.readline()
    line = f.readline()
    nElements = int(line[1:5])
    nSoln = int(line[6:10])
    elements = ['']*nElements
    elIt = 0
    for i in range(math.ceil((nSoln+3)/15)-1):
        f.readline()
    for i in range(math.ceil(nElements/3)):
        els = f.readline()
        elements[elIt] = els[1:25].strip()
        elIt = elIt + 1
        if elIt == nElements:
            break
        elements[elIt] = els[26:50].strip()
        elIt = elIt + 1
        if elIt == nElements:
            break
        elements[elIt] = els[51:75].strip()
        elIt = elIt + 1
        if elIt == nElements:
            break
    print(nElements,nSoln)
    print(elements)

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

i = 0
while i < nElements:
    try:
        index = atomic_number_map.index(elements[i])+1
        print(elements[i],index)
        i = i + 1
    except ValueError:
        print(elements[i]+' not in list')
        elements.remove(elements[i])
        print(elements)
        nElements = nElements - 1
