import PySimpleGUI as sg
import subprocess
import math
import os

sg.theme('GreenTan')

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

menu_def = [['&File', ['&Open', '&Save', '&Exit', 'Properties']],
            ['&Edit', ['Paste', ['Special', 'Normal', ], 'Undo'], ],
            ['&Help', '&About...'], ]

database_frame = sg.Frame(layout=[
[sg.In(size=(30, 1), enable_events=True, key="-FOLDER-"), sg.FolderBrowse('Browse',size=(8,1))],
[sg.Listbox(values=[], enable_events=True, size=(40, 20), key="-FILE LIST-")]], title='Thermodynamic Database', relief=sg.RELIEF_SUNKEN)

system_input_frame = sg.Frame(layout=[
[sg.Text('Temperature')],
[sg.Input(key='-temperature-',size=(16,1)), sg.Combo(['K', 'C', 'F'], default_value='K',key='_tunit_')],
[sg.Text('Pressure')],
[sg.Input(key='-pressure-',size=(16,1)), sg.Combo(['atm', 'Pa', 'bar'], default_value='atm',key='_punit_')]], title='System Inputs', relief=sg.RELIEF_SUNKEN)

layout = [
    [sg.Menu(menu_def, tearoff=True)],
    [sg.Text('THERMOCHIMICA', size=(30,1), justification='center',font=("Helvetica", 25), relief=sg.RELIEF_RIDGE)],
    [database_frame, system_input_frame], ]

window = sg.Window('Thermochimica GUI', layout, default_element_size=(40, 1), grab_anywhere=False)
event, values = window.read()
window.close()

sg.Popup('Results',
         'The results of the calculation shown here.')
