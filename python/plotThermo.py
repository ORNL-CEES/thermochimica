import PySimpleGUI as sg
import subprocess
import math
import os
import json
import matplotlib.pyplot as plt
import numpy as np

timeout = 50
inputSize = 20
buttonSize = 12
windowList = []
popupLocation = [300,0]

class DataWindow:
    def __init__(self):
        windowList.append(self)
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
        self.folder = os.getcwd()
        try:
            file_list = os.listdir(self.folder)
        except:
            file_list = []
        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(self.folder, f))
            and f.lower().endswith(('.json', '.JSON'))
        ]
        self.sgw = sg.Window('Thermochimica output data selection', file_list_column, location = [0,0], finalize=True)
        self.sgw["-FILE LIST-"].update(fnames)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event == "-FOLDER-":
            self.folder = values["-FOLDER-"]
            try:
                file_list = os.listdir(self.folder)
            except:
                file_list = []

            fnames = [
                f
                for f in file_list
                if os.path.isfile(os.path.join(self.folder, f))
                and f.lower().endswith(('.json', '.JSON'))
            ]
            self.sgw["-FILE LIST-"].update(fnames)
        elif event == "-FILE LIST-":  # A file was chosen from the listbox
            try:
                datafile = os.path.join(
                    self.folder, values["-FILE LIST-"][0]
                )
                pWindow = PlotWindow(datafile)
                self.children.append(pWindow)
            except:
                pass

class PlotWindow:
    def __init__(self,datafile):
        windowList.append(self)
        self.datafile = datafile
        self.ykey = [[]]
        self.yen = []
        self.leg = []
        self.ykey2 = [[]]
        self.yen2 = []
        self.leg2 = []
        self.ylab = ''
        self.ylab2 = ''
        f = open(self.datafile,)
        self.data = json.load(f)
        f.close()
        if list(self.data.keys())[0] != '1':
            print('Output does not contain data series')
            exit()
        optionsLayout = [
                          [sg.Text('x-axis')],[sg.Combo(['iteration', 'temperature', 'pressure'], default_value='iteration', key='-xaxis-')],[sg.Checkbox('Log scale',key='-xlog-')],
                          [sg.Text('y-axis')],[sg.Combo(['temperature', 'pressure', 'moles', 'mole fraction', 'chemical potential', 'vapor pressure',
                           'moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase',
                           'element potential', 'integral Gibbs energy', 'functional norm', '# phases'],
                            key='-yaxis-', enable_events=True)],[sg.Checkbox('Log scale',key='-ylog-')],
                          [sg.Text('y-axis 2')],[sg.Combo(['','temperature', 'pressure', 'moles', 'mole fraction', 'chemical potential', 'vapor pressure',
                           'moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase',
                           'element potential', 'integral Gibbs energy', 'functional norm', '# phases'],
                            key='-yaxis2-', enable_events=True, disabled=True)],[sg.Checkbox('Log scale',key='-y2log-')]
                        ]
        plotLayout = [optionsLayout,
                      [sg.Column([[sg.Button('Plot', disabled = True, size = buttonSize)],
                                  [sg.Button('Export Plot', disabled = True, size = buttonSize)]
                                 ],vertical_alignment='t'),
                      sg.Column([[sg.Button('Export Plot Script', disabled = True, size = buttonSize)],
                                  [sg.Button('Plot Settings', size = buttonSize)]
                                 ],vertical_alignment='t')]]
        self.sgw = sg.Window('Thermochimica plot setup', plotLayout, location = [400,0], finalize=True)
        self.children = []
        self.currentPlot = []
        self.figureList = []
        self.exportFormat = 'png'
        self.exportFileName = 'plot'
        self.plotMarker = '.-'
        self.plotMarker2 = '*--'
        self.plotColor = 'colorful'
        self.plotColor2 = 'colorful'
        self.exportDPI = 300
    def close(self):
        for child in self.children:
            child.close()
        for fig in self.figureList:
            plt.close(fig=fig)
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event == '-yaxis-':
            self.ykey = [[]]
            self.yen = []
            self.leg = []
            self.sgw.Element('Plot').Update(disabled = True)
            if values['-yaxis-'] in ['temperature','pressure','integral Gibbs energy','functional norm']:
                try:
                    self.ykey[0].append(values['-yaxis-'])
                    self.yen.append(True)
                    self.leg.append(values['-yaxis-'])
                    if values['-yaxis-'] == 'temperature':
                        self.ylab = 'Temperature [K]'
                    elif values['-yaxis-'] == 'pressure':
                        self.ylab = 'Pressure [atm]'
                    elif values['-yaxis-'] == 'integral Gibbs energy':
                        self.ylab = 'Integral Gibbs Energy [J]'
                    elif values['-yaxis-'] == 'functional norm':
                        self.ylab = 'Functional Norm'
                    self.sgw.Element('Plot').Update(disabled = False)
                    self.sgw.Element('-yaxis2-').Update(disabled = False)
                except:
                    return
            elif values['-yaxis-'] == '# phases':
                try:
                    self.ykey[0].append('# solution phases')
                    self.yen.append(True)
                    self.leg.append('# of Solution Phases')
                    self.ykey.append([])
                    self.ykey[1].append('# pure condensed phases')
                    self.yen.append(True)
                    self.leg.append('# of Pure Condensed Phases')
                    self.ylab = 'Number of Stable Phases'
                    self.sgw.Element('Plot').Update(disabled = False)
                    self.sgw.Element('-yaxis2-').Update(disabled = False)
                except:
                    return
            elif values['-yaxis-'] in ['moles','chemical potential']:
                self.ykey = []
                try:
                    solutionPhases = list(self.data['1']['solution phases'].keys())
                    pureCondensedPhases = list(self.data['1']['pure condensed phases'].keys())
                except:
                    return
                phaseColumns = []
                yi = 0
                for j in solutionPhases:
                    phaseColumns.append([[sg.Text(j)]])
                    if values['-yaxis-'] == 'moles':
                        # total moles of solution phase
                        try:
                            self.ykey.append(['solution phases',j,values['-yaxis-']])
                            self.yen.append(False)
                            phaseColumns[-1].append([sg.Checkbox(self.ykey[yi][-2],key=str(yi))])
                            self.leg.append(j)
                            self.ylab = 'Moles'
                            yi = yi + 1
                        except:
                            continue
                    else:
                        self.ylab = 'Chemical Potential [J]'
                    if self.data['1']['solution phases'][j]['phase model'] in ['SUBG', 'SUBQ']:
                        speciesLabel = 'quadruplets'
                    else:
                        speciesLabel = 'species'
                    for k in list(self.data['1']['solution phases'][j][speciesLabel].keys()):
                        try:
                            self.ykey.append(['solution phases',j,speciesLabel,k,values['-yaxis-']])
                            self.yen.append(False)
                            phaseColumns[-1].append([sg.Checkbox(self.ykey[yi][-2],key=str(yi))])
                            self.leg.append(j+': '+k)
                            yi = yi + 1
                        except:
                            continue
                phaseColumns.append([[sg.Text('Pure Condensed Phases')]])
                for j in pureCondensedPhases:
                    try:
                        self.ykey.append(['pure condensed phases',j,values['-yaxis-']])
                        self.yen.append(False)
                        phaseColumns[-1].append([sg.Checkbox(self.ykey[yi][-2],key=str(yi))])
                        self.leg.append(j)
                        yi = yi + 1
                    except:
                        continue
                phaseSelectLayout = [[]]
                for j in phaseColumns:
                    phaseSelectLayout[0].append(sg.Column(j,vertical_alignment='t'))
                phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
                selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, location = popupLocation, finalize=True)
                while True:
                    event, values = selectWindow.read()
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Accept':
                        for yi in range(len(self.ykey)):
                            self.yen[yi] = values[str(yi)]
                        self.sgw.Element('Plot').Update(disabled = False)
                        self.sgw.Element('-yaxis2-').Update(disabled = False)
                        break
                selectWindow.close()
            elif values['-yaxis-'] in ['moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase']:
                self.ykey = []
                solutionPhases = list(self.data['1']['solution phases'].keys())
                pureCondensedPhases = list(self.data['1']['pure condensed phases'].keys())
                phaseColumns = []
                yi = 0
                for j in solutionPhases:
                    phaseColumns.append([[sg.Text(j)]])
                    if values['-yaxis-'] == 'moles of element in phase':
                        self.ylab = 'Moles of Element in Phase'
                    if values['-yaxis-'] == 'mole fraction of phase by element':
                        self.ylab = 'Mole Fraction of Phase by Element'
                    else:
                        self.ylab = 'Mole Fraction of Element by Phase'
                    for k in list(self.data['1']['solution phases'][j]['elements'].keys()):
                        try:
                            self.ykey.append(['solution phases',j,'elements',k,values['-yaxis-']])
                            self.yen.append(False)
                            phaseColumns[-1].append([sg.Checkbox(self.ykey[yi][-2],key=str(yi))])
                            self.leg.append(j+': '+k)
                            yi = yi + 1
                        except:
                            continue
                for j in pureCondensedPhases:
                    phaseColumns.append([[sg.Text(j)]])
                    for k in list(self.data['1']['pure condensed phases'][j]['elements'].keys()):
                        try:
                            self.ykey.append(['pure condensed phases',j,'elements',k,values['-yaxis-']])
                            self.yen.append(False)
                            phaseColumns[-1].append([sg.Checkbox(self.ykey[yi][-2],key=str(yi))])
                            self.leg.append(j+': '+k)
                            yi = yi + 1
                        except:
                            continue
                phaseSelectLayout = [[]]
                for j in phaseColumns:
                    phaseSelectLayout[0].append(sg.Column(j,vertical_alignment='t'))
                phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
                selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, location = popupLocation, finalize=True)
                while True:
                    event, values = selectWindow.read()
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Accept':
                        for yi in range(len(self.ykey)):
                            self.yen[yi] = values[str(yi)]
                        self.sgw.Element('Plot').Update(disabled = False)
                        self.sgw.Element('-yaxis2-').Update(disabled = False)
                        break
                selectWindow.close()
            elif values['-yaxis-'] == 'mole fraction':
                self.ykey = []
                self.ylab = 'Mole Fraction'
                solutionPhases = list(self.data['1']['solution phases'].keys())
                phaseColumns = []
                yi = 0
                for j in solutionPhases:
                    phaseColumns.append([[sg.Text(j)]])
                    if self.data['1']['solution phases'][j]['phase model'] in ['SUBG', 'SUBQ']:
                        speciesLabel = 'quadruplets'
                    else:
                        speciesLabel = 'species'
                    for k in list(self.data['1']['solution phases'][j][speciesLabel].keys()):
                        try:
                            self.ykey.append(['solution phases',j,speciesLabel,k,values['-yaxis-']])
                            self.yen.append(False)
                            phaseColumns[-1].append([sg.Checkbox(self.ykey[yi][-2],key=str(yi))])
                            self.leg.append(j+': '+k)
                            yi = yi + 1
                        except:
                            continue
                phaseSelectLayout = [[]]
                for j in phaseColumns:
                    phaseSelectLayout[0].append(sg.Column(j,vertical_alignment='t'))
                phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
                selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, location = popupLocation, finalize=True)
                while True:
                    event, values = selectWindow.read()
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Accept':
                        for yi in range(len(self.ykey)):
                            self.yen[yi] = values[str(yi)]
                        self.sgw.Element('Plot').Update(disabled = False)
                        self.sgw.Element('-yaxis2-').Update(disabled = False)
                        break
                selectWindow.close()
            elif values['-yaxis-'] == 'vapor pressure':
                self.ykey = []
                self.ylab = 'Vapor Pressure [atm]'
                solutionPhases = list(self.data['1']['solution phases'].keys())
                phaseColumns = []
                yi = 0
                for j in solutionPhases:
                    if self.data['1']['solution phases'][j]['phase model'] != 'IDMX':
                        break
                    for k in list(self.data['1']['solution phases'][j]['species'].keys()):
                        try:
                            self.ykey.append(['solution phases',j,'species',k,values['-yaxis-']])
                            self.yen.append(False)
                            phaseColumns.append([sg.Checkbox(self.ykey[yi][-2],key=str(yi))])
                            self.leg.append(j+': '+k)
                            yi = yi + 1
                        except:
                            continue
                    break
                phaseSelectLayout = phaseColumns
                phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
                selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, location = popupLocation, finalize=True)
                while True:
                    event, values = selectWindow.read()
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Accept':
                        for yi in range(len(self.ykey)):
                            self.yen[yi] = values[str(yi)]
                        self.sgw.Element('Plot').Update(disabled = False)
                        self.sgw.Element('-yaxis2-').Update(disabled = False)
                        break
                selectWindow.close()
            elif values['-yaxis-'] == 'element potential':
                self.ykey = []
                self.ylab = 'Element Potential [J]'
                elements = list(self.data['1']['elements'].keys())
                for j in elements:
                    try:
                        self.ykey.append(['elements',j,values['-yaxis-']])
                        self.yen.append(True)
                        self.leg.append(j)
                    except:
                        continue
                self.sgw.Element('Plot').Update(disabled = False)
                self.sgw.Element('-yaxis2-').Update(disabled = False)
        elif event == '-yaxis2-':
            self.ykey2 = [[]]
            self.yen2 = []
            self.leg2 = []
            if values['-yaxis2-'] in ['temperature','pressure','integral Gibbs energy','functional norm']:
                self.ykey2[0].append(values['-yaxis2-'])
                self.yen2.append(True)
                self.leg2.append(values['-yaxis2-'])
                if values['-yaxis2-'] == 'temperature':
                    self.ylab2 = 'Temperature [K]'
                elif values['-yaxis2-'] == 'pressure':
                    self.ylab2 = 'Pressure [atm]'
                elif values['-yaxis2-'] == 'integral Gibbs energy':
                    self.ylab2 = 'Integral Gibbs Energy [J]'
                elif values['-yaxis2-'] == 'functional norm':
                    self.ylab2 = 'Functional Norm'
                self.sgw.Element('Plot').Update(disabled = False)
            elif values['-yaxis2-'] == '# phases':
                self.ykey2[0].append('# solution phases')
                self.yen2.append(True)
                self.leg2.append('# of Solution Phases')
                self.ykey2.append([])
                self.ykey2[1].append('# pure condensed phases')
                self.yen2.append(True)
                self.leg2.append('# of Pure Condensed Phases')
                self.ylab2 = 'Number of Stable Phases'
                self.sgw.Element('Plot').Update(disabled = False)
            elif values['-yaxis2-'] in ['moles','chemical potential']:
                self.ykey2 = []
                solutionPhases = list(self.data['1']['solution phases'].keys())
                pureCondensedPhases = list(self.data['1']['pure condensed phases'].keys())
                phaseColumns = []
                yi = 0
                for j in solutionPhases:
                    phaseColumns.append([[sg.Text(j)]])
                    if values['-yaxis2-'] == 'moles':
                        # total moles of solution phase
                        self.ykey2.append(['solution phases',j,values['-yaxis2-']])
                        self.yen2.append(False)
                        phaseColumns[-1].append([sg.Checkbox(self.ykey2[yi][-2],key=str(yi))])
                        self.leg2.append(j)
                        self.ylab2 = 'Moles'
                        yi = yi + 1
                    else:
                        self.ylab2 = 'Chemical Potential [J]'
                    if self.data['1']['solution phases'][j]['phase model'] in ['SUBG', 'SUBQ']:
                        speciesLabel = 'quadruplets'
                    else:
                        speciesLabel = 'species'
                    for k in list(self.data['1']['solution phases'][j][speciesLabel].keys()):
                        self.ykey2.append(['solution phases',j,speciesLabel,k,values['-yaxis2-']])
                        self.yen2.append(False)
                        phaseColumns[-1].append([sg.Checkbox(self.ykey2[yi][-2],key=str(yi))])
                        self.leg2.append(j+': '+k)
                        yi = yi + 1
                phaseColumns.append([[sg.Text('Pure Condensed Phases')]])
                for j in pureCondensedPhases:
                    self.ykey2.append(['pure condensed phases',j,values['-yaxis2-']])
                    self.yen2.append(False)
                    phaseColumns[-1].append([sg.Checkbox(self.ykey2[yi][-2],key=str(yi))])
                    self.leg2.append(j)
                    yi = yi + 1
                phaseSelectLayout = [[]]
                for j in phaseColumns:
                    phaseSelectLayout[0].append(sg.Column(j,vertical_alignment='t'))
                phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
                selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, location = popupLocation, finalize=True)
                while True:
                    event, values = selectWindow.read()
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Accept':
                        for yi in range(len(self.ykey2)):
                            self.yen2[yi] = values[str(yi)]
                        self.sgw.Element('Plot').Update(disabled = False)
                        break
                selectWindow.close()
            elif values['-yaxis2-'] in ['moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase']:
                self.ykey2 = []
                solutionPhases = list(self.data['1']['solution phases'].keys())
                pureCondensedPhases = list(self.data['1']['pure condensed phases'].keys())
                phaseColumns = []
                yi = 0
                for j in solutionPhases:
                    phaseColumns.append([[sg.Text(j)]])
                    if values['-yaxis2-'] == 'moles of element in phase':
                        self.ylab2 = 'Moles of Element in Phase'
                    if values['-yaxis2-'] == 'mole fraction of phase by element':
                        self.ylab2 = 'Mole Fraction of Phase by Element'
                    else:
                        self.ylab2 = 'Mole Fraction of Element by Phase'
                    for k in list(self.data['1']['solution phases'][j]['elements'].keys()):
                        self.ykey2.append(['solution phases',j,'elements',k,values['-yaxis2-']])
                        self.yen2.append(False)
                        phaseColumns[-1].append([sg.Checkbox(self.ykey2[yi][-2],key=str(yi))])
                        self.leg2.append(j+': '+k)
                        yi = yi + 1
                for j in pureCondensedPhases:
                    phaseColumns.append([[sg.Text(j)]])
                    for k in list(self.data['1']['pure condensed phases'][j]['elements'].keys()):
                        self.ykey2.append(['pure condensed phases',j,'elements',k,values['-yaxis2-']])
                        self.yen2.append(False)
                        phaseColumns[-1].append([sg.Checkbox(self.ykey2[yi][-2],key=str(yi))])
                        self.leg2.append(j+': '+k)
                        yi = yi + 1
                phaseSelectLayout = [[]]
                for j in phaseColumns:
                    phaseSelectLayout[0].append(sg.Column(j,vertical_alignment='t'))
                phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
                selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, location = popupLocation, finalize=True)
                while True:
                    event, values = selectWindow.read()
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Accept':
                        for yi in range(len(self.ykey2)):
                            self.yen2[yi] = values[str(yi)]
                        self.sgw.Element('Plot').Update(disabled = False)
                        break
                selectWindow.close()
            elif values['-yaxis2-'] == 'mole fraction':
                self.ykey2 = []
                self.ylab2 = 'Mole Fraction'
                solutionPhases = list(self.data['1']['solution phases'].keys())
                phaseColumns = []
                yi = 0
                for j in solutionPhases:
                    phaseColumns.append([[sg.Text(j)]])
                    if self.data['1']['solution phases'][j]['phase model'] in ['SUBG', 'SUBQ']:
                        speciesLabel = 'quadruplets'
                    else:
                        speciesLabel = 'species'
                    for k in list(self.data['1']['solution phases'][j][speciesLabel].keys()):
                        self.ykey2.append(['solution phases',j,speciesLabel,k,values['-yaxis2-']])
                        self.yen2.append(False)
                        phaseColumns[-1].append([sg.Checkbox(self.ykey2[yi][-2],key=str(yi))])
                        self.leg2.append(j+': '+k)
                        yi = yi + 1
                phaseSelectLayout = [[]]
                for j in phaseColumns:
                    phaseSelectLayout[0].append(sg.Column(j,vertical_alignment='t'))
                phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
                selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, location = popupLocation, finalize=True)
                while True:
                    event, values = selectWindow.read()
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Accept':
                        for yi in range(len(self.ykey2)):
                            self.yen2[yi] = values[str(yi)]
                        self.sgw.Element('Plot').Update(disabled = False)
                        break
                selectWindow.close()
                self.sgw.Element('Plot').Update(disabled = False)
            elif values['-yaxis2-'] == 'vapor pressure':
                self.ykey2 = []
                self.ylab2 = 'Vapor Pressure [atm]'
                solutionPhases = list(self.data['1']['solution phases'].keys())
                phaseColumns = []
                yi = 0
                for j in solutionPhases:
                    if self.data['1']['solution phases'][j]['phase model'] != 'IDMX':
                        break
                    for k in list(self.data['1']['solution phases'][j]['species'].keys()):
                        try:
                            self.ykey2.append(['solution phases',j,'species',k,values['-yaxis2-']])
                            self.yen2.append(False)
                            phaseColumns.append([sg.Checkbox(self.ykey2[yi][-2],key=str(yi))])
                            self.leg2.append(j+': '+k)
                            yi = yi + 1
                        except:
                            continue
                    break
                phaseSelectLayout = phaseColumns
                phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
                selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, location = popupLocation, finalize=True)
                while True:
                    event, values = selectWindow.read()
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Accept':
                        for yi in range(len(self.ykey2)):
                            self.yen2[yi] = values[str(yi)]
                        self.sgw.Element('Plot').Update(disabled = False)
                        self.sgw.Element('-yaxis2-').Update(disabled = False)
                        break
                selectWindow.close()
            elif values['-yaxis2-'] == 'element potential':
                self.ykey2 = []
                self.ylab2 = 'Element Potential [J]'
                elements = list(self.data['1']['elements'].keys())
                for j in elements:
                    self.ykey2.append(['elements',j,values['-yaxis2-']])
                    self.yen2.append(True)
                    self.leg2.append(j)
                self.sgw.Element('Plot').Update(disabled = False)
        elif event == 'Plot':
            x = []
            y = []
            y2 = []
            xkey = values['-xaxis-']
            for yi in range(len(self.ykey)):
                y.append([])
            for yi in range(len(self.ykey2)):
                y2.append([])
            for j in self.data.keys():
                try:
                    for yi in range(len(self.ykey)):
                        if len(self.ykey[yi]) == 1:
                            y[yi].append(self.data[j][self.ykey[yi][0]])
                        elif len(self.ykey[yi]) == 3:
                            y[yi].append(self.data[j][self.ykey[yi][0]][self.ykey[yi][1]][self.ykey[yi][2]])
                        elif len(self.ykey[yi]) == 5:
                            if self.ykey[yi][4] == 'vapor pressure':
                                y[yi].append(self.data[j][self.ykey[yi][0]][self.ykey[yi][1]][self.ykey[yi][2]][self.ykey[yi][3]]['mole fraction']*self.data[j]['pressure'])
                            else:
                                y[yi].append(self.data[j][self.ykey[yi][0]][self.ykey[yi][1]][self.ykey[yi][2]][self.ykey[yi][3]][self.ykey[yi][4]])
                    for yi in range(len(self.ykey2)):
                        if len(self.ykey2[yi]) == 1:
                            y2[yi].append(self.data[j][self.ykey2[yi][0]])
                        elif len(self.ykey2[yi]) == 3:
                            y2[yi].append(self.data[j][self.ykey2[yi][0]][self.ykey2[yi][1]][self.ykey2[yi][2]])
                        elif len(self.ykey2[yi]) == 5:
                            if self.ykey2[yi][4] == 'vapor pressure':
                                y2[yi].append(self.data[j][self.ykey2[yi][0]][self.ykey2[yi][1]][self.ykey2[yi][2]][self.ykey2[yi][3]]['mole fraction']*self.data[j]['pressure'])
                            else:
                                y2[yi].append(self.data[j][self.ykey2[yi][0]][self.ykey2[yi][1]][self.ykey2[yi][2]][self.ykey2[yi][3]][self.ykey2[yi][4]])
                    if xkey == 'iteration':
                        x.append(int(j))
                        xlab = 'Iteration'
                    else:
                        x.append(self.data[j][xkey])
                        if xkey == 'temperature':
                            xlab = 'Temperature [K]'
                        elif xkey == 'pressure':
                            xlab = 'Pressure [atm]'
                except:
                    # do nothing
                    continue
            # Start figure
            fig = plt.figure()
            plt.ion()
            lns=[]
            if True in self.yen2:
                ax = fig.add_axes([0.2, 0.1, 0.65, 0.85])
            else:
                ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(self.yen))))
            for yi in range(len(self.yen)):
                if self.yen[yi]:
                    if self.plotColor == 'colorful':
                        c = next(color)
                    else:
                        c = 'k'
                    lns = lns + ax.plot(x,y[yi],self.plotMarker,c=c,label = self.leg[yi])
            ax.set_xlabel(xlab)
            if values['-xlog-']:
                ax.set_xscale('log')
            ax.set_ylabel(self.ylab)
            if True in self.yen2:
                ax2 = ax.twinx()
                color = iter(plt.cm.rainbow(np.linspace(0, 1, len(self.yen2))))
                for yi in range(len(self.yen2)):
                    if self.yen2[yi]:
                        if self.plotColor2 == 'colorful':
                            c = next(color)
                        else:
                            c = 'k'
                        lns = lns + ax2.plot(x,y2[yi],self.plotMarker2,c=c,label = self.leg2[yi])
                ax2.set_ylabel(self.ylab2)
                if values['-y2log-']:
                    ax2.set_yscale('log')
            labs = [l.get_label() for l in lns]
            if values['-ylog-']:
                ax.set_yscale('log')
            ax.legend(lns, labs, loc=0)
            plt.show()
            plt.pause(0.001)
            self.currentPlot = fig
            self.figureList.append(fig)
            self.sgw.Element('Export Plot').Update(disabled = False)
            self.sgw.Element('Export Plot Script').Update(disabled = False)
        elif event == 'Export Plot Script':
            with open('python/generatedPlotScript.py', 'w') as f:
                f.write('# Thermochimica-generated plot script\n')
                f.write('import matplotlib.pyplot as plt\n')
                f.write('x = '+"{}\n".format(x))
                f.write('y = '+"{}\n".format(y))
                f.write('xlab = \''+xlab+'\'\n')
                f.write('ylab = \''+self.ylab+'\'\n')
                f.write('yen = '+"{}\n".format(self.yen))
                f.write('leg = '+"{}\n".format(self.leg))
                f.write('lns=[]\n')
                f.write('# Start figure\n')
                f.write('fig = plt.figure()\n')
                if True in self.yen2:
                    f.write('ax  = fig.add_axes([0.2, 0.1, 0.65, 0.85])\n')
                else:
                    f.write('ax  = fig.add_axes([0.2, 0.1, 0.75, 0.85])\n')
                f.write('for yi in range(len(yen)):\n')
                f.write('    if yen[yi]:\n')
                f.write('        lns = lns + ax.plot(x,y[yi],\'.-\',label = leg[yi])\n')
                if True in self.yen2:
                    f.write('y2 = '+"{}\n".format(y2))
                    f.write('ylab2 = \''+self.ylab2+'\'\n')
                    f.write('yen2 = '+"{}\n".format(self.yen2))
                    f.write('leg2 = '+"{}\n".format(self.leg2))
                    f.write('ax2 = ax.twinx()\n')
                    f.write('for yi in range(len(yen2)):\n')
                    f.write('    if yen2[yi]:\n')
                    f.write('        lns = lns + ax2.plot(x,y2[yi],\'^--\',label = leg2[yi])\n')
                    f.write('ax2.set_ylabel(ylab2)\n')
                f.write('labs = [l.get_label() for l in lns]\n')
                f.write('ax.legend(lns, labs, loc=0)\n')
                f.write('ax.set_xlabel(xlab)\n')
                f.write('ax.set_ylabel(ylab)\n')
                f.write('plt.show()\n')
        elif event == 'Export Plot':
            self.exportPlot()
        elif event == 'Plot Settings':
            if self.plotMarker == '-':
                line  = True
                point = False
                both  = False
            elif self.plotMarker == '.':
                line  = False
                point = True
                both  = False
            else:
                line  = False
                point = False
                both  = True
            if self.plotMarker2 == '--':
                line2  = True
                point2 = False
                both2  = False
            elif self.plotMarker2 == '*':
                line2  = False
                point2 = True
                both2  = False
            else:
                line2  = False
                point2 = False
                both2  = True
            if self.plotColor == 'colorful':
                colorful = True
                bland    = False
            else:
                colorful = False
                bland    = True
            if self.plotColor2 == 'colorful':
                colorful2 = True
                bland2    = False
            else:
                colorful2 = False
                bland2    = True
            settingsLayout = [[sg.Column([[sg.Text('Marker Style:')],
                                          [sg.Radio('Lines', 'mstyle', default=line,  enable_events=True, key='-mline-')],
                                          [sg.Radio('Points','mstyle', default=point, enable_events=True, key='-mpoint-')],
                                          [sg.Radio('Both',  'mstyle', default=both,  enable_events=True, key='-mboth-')]
                                         ],vertical_alignment='t'),
                               sg.Column([[sg.Text('Marker Style 2:')],
                                          [sg.Radio('Lines', 'mstyle2', default=line2,  enable_events=True, key='-mline2-')],
                                          [sg.Radio('Points','mstyle2', default=point2, enable_events=True, key='-mpoint2-')],
                                          [sg.Radio('Both',  'mstyle2', default=both2,  enable_events=True, key='-mboth2-')]
                                         ],vertical_alignment='t')],
                              [sg.Column([[sg.Text('Plot Colors:')],
                                          [sg.Radio('Colorful', 'mcolor', default=colorful, enable_events=True, key='-mcolorful-')],
                                          [sg.Radio('Black',    'mcolor', default=bland,    enable_events=True, key='-mbland-')]
                                         ],vertical_alignment='t'),
                               sg.Column([[sg.Text('Plot Colors 2:')],
                                          [sg.Radio('Colorful', 'mcolor2', default=colorful2, enable_events=True, key='-mcolorful2-')],
                                          [sg.Radio('Black',    'mcolor2', default=bland2,    enable_events=True, key='-mbland2-')]
                                         ],vertical_alignment='t')],
                             [sg.Text('Export Filename'),sg.Input(key='-filename-',size=(inputSize,1))],
                             [sg.Text('Export Format'),sg.Combo(['png', 'pdf', 'ps', 'eps', 'svg'],default_value='png',key='-format-')],
                             [sg.Text('Export DPI'),sg.Input(key='-dpi-',size=(inputSize,1))],
                             [sg.Button('Accept')]]
            settingsWindow = SettingsWindow(self, settingsLayout)
            self.children.append(settingsWindow)
    def exportPlot(self):
        try:
            self.currentPlot.savefig(f'{self.exportFileName}.{self.exportFormat}', format=self.exportFormat, dpi=self.exportDPI)
        except:
            errorLayout = [[sg.Text('The export failed, try changing plot settings.')],[sg.Button('Continue'), sg.Button('Cancel')]]
            errorWindow = sg.Window('Plot export failed', errorLayout, location = [400,0], finalize=True, keep_on_top = True)
            while True:
                event, values = errorWindow.read(timeout=timeout)
                if event == sg.WIN_CLOSED or event == 'Continue':
                    break
            errorWindow.close()

class SettingsWindow:
    def __init__(self, parent, windowLayout):
        self.parent = parent
        windowList.append(self)
        self.sgw = sg.Window('Plot Settings', windowLayout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED:
            self.close()
        elif event == '-mline-':
            self.parent.plotMarker = '-'
        elif event =='-mpoint-':
            self.parent.plotMarker = '.'
        elif event =='-mboth-':
            self.parent.plotMarker = '.-'
        elif event == '-mline2-':
            self.parent.plotMarker2 = '--'
        elif event =='-mpoint2-':
            self.parent.plotMarker2 = '*'
        elif event =='-mboth2-':
            self.parent.plotMarker2 = '*--'
        elif event =='-mcolorful-':
            self.parent.plotColor = 'colorful'
        elif event =='-mbland-':
            self.parent.plotColor = 'bland'
        elif event =='-mcolorful2-':
            self.parent.plotColor2 = 'colorful'
        elif event =='-mbland2-':
            self.parent.plotColor2 = 'bland'
        elif event =='Accept':
            try:
                self.parent.exportFileName = str(values['-filename-'])
            except:
                pass
            self.parent.exportFormat = values['-format-']
            try:
                tempDPI = int(values['-dpi-'])
                if tempDPI > 0 > 10000:
                    self.parent.exportDPI = int(values['-dpi-'])
            except:
                pass
            self.close()

dataWindow = DataWindow()
while len(windowList) > 0:
    for window in windowList:
        window.read()
