import PySimpleGUI as sg
import subprocess
import math
import os
import json
import matplotlib.pyplot as plt

f = open('thermoout.json',)

data = json.load(f)
f.close()
if list(data.keys())[0] != '1':
    print('Output does not contain data series')
    exit()

optionsLayout = [
                  [sg.Text('x-axis')],[sg.Combo(['iteration', 'temperature', 'pressure'], default_value='iteration', key='-xaxis-')],
                  [sg.Text('y-axis')],[sg.Combo(['temperature', 'pressure', 'moles', 'mole fraction',
                   'chemical potential', 'element potential'], key='-yaxis-', enable_events=True)]
                ]
plotLayout = [optionsLayout,
              [sg.Button('Plot', disabled = True), sg.Exit()]]
plotWindow = sg.Window('Thermochimica plot setup', plotLayout, finalize=True)
while True:
    event, values = plotWindow.read()
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    elif event == '-yaxis-':
        ykey = [[]]
        yen = []
        plotWindow.Element('Plot').Update(disabled = True)
        if values['-yaxis-'] == 'temperature' or values['-yaxis-'] == 'pressure':
            ykey[0].append(values['-yaxis-'])
            yen.append(True)
            plotWindow.Element('Plot').Update(disabled = False)
        elif values['-yaxis-'] == 'moles' or values['-yaxis-'] == 'chemical potential':
            ykey = []
            solutionPhases = list(data['1']['solution phases'].keys())
            pureCondensedPhases = list(data['1']['pure condensed phases'].keys())
            phaseColumns = []
            yi = 0
            for j in solutionPhases:
                phaseColumns.append([[sg.Text(j)]])
                if data['1']['solution phases'][j]['phase model'] == 'SUBG' or data['1']['solution phases'][j]['phase model'] == 'SUBQ':
                    speciesLabel = 'quadruplets'
                else:
                    speciesLabel = 'species'
                for k in list(data['1']['solution phases'][j][speciesLabel].keys()):
                    ykey.append(['solution phases',j,speciesLabel,k,values['-yaxis-']])
                    yen.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey[yi][-2],key=str(yi))])
                    yi = yi + 1
            phaseColumns.append([[sg.Text('Pure Condensed Phases')]])
            for j in pureCondensedPhases:
                ykey.append(['pure condensed phases',j,values['-yaxis-']])
                yen.append(False)
                phaseColumns[-1].append([sg.Checkbox(ykey[yi][-2],key=str(yi))])
                yi = yi + 1
            phaseSelectLayout = [[]]
            for j in phaseColumns:
                phaseSelectLayout[0].append(sg.Column(j,vertical_alignment='t'))
            phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
            selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, finalize=True)
            while True:
                event, values = selectWindow.read()
                if event == sg.WIN_CLOSED or event == 'Cancel':
                    break
                elif event == 'Accept':
                    for yi in range(len(ykey)):
                        yen[yi] = values[str(yi)]
                    plotWindow.Element('Plot').Update(disabled = False)
                    break
            selectWindow.close()
        elif values['-yaxis-'] == 'mole fraction':
            ykey = []
            solutionPhases = list(data['1']['solution phases'].keys())
            phaseColumns = []
            yi = 0
            for j in solutionPhases:
                phaseColumns.append([[sg.Text(j)]])
                if data['1']['solution phases'][j]['phase model'] == 'SUBG' or data['1']['solution phases'][j]['phase model'] == 'SUBQ':
                    speciesLabel = 'quadruplets'
                else:
                    speciesLabel = 'species'
                for k in list(data['1']['solution phases'][j][speciesLabel].keys()):
                    ykey.append(['solution phases',j,speciesLabel,k,values['-yaxis-']])
                    yen.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey[yi][-2],key=str(yi))])
                    yi = yi + 1
                phaseSelectLayout = [[]]
                for j in phaseColumns:
                    phaseSelectLayout[0].append(sg.Column(j,vertical_alignment='t'))
                phaseSelectLayout.append([sg.Button('Accept'), sg.Button('Cancel')])
                selectWindow = sg.Window('Thermochimica species selection', phaseSelectLayout, finalize=True)
                while True:
                    event, values = selectWindow.read()
                    if event == sg.WIN_CLOSED or event == 'Cancel':
                        break
                    elif event == 'Accept':
                        for yi in range(len(ykey)):
                            yen[yi] = values[str(yi)]
                        plotWindow.Element('Plot').Update(disabled = False)
                        break
                selectWindow.close()
            plotWindow.Element('Plot').Update(disabled = False)
        elif values['-yaxis-'] == 'element potential':
            ykey = []
            elements = list(data['1']['elements'].keys())
            for j in elements:
                ykey.append(['elements',j,values['-yaxis-']])
                yen.append(True)
            plotWindow.Element('Plot').Update(disabled = False)
    elif event == 'Plot':
        x = []
        y = []
        xkey = values['-xaxis-']
        for yi in range(len(ykey)):
            y.append([])
        for j in data.keys():
            if xkey == 'iteration':
                x.append(int(j))
            else:
                x.append(data[j][xkey])
            for yi in range(len(ykey)):
                if len(ykey[yi]) == 1:
                    y[yi].append(data[j][ykey[yi][0]])
                elif len(ykey[yi]) == 3:
                    y[yi].append(data[j][ykey[yi][0]][ykey[yi][1]][ykey[yi][2]])
                elif len(ykey[yi]) == 5:
                    y[yi].append(data[j][ykey[yi][0]][ykey[yi][1]][ykey[yi][2]][ykey[yi][3]][ykey[yi][4]])
        # Start figure
        fig = plt.figure()
        ax  = fig.add_axes([0.12, 0.1, 0.85, 0.85])
        for yi in range(len(ykey)):
            if yen[yi]:
                ax.plot(x,y[yi],'.-')
        plt.show()
plotWindow.close()
