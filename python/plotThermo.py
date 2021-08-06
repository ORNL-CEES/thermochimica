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
                   'chemical potential', 'element potential', 'integral Gibbs energy', 'functional norm', '# phases'],
                    key='-yaxis-', enable_events=True)]
                ]
plotLayout = [optionsLayout,
              [sg.Button('Plot', disabled = True), sg.Button('Export Plot Script', disabled = True)]]
plotWindow = sg.Window('Thermochimica plot setup', plotLayout, finalize=True)
while True:
    event, values = plotWindow.read()
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    elif event == '-yaxis-':
        ykey = [[]]
        yen = []
        leg = []
        plotWindow.Element('Plot').Update(disabled = True)
        if values['-yaxis-'] in ['temperature','pressure','integral Gibbs energy','functional norm']:
            ykey[0].append(values['-yaxis-'])
            yen.append(True)
            leg.append(values['-yaxis-'])
            if values['-yaxis-'] == 'temperature':
                ylab = 'Temperature [K]'
            elif values['-yaxis-'] == 'pressure':
                ylab = 'Pressure [atm]'
            elif values['-yaxis-'] == 'integral Gibbs energy':
                ylab = 'Integral Gibbs Energy [J]'
            elif values['-yaxis-'] == 'functional norm':
                ylab = 'Functional Norm'
            plotWindow.Element('Plot').Update(disabled = False)
        elif values['-yaxis-'] == '# phases':
            ykey[0].append('# solution phases')
            yen.append(True)
            leg.append('# of Solution Phases')
            ykey.append([])
            ykey[1].append('# pure condensed phases')
            yen.append(True)
            leg.append('# of Pure Condensed Phases')
            ylab = 'Number of Stable Phases'
            plotWindow.Element('Plot').Update(disabled = False)
        elif values['-yaxis-'] in ['moles','chemical potential']:
            ykey = []
            solutionPhases = list(data['1']['solution phases'].keys())
            pureCondensedPhases = list(data['1']['pure condensed phases'].keys())
            phaseColumns = []
            yi = 0
            for j in solutionPhases:
                phaseColumns.append([[sg.Text(j)]])
                if values['-yaxis-'] == 'moles':
                    # total moles of solution phase
                    ykey.append(['solution phases',j,values['-yaxis-']])
                    yen.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey[yi][-2],key=str(yi))])
                    leg.append(j)
                    ylab = 'Moles'
                    yi = yi + 1
                else:
                    ylab = 'Chemical Potential [J]'
                if data['1']['solution phases'][j]['phase model'] in ['SUBG', 'SUBQ']:
                    speciesLabel = 'quadruplets'
                else:
                    speciesLabel = 'species'
                for k in list(data['1']['solution phases'][j][speciesLabel].keys()):
                    ykey.append(['solution phases',j,speciesLabel,k,values['-yaxis-']])
                    yen.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey[yi][-2],key=str(yi))])
                    leg.append(j+': '+k)
                    yi = yi + 1
            phaseColumns.append([[sg.Text('Pure Condensed Phases')]])
            for j in pureCondensedPhases:
                ykey.append(['pure condensed phases',j,values['-yaxis-']])
                yen.append(False)
                phaseColumns[-1].append([sg.Checkbox(ykey[yi][-2],key=str(yi))])
                leg.append(j)
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
            ylab = 'Mole Fraction'
            solutionPhases = list(data['1']['solution phases'].keys())
            phaseColumns = []
            yi = 0
            for j in solutionPhases:
                phaseColumns.append([[sg.Text(j)]])
                if data['1']['solution phases'][j]['phase model'] in ['SUBG', 'SUBQ']:
                    speciesLabel = 'quadruplets'
                else:
                    speciesLabel = 'species'
                for k in list(data['1']['solution phases'][j][speciesLabel].keys()):
                    ykey.append(['solution phases',j,speciesLabel,k,values['-yaxis-']])
                    yen.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey[yi][-2],key=str(yi))])
                    leg.append(j+': '+k)
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
            ylab = 'Element Potential [J]'
            elements = list(data['1']['elements'].keys())
            for j in elements:
                ykey.append(['elements',j,values['-yaxis-']])
                yen.append(True)
                leg.append(j)
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
                xlab = 'Iteration'
            else:
                x.append(data[j][xkey])
                if xkey == 'temperature':
                    xlab = 'Temperature [K]'
                elif xkey == 'pressure':
                    xlab = 'Pressure [atm]'
            for yi in range(len(ykey)):
                if len(ykey[yi]) == 1:
                    y[yi].append(data[j][ykey[yi][0]])
                elif len(ykey[yi]) == 3:
                    y[yi].append(data[j][ykey[yi][0]][ykey[yi][1]][ykey[yi][2]])
                elif len(ykey[yi]) == 5:
                    y[yi].append(data[j][ykey[yi][0]][ykey[yi][1]][ykey[yi][2]][ykey[yi][3]][ykey[yi][4]])
        plotWindow.Element('Export Plot Script').Update(disabled = False)
        # Start figure
        fig = plt.figure()
        ax  = fig.add_axes([0.2, 0.1, 0.75, 0.85])
        for yi in range(len(yen)):
            if yen[yi]:
                ax.plot(x,y[yi],'.-',label = leg[yi])
        ax.legend()
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        plt.show()
    elif event == 'Export Plot Script':
        with open('python/generatedPlotScript.py', 'w') as f:
            f.write('# Thermochimica-generated plot script\n')
            f.write('import matplotlib.pyplot as plt\n')
            f.write('x = '+"{}\n".format(x))
            f.write('y = '+"{}\n".format(y))
            f.write('xlab = \''+xlab+'\'\n')
            f.write('ylab = \''+ylab+'\'\n')
            f.write('yen = '+"{}\n".format(yen))
            f.write('leg = '+"{}\n".format(leg))
            f.write('# Start figure\n')
            f.write('fig = plt.figure()\n')
            f.write('ax  = fig.add_axes([0.2, 0.1, 0.75, 0.85])\n')
            f.write('for yi in range(len(yen)):\n')
            f.write('    if yen[yi]:\n')
            f.write('        ax.plot(x,y[yi],\'.-\',label = leg[yi])\n')
            f.write('ax.legend()\n')
            f.write('ax.set_xlabel(xlab)\n')
            f.write('ax.set_ylabel(ylab)\n')
            f.write('plt.show()\n')
plotWindow.close()
