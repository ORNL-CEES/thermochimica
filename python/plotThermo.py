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
                  [sg.Text('y-axis')],[sg.Combo(['temperature', 'pressure', 'moles', 'mole fraction', 'chemical potential',
                   'moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase',
                   'element potential', 'integral Gibbs energy', 'functional norm', '# phases'],
                    key='-yaxis-', enable_events=True)],
                  [sg.Text('y-axis')],[sg.Combo(['temperature', 'pressure', 'moles', 'mole fraction', 'chemical potential',
                   'moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase',
                   'element potential', 'integral Gibbs energy', 'functional norm', '# phases'],
                    key='-yaxis2-', enable_events=True, disabled=True)]
                ]
plotLayout = [optionsLayout,
              [sg.Button('Plot', disabled = True), sg.Button('Export Plot Script', disabled = True)]]
plotWindow = sg.Window('Thermochimica plot setup', plotLayout, finalize=True)

ykey = [[]]
yen = []
leg = []
ykey2 = [[]]
yen2 = []
leg2 = []

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
            plotWindow.Element('-yaxis2-').Update(disabled = False)
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
            plotWindow.Element('-yaxis2-').Update(disabled = False)
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
                    plotWindow.Element('-yaxis2-').Update(disabled = False)
                    break
            selectWindow.close()
        elif values['-yaxis-'] in ['moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase']:
            ykey = []
            solutionPhases = list(data['1']['solution phases'].keys())
            pureCondensedPhases = list(data['1']['pure condensed phases'].keys())
            phaseColumns = []
            yi = 0
            for j in solutionPhases:
                phaseColumns.append([[sg.Text(j)]])
                if values['-yaxis-'] == 'moles of element in phase':
                    ylab = 'Moles of Element in Phase'
                if values['-yaxis-'] == 'mole fraction of phase by element':
                    ylab = 'Mole Fraction of Phase by Element'
                else:
                    ylab = 'Mole Fraction of Element by Phase'
                for k in list(data['1']['solution phases'][j]['elements'].keys()):
                    ykey.append(['solution phases',j,'elements',k,values['-yaxis-']])
                    yen.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey[yi][-2],key=str(yi))])
                    leg.append(j+': '+k)
                    yi = yi + 1
            for j in pureCondensedPhases:
                phaseColumns.append([[sg.Text(j)]])
                for k in list(data['1']['pure condensed phases'][j]['elements'].keys()):
                    ykey.append(['pure condensed phases',j,'elements',k,values['-yaxis-']])
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
                    plotWindow.Element('-yaxis2-').Update(disabled = False)
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
                    plotWindow.Element('-yaxis2-').Update(disabled = False)
                    break
            selectWindow.close()
        elif values['-yaxis-'] == 'element potential':
            ykey = []
            ylab = 'Element Potential [J]'
            elements = list(data['1']['elements'].keys())
            for j in elements:
                ykey.append(['elements',j,values['-yaxis-']])
                yen.append(True)
                leg.append(j)
            plotWindow.Element('Plot').Update(disabled = False)
            plotWindow.Element('-yaxis2-').Update(disabled = False)
    elif event == '-yaxis2-':
        ykey2 = [[]]
        yen2 = []
        leg2 = []
        if values['-yaxis2-'] in ['temperature','pressure','integral Gibbs energy','functional norm']:
            ykey2[0].append(values['-yaxis2-'])
            yen2.append(True)
            leg2.append(values['-yaxis2-'])
            if values['-yaxis2-'] == 'temperature':
                ylab2 = 'Temperature [K]'
            elif values['-yaxis2-'] == 'pressure':
                ylab2 = 'Pressure [atm]'
            elif values['-yaxis2-'] == 'integral Gibbs energy':
                ylab2 = 'Integral Gibbs Energy [J]'
            elif values['-yaxis2-'] == 'functional norm':
                ylab2 = 'Functional Norm'
            plotWindow.Element('Plot').Update(disabled = False)
        elif values['-yaxis2-'] == '# phases':
            ykey2[0].append('# solution phases')
            yen2.append(True)
            leg2.append('# of Solution Phases')
            ykey2.append([])
            ykey2[1].append('# pure condensed phases')
            yen2.append(True)
            leg2.append('# of Pure Condensed Phases')
            ylab2 = 'Number of Stable Phases'
            plotWindow.Element('Plot').Update(disabled = False)
        elif values['-yaxis2-'] in ['moles','chemical potential']:
            ykey2 = []
            solutionPhases = list(data['1']['solution phases'].keys())
            pureCondensedPhases = list(data['1']['pure condensed phases'].keys())
            phaseColumns = []
            yi = 0
            for j in solutionPhases:
                phaseColumns.append([[sg.Text(j)]])
                if values['-yaxis2-'] == 'moles':
                    # total moles of solution phase
                    ykey2.append(['solution phases',j,values['-yaxis2-']])
                    yen2.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey2[yi][-2],key=str(yi))])
                    leg2.append(j)
                    ylab2 = 'Moles'
                    yi = yi + 1
                else:
                    ylab2 = 'Chemical Potential [J]'
                if data['1']['solution phases'][j]['phase model'] in ['SUBG', 'SUBQ']:
                    speciesLabel = 'quadruplets'
                else:
                    speciesLabel = 'species'
                for k in list(data['1']['solution phases'][j][speciesLabel].keys()):
                    ykey2.append(['solution phases',j,speciesLabel,k,values['-yaxis2-']])
                    yen2.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey2[yi][-2],key=str(yi))])
                    leg2.append(j+': '+k)
                    yi = yi + 1
            phaseColumns.append([[sg.Text('Pure Condensed Phases')]])
            for j in pureCondensedPhases:
                ykey2.append(['pure condensed phases',j,values['-yaxis2-']])
                yen2.append(False)
                phaseColumns[-1].append([sg.Checkbox(ykey2[yi][-2],key=str(yi))])
                leg2.append(j)
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
                    for yi in range(len(ykey2)):
                        yen2[yi] = values[str(yi)]
                    plotWindow.Element('Plot').Update(disabled = False)
                    break
            selectWindow.close()
        elif values['-yaxis2-'] in ['moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase']:
            ykey2 = []
            solutionPhases = list(data['1']['solution phases'].keys())
            pureCondensedPhases = list(data['1']['pure condensed phases'].keys())
            phaseColumns = []
            yi = 0
            for j in solutionPhases:
                phaseColumns.append([[sg.Text(j)]])
                if values['-yaxis2-'] == 'moles of element in phase':
                    ylab2 = 'Moles of Element in Phase'
                if values['-yaxis2-'] == 'mole fraction of phase by element':
                    ylab2 = 'Mole Fraction of Phase by Element'
                else:
                    ylab2 = 'Mole Fraction of Element by Phase'
                for k in list(data['1']['solution phases'][j]['elements'].keys()):
                    ykey2.append(['solution phases',j,'elements',k,values['-yaxis2-']])
                    yen2.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey2[yi][-2],key=str(yi))])
                    leg2.append(j+': '+k)
                    yi = yi + 1
            for j in pureCondensedPhases:
                phaseColumns.append([[sg.Text(j)]])
                for k in list(data['1']['pure condensed phases'][j]['elements'].keys()):
                    ykey2.append(['pure condensed phases',j,'elements',k,values['-yaxis2-']])
                    yen2.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey2[yi][-2],key=str(yi))])
                    leg2.append(j+': '+k)
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
                    for yi in range(len(ykey2)):
                        yen2[yi] = values[str(yi)]
                    plotWindow.Element('Plot').Update(disabled = False)
                    break
            selectWindow.close()
        elif values['-yaxis2-'] == 'mole fraction':
            ykey2 = []
            ylab2 = 'Mole Fraction'
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
                    ykey2.append(['solution phases',j,speciesLabel,k,values['-yaxis2-']])
                    yen2.append(False)
                    phaseColumns[-1].append([sg.Checkbox(ykey2[yi][-2],key=str(yi))])
                    leg2.append(j+': '+k)
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
                    for yi in range(len(ykey2)):
                        yen2[yi] = values[str(yi)]
                    plotWindow.Element('Plot').Update(disabled = False)
                    break
            selectWindow.close()
            plotWindow.Element('Plot').Update(disabled = False)
        elif values['-yaxis2-'] == 'element potential':
            ykey2 = []
            ylab2 = 'Element Potential [J]'
            elements = list(data['1']['elements'].keys())
            for j in elements:
                ykey2.append(['elements',j,values['-yaxis2-']])
                yen2.append(True)
                leg2.append(j)
            plotWindow.Element('Plot').Update(disabled = False)
    elif event == 'Plot':
        x = []
        y = []
        y2 = []
        xkey = values['-xaxis-']
        for yi in range(len(ykey)):
            y.append([])
        for yi in range(len(ykey2)):
            y2.append([])
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
            for yi in range(len(ykey2)):
                if len(ykey2[yi]) == 1:
                    y2[yi].append(data[j][ykey2[yi][0]])
                elif len(ykey2[yi]) == 3:
                    y2[yi].append(data[j][ykey2[yi][0]][ykey2[yi][1]][ykey2[yi][2]])
                elif len(ykey2[yi]) == 5:
                    y2[yi].append(data[j][ykey2[yi][0]][ykey2[yi][1]][ykey2[yi][2]][ykey2[yi][3]][ykey2[yi][4]])
        plotWindow.Element('Export Plot Script').Update(disabled = False)
        # Start figure
        fig = plt.figure()
        lns=[]
        if True in yen2:
            ax = fig.add_axes([0.2, 0.1, 0.65, 0.85])
        else:
            ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])
        for yi in range(len(yen)):
            if yen[yi]:
                lns = lns + ax.plot(x,y[yi],'.-',label = leg[yi])
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        if True in yen2:
            ax2 = ax.twinx()
            for yi in range(len(yen2)):
                if yen2[yi]:
                    lns = lns + ax2.plot(x,y2[yi],'^--',label = leg2[yi])
            ax2.set_ylabel(ylab2)
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, loc=0)
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
            f.write('lns=[]\n')
            f.write('# Start figure\n')
            f.write('fig = plt.figure()\n')
            if True in yen2:
                f.write('ax  = fig.add_axes([0.2, 0.1, 0.65, 0.85])\n')
            else:
                f.write('ax  = fig.add_axes([0.2, 0.1, 0.75, 0.85])\n')
            f.write('for yi in range(len(yen)):\n')
            f.write('    if yen[yi]:\n')
            f.write('        lns = lns + ax.plot(x,y[yi],\'.-\',label = leg[yi])\n')
            if True in yen2:
                f.write('y2 = '+"{}\n".format(y2))
                f.write('ylab2 = \''+ylab2+'\'\n')
                f.write('yen2 = '+"{}\n".format(yen2))
                f.write('leg2 = '+"{}\n".format(leg2))
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
plotWindow.close()
