import PySimpleGUI as sg
import matplotlib.pyplot as plt
import thermoTools
import thermoToolsGUI

class PlotWindow:
    def __init__(self,datafile):
        windowList.append(self)
        self.datafile = datafile
        self.ykey = [[]]
        self.yen = []
        self.leg = []
        self.showLeg = True
        self.ykey2 = [[]]
        self.yen2 = []
        self.leg2 = []
        self.showLeg2 = True
        self.xkey = []
        self.xlog = False
        self.ylog = False
        self.ylog2 = False
        self.readDatabase()
        optionsLayout = [
                          [sg.Text('x-axis')],[sg.Combo(['iteration', 'temperature', 'pressure'], default_value='iteration', key='-xaxis-')],[sg.Checkbox('Log scale',key='-xlog-')],
                          [sg.Text('y-axis')],[sg.Combo(['temperature', 'pressure', 'moles', 'mole fraction', 'chemical potential', 'driving force', 'vapor pressure',
                           'moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase','mole fraction of endmembers',
                           'moles of elements', 'element potential', 'integral Gibbs energy', 'functional norm', 'GEM iterations', '# phases', 'heat capacity','enthalpy','entropy'],
                            key='-yaxis-', enable_events=True)],[sg.Checkbox('Log scale',key='-ylog-')],
                          [sg.Text('y-axis 2')],[sg.Combo(['','temperature', 'pressure', 'moles', 'mole fraction', 'chemical potential', 'driving force', 'vapor pressure',
                           'moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase','mole fraction of endmembers',
                           'moles of elements', 'element potential', 'integral Gibbs energy', 'functional norm', 'GEM iterations', '# phases', 'heat capacity','enthalpy','entropy'],
                            key='-yaxis2-', enable_events=True)],[sg.Checkbox('Log scale',key='-ylog2-')]
                        ]
        plotLayout = [optionsLayout,
                      [sg.Column([[sg.Button('Plot', size = thermoToolsGUI.buttonSize)],
                                  [sg.Button('Plot Settings', size = thermoToolsGUI.buttonSize)],
                                  [sg.Button('Refresh Data', size = thermoToolsGUI.buttonSize)]
                                 ],vertical_alignment='t'),
                      sg.Column([[sg.Button('Export Plot Script', disabled = True, size = thermoToolsGUI.buttonSize)],
                                 [sg.Button('Export Plot', disabled = True, size = thermoToolsGUI.buttonSize)]                                  
                                 ],vertical_alignment='t')]]
        self.sgw = sg.Window('Thermochimica plot setup', plotLayout, location = [500,0], finalize=True)
        self.children = []
        self.yWindow = []
        self.yWindow2 = []
        self.currentPlot = []
        self.figureList = []
        self.exportFormat = 'png'
        self.exportFileName = 'plot'
        self.plotMarker = '.-'
        self.plotMarker2 = '*--'
        self.plotColor = 'colorful'
        self.plotColor2 = 'colorful'
        self.exportDPI = 300
        self.plotScriptFilename = 'python/generatedPlotScript.py'
    def close(self):
        for child in self.children:
            child.close()
        for fig in self.figureList:
            plt.close(fig=fig)
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event == '-yaxis-':
            self.ykey = []
            self.yen = []
            self.leg = []
            self.set_y_axis(values['-yaxis-'],self.ykey,self.yen,self.leg,self.yWindow)
        elif event == '-yaxis2-':
            self.ykey2 = []
            self.yen2 = []
            self.leg2 = []
            self.set_y_axis(values['-yaxis2-'],self.ykey2,self.yen2,self.leg2,self.yWindow2,offset=500)
        elif event == 'Plot':
            self.xkey  = values['-xaxis-']
            self.xlog  = values['-xlog-']
            self.ylog  = values['-ylog-']
            self.ylog2 = values['-ylog2-']
            self.makePlot()
        elif event == 'Export Plot Script':
            self.exportPlotScript()
        elif event == 'Export Plot':
            self.exportPlot()
        elif event == 'Plot Settings':
            settingsWindow = SettingsWindow(self)
            self.children.append(settingsWindow)
        elif event == 'Refresh Data':
            self.readDatabase()
    def set_y_axis(self,value,ykey,yen,leg,yWindow,offset=0):
        for window in yWindow:
            window.close()
        if value in ['temperature','pressure','integral Gibbs energy','functional norm','GEM iterations','heat capacity','enthalpy','entropy']:
            try:
                ykey.append([value])
                yen.append(True)
                leg.append(value)
            except:
                return
        elif value == '# phases':
            try:
                ykey.append(['# solution phases'])
                yen.append(True)
                leg.append('# of Solution Phases')
                ykey.append(['# pure condensed phases'])
                yen.append(True)
                leg.append('# of Pure Condensed Phases')
            except:
                return
        elif value in ['moles','chemical potential']:
            try:
                solutionPhases = list(self.data['1']['solution phases'].keys())
                pureCondensedPhases = list(self.data['1']['pure condensed phases'].keys())
            except:
                return
            phaseOptions = {}
            yi = 0
            for phase in solutionPhases:
                phaseOptions[phase] = []
                if value == 'moles':
                    # total moles of solution phase
                    try:
                        if self.data['1']['solution phases'][phase]['phase model'] in ['SUBG', 'SUBQ']:
                            ykey.append(['solution phases',phase,'moles of endmembers'])
                            yen.append(False)
                            phaseOptions[phase].append(['Moles of Endmembers',yi])
                            leg.append(phase)
                            yi = yi + 1
                        ykey.append(['solution phases',phase,value])
                        yen.append(False)
                        phaseOptions[phase].append(['Moles',yi])
                        leg.append(phase)
                        yi = yi + 1
                    except:
                        continue
                if self.data['1']['solution phases'][phase]['phase model'] in ['SUBG', 'SUBQ']:
                    speciesLabel = 'quadruplets'
                else:
                    speciesLabel = 'species'
                for k in list(self.data['1']['solution phases'][phase][speciesLabel].keys()):
                    try:
                        ykey.append(['solution phases',phase,speciesLabel,k,value])
                        yen.append(False)
                        phaseOptions[phase].append([ykey[yi][-2],yi])
                        leg.append(phase+': '+k)
                        yi = yi + 1
                    except:
                        continue
            for phase in pureCondensedPhases:
                phaseOptions[phase] = []
                try:
                    ykey.append(['pure condensed phases',phase,value])
                    yen.append(False)
                    phaseOptions[phase].append([ykey[yi][-2],yi])
                    leg.append(phase)
                    yi = yi + 1
                except:
                    continue
            selectWindow = SelectionWindow(phaseOptions,yen,offset=offset)
            self.children.append(selectWindow)
            yWindow.append(selectWindow)
        elif value in ['driving force']:
            try:
                solutionPhases = list(self.data['1']['solution phases'].keys())
                pureCondensedPhases = list(self.data['1']['pure condensed phases'].keys())
            except:
                return
            phaseOptions = {'Phases':[]}
            yi = 0
            for phase in solutionPhases:
                ykey.append(['solution phases',phase,value])
                yen.append(False)
                leg.append(phase)
                phaseOptions['Phases'].append([phase,yi])
                yi = yi + 1
            for phase in pureCondensedPhases:
                try:
                    ykey.append(['pure condensed phases',phase,value])
                    yen.append(False)
                    phaseOptions['Phases'].append([phase,yi])
                    leg.append(phase)
                    yi = yi + 1
                except:
                    continue
            selectWindow = SelectionWindow(phaseOptions,yen,offset=offset)
            self.children.append(selectWindow)
            yWindow.append(selectWindow)
        elif value in ['moles of element in phase', 'mole fraction of phase by element', 'mole fraction of element by phase']:
            solutionPhases = list(self.data['1']['solution phases'].keys())
            pureCondensedPhases = list(self.data['1']['pure condensed phases'].keys())
            phaseOptions = {}
            yi = 0
            for phase in solutionPhases:
                phaseOptions[phase] = []
                for element in list(self.data['1']['solution phases'][phase]['elements'].keys()):
                    try:
                        ykey.append(['solution phases',phase,'elements',element,value])
                        yen.append(False)
                        phaseOptions[phase].append([ykey[yi][-2],yi])
                        leg.append(phase+': '+element)
                        yi = yi + 1
                    except:
                        continue
            for phase in pureCondensedPhases:
                phaseOptions[phase] = []
                for element in list(self.data['1']['pure condensed phases'][phase]['elements'].keys()):
                    try:
                        ykey.append(['pure condensed phases',phase,'elements',element,value])
                        yen.append(False)
                        phaseOptions[phase].append([ykey[yi][-2],yi])
                        leg.append(phase+': '+element)
                        yi = yi + 1
                    except:
                        continue
            selectWindow = SelectionWindow(phaseOptions,yen,offset=offset)
            self.children.append(selectWindow)
            yWindow.append(selectWindow)
        elif value == 'mole fraction':
            solutionPhases = list(self.data['1']['solution phases'].keys())
            phaseOptions = {}
            yi = 0
            for phase in solutionPhases:
                phaseOptions[phase] = []
                if self.data['1']['solution phases'][phase]['phase model'] in ['SUBG', 'SUBQ']:
                    speciesLabel = 'quadruplets'
                else:
                    speciesLabel = 'species'
                for species in list(self.data['1']['solution phases'][phase][speciesLabel].keys()):
                    try:
                        ykey.append(['solution phases',phase,speciesLabel,species,value])
                        yen.append(False)
                        phaseOptions[phase].append([ykey[yi][-2],yi])
                        leg.append(phase+': '+species)
                        yi = yi + 1
                    except:
                        continue
            selectWindow = SelectionWindow(phaseOptions,yen,offset=offset)
            self.children.append(selectWindow)
            yWindow.append(selectWindow)
        elif value == 'mole fraction of endmembers':
            solutionPhases = list(self.data['1']['solution phases'].keys())
            phaseOptions = {}
            yi = 0
            for phase in solutionPhases:
                if self.data['1']['solution phases'][phase]['phase model'] in ['SUBG', 'SUBQ']:
                    phaseOptions[phase] = []
                    for endmember in list(self.data['1']['solution phases'][phase]['endmembers'].keys()):
                        try:
                            ykey.append(['solution phases',phase,'endmembers',endmember,'mole fraction'])
                            yen.append(False)
                            phaseOptions[phase].append([ykey[yi][-2],yi])
                            leg.append(phase+': '+endmember)
                            yi = yi + 1
                        except:
                            continue
                else:
                    continue
            selectWindow = SelectionWindow(phaseOptions,yen,offset=offset)
            self.children.append(selectWindow)
            yWindow.append(selectWindow)
        elif value == 'vapor pressure':
            solutionPhases = list(self.data['1']['solution phases'].keys())
            phaseOptions = {'Vapor Pressures':[]}
            yi = 0
            for phase in solutionPhases:
                if self.data['1']['solution phases'][phase]['phase model'] != 'IDMX':
                    break
                for species in list(self.data['1']['solution phases'][phase]['species'].keys()):
                    try:
                        ykey.append(['solution phases',phase,'species',species,value])
                        yen.append(False)
                        phaseOptions['Vapor Pressures'].append([ykey[yi][-2],yi])
                        leg.append(phase+': '+species)
                        yi = yi + 1
                    except:
                        continue
                break
            selectWindow = SelectionWindow(phaseOptions,yen,offset=offset)
            self.children.append(selectWindow)
            yWindow.append(selectWindow)
        elif value == 'moles of elements':
            elements = list(self.data['1']['elements'].keys())
            elementOptions = {'Elements':[]}
            yi = 0
            for element in elements:
                try:
                    ykey.append(['elements',element,'moles'])
                    yen.append(False)
                    elementOptions['Elements'].append([element,yi])
                    leg.append(element)
                    yi = yi + 1
                except:
                    continue
            selectWindow = SelectionWindow(elementOptions,yen,offset=offset)
            self.children.append(selectWindow)
            yWindow.append(selectWindow)
        elif value == 'element potential':
            elements = list(self.data['1']['elements'].keys())
            elementOptions = {'Elements':[]}
            yi = 0
            for element in elements:
                try:
                    ykey.append(['elements',element,value])
                    yen.append(False)
                    elementOptions['Elements'].append([element,yi])
                    leg.append(element)
                    yi = yi + 1
                except:
                    continue
            selectWindow = SelectionWindow(elementOptions,yen,offset=offset)
            self.children.append(selectWindow)
            yWindow.append(selectWindow)
    def makePlot(self):
        # Select data
        yused, legend, yused2, legend2 = thermoTools.selectData(self.yen,self.ykey,self.leg,yen2=self.yen2,ykey2=self.ykey2,leg2=self.leg2)
        # Check if legends should be displayed
        if not self.showLeg:
            legend = None
        if not self.showLeg2:
            legend2 = None
        # Call plotter
        self.currentPlot = thermoTools.makePlot(self.datafile,self.xkey,yused,legend=legend,yused2=yused2,legend2=legend2,plotColor=self.plotColor,plotColor2=self.plotColor2,plotMarker=self.plotMarker,plotMarker2=self.plotMarker2,xlog=self.xlog,ylog=self.ylog,ylog2=self.ylog2,interactive=True)
        
        self.figureList.append(self.currentPlot)

        # Update buttons
        self.sgw.Element('Export Plot').Update(disabled = False)
        self.sgw.Element('Export Plot Script').Update(disabled = False)
    def exportPlotScript(self):
        # Select data
        yused, legend, yused2, legend2 = thermoTools.selectData(self.yen,self.ykey,self.leg,yen2=self.yen2,ykey2=self.ykey2,leg2=self.leg2)
        # Check if legends should be displayed
        if not self.showLeg:
            legend = None
        if not self.showLeg2:
            legend2 = None
        # Call plot exporter
        thermoTools.exportPlotScript(self.plotScriptFilename,self.datafile,self.xkey,yused,legend=legend,yused2=yused2,legend2=legend2,plotColor=self.plotColor,plotColor2=self.plotColor2,plotMarker=self.plotMarker,plotMarker2=self.plotMarker2,xlog=self.xlog,ylog=self.ylog,ylog2=self.ylog2)
    def exportPlot(self):
        try:
            self.currentPlot.savefig(f'outputs/{self.exportFileName}.{self.exportFormat}', format=self.exportFormat, dpi=self.exportDPI)
        except:
            errorLayout = [[sg.Text('The export failed, try changing plot settings.')],[sg.Button('Continue'), sg.Button('Cancel')]]
            errorWindow = sg.Window('Plot export failed', errorLayout, location = [400,0], finalize=True, keep_on_top = True)
            while True:
                event, values = errorWindow.read(timeout=thermoToolsGUI.timeout)
                if event == sg.WIN_CLOSED or event == 'Continue':
                    break
            errorWindow.close()
    def readDatabase(self):
        self.data = thermoTools.readDatabase(self.datafile)
        if list(self.data.keys())[0] != '1':
            print('Output does not contain data series')
            exit()

class SelectionWindow:
    def __init__(self, options, yen, offset = 0):
        self.options = options
        self.yen = yen
        self.selectables = []
        self.selected = []
        self.drop_selection = ''
        self.makeLayout(offset)
        windowList.append(self)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def makeLayout(self,offset):
        dropdown_options = list(self.options.keys())
        try:
            self.drop_selection = dropdown_options[0]
        except IndexError:
            self.drop_selection = None
        selectables_column = [
            [
                sg.Text('Unselected Options'),
                sg.Combo(dropdown_options, default_value=self.drop_selection, key='-drop-', enable_events=True)
            ],
            [
                sg.Listbox(
                    values=[], enable_events=False, size=(40, 20), key="-selectables-", select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED
                )
            ],
            [sg.Button('Add Selected Options'),sg.Button('Add All')]
        ]
        selected_column = [
            [
                sg.Text('Selected Options')
            ],
            [
                sg.Listbox(
                    values=[], enable_events=False, size=(40, 20), key="-selected-", select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED
                )
            ],
            [sg.Button('Remove Selected Options'),sg.Button('Remove All')]
        ]
        selectionLayout = [[sg.Column(selectables_column,vertical_alignment='t'),sg.Column(selected_column,vertical_alignment='t')],[sg.Button('Exit')]]
        self.sgw = sg.Window('Plot Selection', selectionLayout, location = [800,offset], finalize=True)
        self.updateSelectables()
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event == '-drop-':
            self.drop_selection = values['-drop-']
            self.updateSelectables()
        elif event == 'Add Selected Options':
            for option in values['-selectables-']:
                o = option.copy()
                o.insert(0,self.drop_selection)
                self.selected.append(o)
                self.options[self.drop_selection].remove(option)
            self.updateSelected()
            self.updateSelectables()
        elif event == 'Remove Selected Options':
            for option in values['-selected-']:
                self.selected.remove(option)
                drop = option[0]
                self.options[drop].append([option[1],option[2]])
            self.updateSelected()
            self.updateSelectables()
        elif event == 'Add All':
            drop_options = list(self.options.keys())
            for drop_option in drop_options:
                options = list(self.options[drop_option])
                for option in options:
                    o = option.copy()
                    o.insert(0,drop_option)
                    self.selected.append(o)
                    self.options[drop_option].remove(option)
            self.updateSelected()
            self.updateSelectables()
        elif event == 'Remove All':
            options = list(self.selected)
            for option in options:
                self.selected.remove(option)
                drop = option[0]
                self.options[drop].append([option[1],option[2]])
            self.updateSelected()
            self.updateSelectables()
    def updateSelectables(self):
        if self.drop_selection:
            self.selectables = self.options[self.drop_selection]
            self.selectables.sort(key=lambda x: x[1])
            self.sgw['-selectables-'].update(self.selectables)
    def updateSelected(self):
        self.selected.sort(key=lambda x: x[2])
        self.sgw['-selected-'].update(self.selected)
        # Set enabled array status in parent
        for yi in range(len(self.yen)):
            self.yen[yi] = False
        for selected in self.selected:
            self.yen[selected[2]] = True

class SettingsWindow:
    def __init__(self, parent):
        self.parent = parent
        self.makeLayout()
        windowList.append(self)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def makeLayout(self):
        if self.parent.plotMarker == '-':
            line  = True
            point = False
            both  = False
        elif self.parent.plotMarker == '.':
            line  = False
            point = True
            both  = False
        else:
            line  = False
            point = False
            both  = True
        if self.parent.plotMarker2 == '--':
            line2  = True
            point2 = False
            both2  = False
        elif self.parent.plotMarker2 == '*':
            line2  = False
            point2 = True
            both2  = False
        else:
            line2  = False
            point2 = False
            both2  = True
        if self.parent.plotColor == 'colorful':
            colorful = True
            bland    = False
        else:
            colorful = False
            bland    = True
        if self.parent.plotColor2 == 'colorful':
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
                            [sg.Column([
                                        [sg.Checkbox('Show Legend', default=self.parent.showLeg, key='-showLegend-')]
                                        ],vertical_alignment='t'),
                            sg.Column([
                                        [sg.Checkbox('Show Legend 2', default=self.parent.showLeg2, key='-showLegend2-')]
                                        ],vertical_alignment='t')],
                            [sg.Text('Export Filename'),sg.Input(key='-filename-',size=(thermoToolsGUI.inputSize,1))],
                            [sg.Text('Export Format'),sg.Combo(['png', 'pdf', 'ps', 'eps', 'svg'],default_value='png',key='-format-')],
                            [sg.Text('Export DPI'),sg.Input(key='-dpi-',size=(thermoToolsGUI.inputSize,1))],
                            [sg.Button('Accept')]]
        self.sgw = sg.Window('Plot Settings', settingsLayout, location = [400,0], finalize=True)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
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
            self.parent.showLeg  = values['-showLegend-']
            self.parent.showLeg2 = values['-showLegend2-']
            try:
                if str(values['-filename-']) != '':
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

windowList = []
dataWindow = thermoToolsGUI.DataWindow(windowList,PlotWindow,thermoToolsGUI.JSONParse,ext='.json',rootDir='outputs')
while len(windowList) > 0:
    for window in windowList:
        window.read()
