import PySimpleGUI as sg
import numpy as np
import thermoToolsGUI
import sys

class pdPoint:
    def __init__(self,elements,temperature,concentration,phases,phaseConcentrations,energy,iterations):
        self.t = temperature
        self.runConcentration = concentration
        self.phaseConcentrations = phaseConcentrations
        self.phases = phases
        self.details = f'Temperature = {temperature:6.2f}\n'
        for elem,conc in zip(elements,self.runConcentration):
            self.details = self.details + f'Moles of {elem} = {conc:9.8f}\n'
        for phase,conc in zip(self.phases,self.phaseConcentrations):
            self.details = self.details + f'{phase} at {conc:5.4f}\n'
        self.details = self.details + f'Integral Gibbs Energy = {energy:.2f}\n'
        self.details = self.details + f'Number of GEM iterations = {iterations}'
        self.suppressed = False

class InspectWindow:
    def __init__(self,parent,endMember2,phases,windowList):
        self.parent = parent
        windowList.append(self)
        self.windowList = windowList
        dataColumn = [
            [sg.Text('Data Points')],
            [sg.Listbox(values=[], enable_events=True, size=(30, 50), key='-dataList-')]
        ]
        outputColumn = [
            [sg.Text('Calculation Details')],
            [sg.Multiline(key='-details-', size=(50,10), no_scrollbar=True)],
            [sg.Text(key = '-status-')],
            [sg.Button('Toggle Active/Suppressed Status', disabled = True)],
            [sg.Text('Filter points', font='underline')],
            [sg.Text('Temperature Range:')],
            [sg.Input(key='-tfilterlow-',size=(thermoToolsGUI.inputSize,1)),sg.Input(key='-tfilterhi-',size=(thermoToolsGUI.inputSize,1))],
            [sg.Text(f'{endMember2} Concentration Range:')],
            [sg.Input(key='-xfilterlow-',size=(thermoToolsGUI.inputSize,1)),sg.Input(key='-xfilterhi-',size=(thermoToolsGUI.inputSize,1))],
            [sg.Text('Contains Phases:')],
            [sg.Combo(['']+phases, key = '-pfilter1-'),sg.Combo(['']+phases, key = '-pfilter2-')],
            [sg.Text('Active/Suppressed Status:')],
            [sg.Combo(['','Active','Suppressed'], key = '-activefilter-')],
            [sg.Button('Apply Filter')]
        ]
        self.data = [[i, f'{point.t:6.2f} K {point.phaseConcentrations[0]:4.3f} {point.phaseConcentrations[1]:4.3f}'] for i,point in enumerate(self.parent.calculation.pdPoints)]
        self.sgw = sg.Window('Data inspection',
            [[sg.Pane([
                sg.Column(dataColumn, element_justification='l', expand_x=True, expand_y=True),
                sg.Column(outputColumn, element_justification='c', expand_x=True, expand_y=True)
            ], orientation='h', k='-PANE-')]],
            location = [0,0], finalize=True)
        self.sgw['-dataList-'].update(self.data)
        self.children = []
        self.index = -1
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in self.windowList:
            self.windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=thermoToolsGUI.timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event == '-dataList-':
            self.index = values['-dataList-'][0][0]
            self.point = self.parent.calculation.pdPoints[self.index]
            self.sgw['-details-'].update(self.point.details)
            self.sgw['Toggle Active/Suppressed Status'].update(disabled = False)
            self.sgw['-status-'].update(f'{"Suppressed" if self.point.suppressed else "Active"}')
        elif event == 'Toggle Active/Suppressed Status':
            if self.index >= 0:
                self.point.suppressed = not(self.point.suppressed)
                self.parent.macro.append(f'macroPD.suppressed[{self.index}] = not(macroPD.suppressed[{self.index}])')
                self.sgw['-status-'].update(f'{"Suppressed" if self.point.suppressed else "Active"}')
        elif event == 'Apply Filter':
            tlo = -np.Inf
            thi  = np.Inf
            xlo = -np.Inf
            xhi  = np.Inf
            try:
                tlo = float(values['-tfilterlow-'])
            except:
                pass
            try:
                thi = float(values['-tfilterhi-'])
            except:
                pass
            try:
                xlo = float(values['-xfilterlow-'])
            except:
                pass
            try:
                xhi = float(values['-xfilterhi-'])
            except:
                pass
            self.data = []
            for i,p in enumerate(self.parent.calculation.pdPoints):
                # Check temperature
                tfilt = tlo <= p.t and thi >= p.t
                # Check concentration
                xfilt = (xlo <= p.phaseConcentrations[0] and xhi >= p.phaseConcentrations[0]) or (xlo <= p.phaseConcentrations[1] and xhi >= p.phaseConcentrations[1])
                # Check phases present
                pfilt = (values['-pfilter1-'] == '' or values['-pfilter1-'] == p.phases[0] or values['-pfilter1-'] == p.phases[1]) and (values['-pfilter2-'] == '' or values['-pfilter2-'] == p.phases[0] or values['-pfilter2-'] == p.phases[1])
                # Check active/suppressed status
                afilt = (values['-activefilter-'] == '') or ((values['-activefilter-'] == 'Suppressed') == p.suppressed)
                # If all filters pass, add to display list
                if tfilt and xfilt and pfilt and afilt:
                    self.data.append([i, f'{p.t:6.2f} K {p.phaseConcentrations[0]:4.3f} {p.phaseConcentrations[1]:4.3f}'])
            self.sgw['-dataList-'].update(self.data)

def runMacro(calc):
    if 'macroPhaseDiagram' in sys.modules:
        del sys.modules['macroPhaseDiagram']
    import macroPhaseDiagram
    calc.calculation = macroPhaseDiagram.macroPD
    calc.calculation.active = True
    calc.calculation.interactivePlot = True
    enableButtons(calc)
    
def enableButtons(calc):
    calc.sgw.Element('Refine').Update(disabled = False)
    calc.sgw.Element('Auto Refine').Update(disabled = False)
    calc.sgw.Element('Auto Smoothen').Update(disabled = False)
    calc.sgw.Element('Add Label').Update(disabled = False)
    calc.sgw.Element('Auto Label').Update(disabled = False)
    calc.sgw.Element('Plot').Update(disabled = False)
    calc.sgw.Element('Undo').Update(disabled = False)
    calc.sgw.Element('Inspect').Update(disabled = False)
    calc.sgw.Element('Export Diagram Data').Update(disabled = False)
    calc.sgw.Element('Export Plot').Update(disabled = False)
