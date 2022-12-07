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