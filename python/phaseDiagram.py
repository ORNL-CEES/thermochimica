import math
import json
import matplotlib.pyplot as plt

f = open('pd-ru-out.json',)
data = json.load(f)
f.close()
if list(data.keys())[0] != '1':
    print('Output does not contain data series')
    exit()

phases = list(data['1']['solution phases'].keys())
for phaseName in list(data['1']['pure condensed phases'].keys()):
    phases.append(phaseName)
print(phases)
elements = list(data['1']['elements'].keys())
print(elements)
elx = 'Ru'

ts = []
x1 = []
x2 = []
p1 = []
p2 = []
for i in list(data.keys()):
    # print(data[i]['temperature'])
    if (data[i]['# solution phases'] + data[i]['# pure condensed phases']) == 2:
        ts.append(data[i]['temperature'])
        boundPhases = []
        boundComps = []
        for phaseName in list(data[i]['solution phases'].keys()):
            if (data[i]['solution phases'][phaseName]['moles'] > 0):
                boundPhases.append(phaseName)
                boundComps.append(data[i]['solution phases'][phaseName]['elements'][elx]['mole fraction of phase by element'])
        for phaseName in list(data[i]['pure condensed phases'].keys()):
            if (data[i]['pure condensed phases'][phaseName]['moles'] > 0):
                boundPhases.append(phaseName)
                boundComps.append(data[i]['pure condensed phases'][phaseName]['elements'][elx]['mole fraction of phase by element'])
        x1.append(boundComps[0])
        x2.append(boundComps[1])
        p1.append(boundPhases[0])
        p2.append(boundPhases[1])

# Start figure
fig = plt.figure()
ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])
ax.plot(x1,ts,'.')
ax.plot(x2,ts,'.')
plt.show()
