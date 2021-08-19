import json
import matplotlib.pyplot as plt
import numpy as np

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

boundaries = []
b = []
for i in range(len(p1)):
    # If a miscibility gap label has been used unnecessarily, remove it
    if p1[i].find('#2') > 0:
        if not(p1[i][0:p1[i].find('#2')] == p2[i]):
            p1[i] = p1[i][0:p1[i].find('#2')]
    if p2[i].find('#2') > 0:
        if not(p2[i][0:p2[i].find('#2')] == p1[i]):
            p2[i] = p2[i][0:p2[i].find('#2')]

    repeat = False
    for j in range(len(boundaries)):
        if (boundaries[j][0] == p1[i]) and (boundaries[j][1] == p2[i]):
            b.append(j)
            repeat = True
    if not(repeat):
        boundaries.append([p1[i],p2[i]])
        b.append(len(boundaries)-1)
print(boundaries)
# print(b)
# Start figure
fig = plt.figure()
ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])
for j in range(len(boundaries)):
    inds = [i for i, k in enumerate(b) if k == j]
    ax.plot(np.array(x1)[inds],np.array(ts)[inds],'.')
    ax.plot(np.array(x2)[inds],np.array(ts)[inds],'.')
ax.set_xlim(0,1)
plt.show()
