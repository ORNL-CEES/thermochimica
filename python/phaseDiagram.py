import json
import matplotlib.pyplot as plt
import numpy as np

def processPhaseDiagramData(fname,elx):
    f = open(fname,)
    data = json.load(f)
    f.close()
    if list(data.keys())[0] != '1':
        print('Output does not contain data series')
        exit()
    ts = []
    x1 = []
    x2 = []
    p1 = []
    p2 = []
    mint = 1e6
    maxt = 0
    for i in list(data.keys()):
        mint = min(mint,data[i]['temperature'])
        maxt = max(maxt,data[i]['temperature'])
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
    return ts, x1, x2, p1, p2, mint, maxt

fname = 'pd-ru-out.json'
elx = 'Ru'

ts = []
x1 = []
x2 = []
p1 = []
p2 = []
mint = 1e6
maxt = 0

tsr, x1r, x2r, p1r, p2r, mintr, maxtr = processPhaseDiagramData(fname,elx)
ts = [item for sublist in [ts,tsr] for item in sublist]
x1 = [item for sublist in [x1,x1r] for item in sublist]
x2 = [item for sublist in [x2,x2r] for item in sublist]
p1 = [item for sublist in [p1,p1r] for item in sublist]
p2 = [item for sublist in [p2,p2r] for item in sublist]
mint = min(mint,mintr)
maxt = max(maxt,maxtr)

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
    minj = np.argmin(np.array(ts)[inds])
    maxj = np.argmax(np.array(ts)[inds])
    if (np.array(ts)[inds][minj] > mint):
        ax.plot([np.array(x1)[inds][minj],np.array(x2)[inds][minj]],[np.array(ts)[inds][minj],np.array(ts)[inds][minj]],'k-')
    if (np.array(ts)[inds][maxj] < maxt):
        ax.plot([np.array(x1)[inds][maxj],np.array(x2)[inds][maxj]],[np.array(ts)[inds][maxj],np.array(ts)[inds][maxj]],'k-')
ax.set_xlim(0,1)
plt.show()
