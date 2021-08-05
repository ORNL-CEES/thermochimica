import PySimpleGUI as sg
import subprocess
import math
import os
import json
import matplotlib.pyplot as plt

f = open('thermoout.json',)

data = json.load(f)
f.close()

xkey = 'temperature'
ykey = [[]]
ykey[0].append('solution phases')
ykey[0].append('gas_ideal')
ykey[0].append('moles')
ykey.append([])
ykey[1].append('pure condensed phases')
ykey[1].append('C_Graphite(s)')
ykey[1].append('moles')
ykey.append([])
ykey[2].append('pressure')
ykey.append([])
ykey[3].append('solution phases')
ykey[3].append('gas_ideal')
ykey[3].append('species')
ykey[3].append('CO')
ykey[3].append('mole fraction')
x = []
y = []
for yi in range(len(ykey)):
    y.append([])
print(len(ykey))
for j in data.keys():
    x.append(data[j][xkey])
    for yi in range(len(ykey)):
        if len(ykey[yi]) == 1:
            y[yi].append(data[j][ykey[yi][0]])
        elif len(ykey[yi]) == 3:
            y[yi].append(data[j][ykey[yi][0]][ykey[yi][1]][ykey[yi][2]])
        elif len(ykey[yi]) == 5:
            y[yi].append(data[j][ykey[yi][0]][ykey[yi][1]][ykey[yi][2]][ykey[yi][3]][ykey[yi][4]])
print(x)
print(y)
# Start figure
fig = plt.figure()
ax  = fig.add_axes([0.12, 0.1, 0.85, 0.85])
for yi in range(len(ykey)):
    ax.plot(x,y[yi],'.-')
plt.show()
