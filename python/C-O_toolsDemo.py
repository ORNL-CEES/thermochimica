import thermoTools

scriptName = 'inputs/co-script-demo.ti'
datafile   = 'data/C-O.dat'
jsonName   = 'c-o.json'
elements   = ['C','O']
tstart     = 300
tend       = 1500
ntstep     = 12
pstart     = 1
pend       = 1
npstep     = 1
masses     = [1,1.5]

thermoTools.WriteInputScript(scriptName,datafile,elements,tstart,tend,ntstep,pstart,pend,npstep,masses,heatCapacity=True)
thermoTools.RunInputScript(scriptName,jsonName=jsonName)

xkey = 'temperature'
yused = [['solution phases', 'gas_ideal','moles'],['pure condensed phases', 'C_Graphite(s)','moles']]
y2used = [['heat capacity']]
leg = ['Ideal Gas','Graphite']

thermoTools.makePlot(jsonName,xkey,yused,legend=leg,yused2=y2used,plotColor2='bland')