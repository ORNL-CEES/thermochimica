import pd

testcalc = pd.diagram([],'/media/max/data/thermochimicastuff/thermochimica/data/Kaye_NobleMetals.dat', False)
testcalc.run(10,10,1,'K','atm',0,1,300,3000,'Pd','Ru','mole fraction')

functions = [testcalc.refinery,testcalc.autoSmooth,testcalc.makePlot]

[f() for f in functions]
