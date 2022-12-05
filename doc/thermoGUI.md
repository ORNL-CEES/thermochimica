# `thermoGUI`
`thermoGUI` is used to run equilibrium calculations using Thermochimica.  It can be launched by using the following command:

```bash
./scripts/thermoGUI.sh
```

## Simple calculations
The default `thermoGUI` calculation window opened when the database `Kaye_NobleMetals.dat` is selected is shown below.

![Default `thermoGUI` calculation window](/doc/images/thermoGUI-default.png)

In the default configuration, a single Thermochimica calculation can be run. The system conditions should be specified in the appropriate input boxes, and units can be selected from the drop-down menus. Values can be entered in decimal (`1234.5`) or scientific (`1.2345e3`) notation. All boxes have a default value if none is entered (or an invalid entry is supplied). For temperature, this is 300, for pressure, 1, and for each composition, 0. Therefore at a minimum, one element composition must be entered manually.

When a calculation is run (by pressing `Run`), the results are displayed in an output window as shown below. This output corresponds to the terminal output from Thermochimica calculations with `iPrintResultsMode = 2`.

![Default `thermoGUI` output window](/doc/images/thermoGUI-output-default.png)

There are three checkboxes towards the bottom of the `thermoGUI` calculation window, `Save JSON`, `Calculate heat capacity, entropy, and enthalpy`, and `Use fuzzy stoichiometry`. When `Save JSON` is selected, a JSON database containing all calculation outputs will be generated and placed in the `outputs` directory. The default name for this file is `thermoout.json`, but a different name can be entered by using the `Set name` button. Note the extension `.json` will be appended automatically to the entered name. Selecting `Calculate heat capacity, entropy, and enthalpy` requests Thermochimica to compute these values, which will be both appended to the text output and included in the exported JSON database, if selected. When `Use fuzzy stoichiometry` is selected, noise will be applied to the species stoichiometries at input. This is useful for calculations on under-determined systems (including some molten salt calculations). Generally this option should not be used unless the calculation fails to converge to a reasonable assemblage consistently without it.

## Loops over temperature and pressure
If the radio buttons under `Temperature range:` or `Pressure range:` are set to `Enabled`, additional input boxes are revealed, as shown below.

![`thermoGUI` window with temperature loop input.](/doc/images/thermoGUI-tloop.png)

Now, the original `Temperature` or `Pressure` input is taken as a start point, an endpoint is supplied in the new `End Temperature` or `End Pressure` input box, and a number of steps is supplied in the `# of steps` box. Note that the start point is not counted as a step, so the step size is `(end - start)/(# of steps)`. For example, setting `Temperature = 400`, `End Temperature = 600`, and `# of steps = 2` results in calculations at 400 K, 500 K, and 600 K. The default value for `# of steps` is 10.

Loops over temperature and pressure may be used simultaneously. If both are set to `Enabled`, then the loops will be nested. Alternatively, under `Pressure range:` there is an option `Enabled, step with temperature`, which constructs a single loop in which both temperature and pressure are incremented from one endpoint to the other.

## Loops over composition
Similarly, a loop from one composition endpoint to another can be run by selecting `Enabled` under `Composition range:`. This will reveal a second composition input column, as shown below. At this time, it will also set `Temperature range:` and `Pressure range:` to `Disabled`, as simultaneous loops over composition and temperature and pressure are not currently supported.

![`thermoGUI` window with temperature loop input.](/doc/images/thermoGUI-mloop.png)

In the example calculation loop shown above, the composition will be varied linearly between pure Pd and pure Ru.

## Saving, editing, and re-running calculations
`thermoGUI` works by writing input files in the Thermochimica input script mode, then calling `bin/InputScriptMode` (for simple calculations as well as temperature and/or pressure looping calculations) or `bin/RunCalculationList` (for calculations looping over composition). The input script is saved to `inputs/pythonInput.ti`. If one wishes to save these calculations settings for later, simply rename this file so it won't be overwritten. It can then be edited and re-used following the instructions for the Thermochimica input script mode.