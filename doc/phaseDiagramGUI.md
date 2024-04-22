# `binaryPhaseDiagramGUI`

This interface is used to generate **binary** phase diagrams. It includes the following functions:

- [set up](#calculation-setup),
- [refinement](#diagram-refinement),
- [labelling](#labels),
- [undo](#undo),
- [plot settings](#plot-settings),
- [figure export](#export-plot),
- [data inspection](#inspect),
- [diagram overlay](#overlaying-diagrams),
- [comparison with experimental data](#experimental-data), and
- [generation and use of macros](#macros).

It can be launched by using the following command:

```bash
./scripts/phaseDiagramGUI.sh
```

The `binaryPhaseDiagramGUI` window is shown below.

![`binaryPhaseDiagramGUI` default window](/doc/images/binaryPhaseDiagramGUI-default.png)

## Calculation Setup

The `Element 1` and `Element 2` dropdown menus are used to select the endmembers of the binary phase diagram. In this GUI, these must be pure elements.

The input boxes labelled `Minimum Temperature` and `Maximum Temperature` are used to set the limits of the diagram. If left blank, the default values for are 300 and 1000. The dropdown menu below is used to select the `Temperature unit`. Similarly, the `Pressure` input box can be used to set the pressure for which the diagram is calculated, with a default value of 1, and unit selected from the `Pressure unit` dropdown.

The `Initial grid density` input box is used to select the number of points along each axis used to populate the phase diagram. As explained in the [Diagram Refinement](#diagram-refinement) section, the phase diagram can be progressively improved from the initial calculation, and therefore it is not important to modify this setting from its default value of 10.

As in [`thermoGUI`](/doc/thermoGUI.md), the `Use fuzzy stoichiometry` option is available. This setting is typically **not recommended** for binary phase diagrams.

Once `Run` is pressed, the settings described above will be loaded, the initial phase diagram produced, and most features will be unlocked. From this point forward, all further operations will be applied to this diagram, with the settings selected at this time. The only exception is the `Run` button itself, which will close any active diagrams and start fresh with the current settings.

An example configuration for the Pd - Ru phase diagram using the [`Kaye_NobleMetals.dat`](/data/Kaye_NobleMetals.dat) database is shown below, followed by the resulting initial phase diagram.

![`binaryPhaseDiagramGUI` example setup](/doc/images/binaryPhaseDiagramGUI-example-setup.png)

![`binaryPhaseDiagramGUI` example initial diagram](/doc/images/pd-ru-phaseDiagram-setup.png)

## Diagram Refinement

The initial phase diagram will likely be quite sparsely populated, with rough phase boundary edges and missing regions. Three functions are available to improve the phase diagram: (manual) `Refine`, `Auto Refine` and `Auto Smoothen`. Each performs new Thermochimica calculations with the aim of adding more points to the phase diagram, and thereby increasing the coverage of the phase diagram and smoothness of the phase boundaries.

### `Refine`

The first option is to manually refine the phase diagram by selecting regions in which to perform additional calculations. The window shown below is used for entry of the limits of this region. The input mirrors that used to set up the calculation, except different numbers of grid points along the two axes may be supplied.

![`binaryPhaseDiagramGUI` `Refine` window](/doc/images/phaseDiagramGUI-refine.png)

`Refine` has largely been superceded by the automated routines `Auto Refine` and `Auto Smoothen`. It is recommended for use when these routines fail to resolve some area or feature of a phase diagram.

An example refinement and the resulting improved phase diagram are shown below.

![`binaryPhaseDiagramGUI` example setup](/doc/images/phaseDiagramGUI-example-refine.png)

![`binaryPhaseDiagramGUI` example refined diagram](/doc/images/pd-ru-phaseDiagram-refine.png)

### `Auto Refine`

The suggested first step to improve a rough phase diagram is to use `Auto Refine`. This routine uses a geometric analysis of the existing regions of the phase diagram to determine where the diagram is missing phase information, and efficiently deploys calculations in those areas. The density of the refinement grid is set by an internal parameter, which is increased with each use of `Auto Refine`. Thus, successive calls to `Auto Refine` will continue to improve the resolution of the phase diagram.

Some phase regions may be skipped and the terminal output should indicate this.

The example phase diagram is shown after one call to `Auto Refine` below.

![`binaryPhaseDiagramGUI` example auto-refined diagram](/doc/images/pd-ru-phaseDiagram-autorefine.png)

### `Auto Smoothen`

`Auto Smoothen` is similar in operation to `Auto Refine`: it also automatically determines where calculations should be performed, and uses a ratcheting internal parameter to progressively increase the diagram resolution. However, whereas `Auto Refine` searches for regions of missing phase information to populate, `Auto Smoothen` increases the smoothness of phase boundaries by performing calculations within known two-phase regions.

The example phase diagram is shown after one call to `Auto Smoothen` below.

![`binaryPhaseDiagramGUI` example auto-smoothened diagram](/doc/images/pd-ru-phaseDiagram-autosmoothen.png)

`Auto Smoothen` also performs a second function, which is to detect overlapping phase regions. This is required in the somewhat uncommon case of a two-phase region existing at low temperature, disappears, and then reappears as temperature is increased. Consider the following Mo - Ru phase diagram, again using the [`Kaye_NobleMetals.dat`](/data/Kaye_NobleMetals.dat) database. After `Auto Refine`, it has multiple overlapping two-phase regions.

![`binaryPhaseDiagramGUI` example phase diagram with overlapping two-phase regions](/doc/images/mo-ru-phaseDiagram-overlap.png)

After calling `Auto Smoothen`, this issue is resolved as shown below.

![`binaryPhaseDiagramGUI` example phase diagram with resolved two-phase regions](/doc/images/mo-ru-phaseDiagram-resolved.png)

This detection of overlapping regions is done by using heuristic analysis of the spacing between consecutive points on a phase boundary line. This analysis fails if phase boundaries are sparsely populated, which is why it is not performed until `Auto Smoothen` is called.

## Labels

Labels of phase regions can be added and removed manually, as well as automatically generated.

### `Add Label`

Pressing the `Add Label` button opens a very simple dialogue box, as shown below, which requests the location at which a label will be added. Pressing the `Add Label` button in this dialogue performs a new Thermochimica calculation at the specified point, adds phase boundary data corresponding to that point to the existing phase diagram, and adds text corresponding to the equilibrium phases determined at that point.

![`binaryPhaseDiagramGUI` `Add Label` window](/doc/images/phaseDiagramGUI-addlabel.png)

An example label point and the result on the phase diagram are shown below.

![`binaryPhaseDiagramGUI` example add label](/doc/images/phaseDiagramGUI-example-addlabel.png)

![`binaryPhaseDiagramGUI` example diagram with label](/doc/images/pd-ru-phaseDiagram-label.png)

### `Auto Label`

`Auto Label` attempts to create one label per phase region on the phase diagram. This routine is different from `Add Label` in that it does not run any new Thermochimica calculations, so no new phase region data is added to the diagram. Instead, the geometric analysis used for [`Auto Refine`](#auto-refine) is used to detect the regions. The algorithm also attempts to determine the midpoints of the regions, so that the labels are placed in convenient positions.

As with `Auto Refine`, some phase regions may be skipped and the terminal output should indicate this.

The example phase diagram with automatically-generated labels is shown below.

![`binaryPhaseDiagramGUI` example diagram with automatic labels](/doc/images/pd-ru-phaseDiagram-autolabel.png)

### `Remove Label`

When there are labels on the current phase diagram, pressing `Remove Label` opens a dialogue box with a list of the labels. The label text as well as the concentration (of element 2) and temperature coordinates are displayed. Under the `Remove Label?` heading, a checkbox is present corresponding to each label. When `Remove Label(s)` is pressed, all the labels for which the box is checked are removed, and the dialogue window closed automatically.

An example dialogue window for the Pd - Ru phase diagram with automatically generated labels is shown below.

![`binaryPhaseDiagramGUI` example remove label](/doc/images/phaseDiagramGUI-example-removelabel.png)

## `Undo`

The most recent `Refine`, `Auto Refine`, `Auto Smoothen`, `Label`, `Auto Label`, or `Remove Label` operation can be undone by pressing `Undo`. Note that only one level of history is stored at this time, so subsequent `Undo` operations are not possible.

## `Plot Settings`

When `Plot Settings` is pressed, the following settings window is opened.

![`binaryPhaseDiagramGUI` plot settings window](/doc/images/phaseDiagramGUI-plotsettings.png)

A few plot options are available:

- `Marker Style`:
  - `Lines`: display only lines connecting data points
  - `Points`: display only closed circles (•) at data points
  - `Both`: display both lines and closed circles
- `Plot Colors`:
  - `Colorful`: the color for the outline of each two-phase region is determined by taking equal intervals in the `rainbow` colorspace
  - `Black`: all lines are set to black
- `Experimental Data Colors`:
  - `Colorful`: the color for each loaded experimental data series is determined by taking equal intervals in the `rainbow` colorspace
  - `Black`: all loaded experimental data series are set to black
- `Show`:
  - `Experimental Data`: toggles whether loaded experimental data is displayed
  - `Loaded Diagram`: toggles whether a loaded phase diagram is displayed

`Auto-Label Settings` presents two toggles for the [`Auto Label`](#auto-label) feature, which determine whether single-phase regions and two-phase regions given labels when using `Auto Label`. Note that these only apply to subsequent uses of `Auto Label` and will not alter existing labels.

The following settings are available to configure exported figures:

- `Export Filename`: name of exported figure (default is `thermochimicaPhaseDiagram`), which will be saved to the `outputs` directory (note that the extension corresponding to the format selected in `Export Format` will be appended)
- `Export Format`: image file format may be selected from the following:
  - `png`: portable network graphic
  - `pdf`: portable document format
  - `ps`: postscript image
  - `eps`: encapsulated postscript
  - `svg`: scalable vector graphic
- `Export DPI`: set resolution in dots per inch (default is 300)

## `Export Plot`

Press `Export Plot` to save the current phase diagram to the location specified using `Export Filename` in the settings menu.

## `Inspect`

The `Inspect` feature allows a user to examine the underlying equilibrium calculations used to construct the phase diagram. This can be useful for debugging purposes, and also allows spurious calculations to be suppressed from a phase diagram. An example of the inspection window corresponding to the Pd - Ru phase diagram is shown below.

![`binaryPhaseDiagramGUI` `Inspect` window example for Pd - Ru system](/doc/images/phaseDiagramGUI-example-inspect.png)

The column on the left-hand side lists all equilibrium calculations used in the phase diagram. These are listed with an index, followed by the temperature and concentrations of `Element 2` in the two phases present. When selected, information and options for a calculation are displayed on the right-hand side of the `Inspect` window.

Under `Calculation Details`, the temperature and concentrations of the two elements are displayed, followed by the calculated compositions of the two phases determined to be present at equilibrium, and the Gibbs energy and total number of Gibbs energy minimization (GEM) iterations required for the calculation. For each calculation, it is the calculated compositions `Phase 1` and `Phase 2` that appear on the phase diagram.

The `Toggle Active/Suppressed Status` button sets whether the calculation is included in the phase diagram. The current status is indicated above this button (`Active` in the example figure above). By default, all calculations are set to `Active`.

Below `Filter Points` there are options for filtering the points listed in the left-hand column. This can be useful for locating a calculation of interest. The boxes below `Temperature Range` and `<Element 2> Concentration Range` (`Ru Concentration Range` in the example) are used to set a region of the phase diagram to consider. Note that the concentration used to determine inclusion is the location of the plotted points (i.e. `Phase 1` or `Phase 2` under `Calculation Details`), rather than the total concentration for the calculation (i.e. `Moles of <Element 2>` under `Calculation Details`).

There are two dropdowns below `Contains Phases`, and the various phases present in the phase diagram may be selected from these. Zero, one, or two phases may be selected. If two phases are selected, **both** phases must be present for a calculation to be listed.

There is also a dropdown to filter by active or suppressed status. When this menu is empty, points with either are shown.

An example of filtered data is shown below.

![`binaryPhaseDiagramGUI` `Inspect` window example with filter for Pd - Ru system](/doc/images/phaseDiagramGUI-example-inspect-filtered.png)

## Overlaying Diagrams

Once a phase diagram has been created, the data corresponding to that diagram can be saved by pressing `Export Diagram Data`. A dialogue window (shown below) will open to ask for a name for the file to which the data will be saved. The default name is `savedDiagram`. The extension is `pkl`, as the saved file is a [`pickle`](https://docs.python.org/3/library/pickle.html). The file will be saved to the `outputs` directory.

![`binaryPhaseDiagramGUI` `Export Diagram Data` window](/doc/images/phaseDiagramGUI-exportdiagramdata.png)

A diagram saved by this method can be loaded with the `Load Diagram` button. This opens a file selection window (opened by default in the `outputs` directory) that lists all `pkl` files. When a file is selected, the diagram saved in that file will be added as an overlay to the current phase diagram.

For example, consider a case in which the example Pd - Ru phase diagram was saved after the initial diagram setup to a file `pd-ru-example.pkl`. The `Load Diagram` window will appear as the below.

![`binaryPhaseDiagramGUI` example `Load Diagram` window](/doc/images/phaseDiagramGUI-example-loaddata.png)

Once `Auto Refine` and `Auto Smoothen` have been applied to the current phase diagram (as described in [Diagram Refinement](#diagram-refinement)), the resulting phase diagram with overlay will look like the figure below.

![`binaryPhaseDiagramGUI` example diagram with overlay](/doc/images/pd-ru-phaseDiagram-overlay.png)

## Experimental Data

Experimental data points may be added to the phase diagram for the purposes of comparing the computed diagram to experimental data. Experimental data is expected in a simple two-column CSV format, with concentration of `Element 2` in the first column, and temperature in the second column. [An example CSV file is available here](/doc/examples/pd-ru-example_experiment1.csv).

When the `Add Data` button is pressed, a file selection window is opened, which lists `csv` files. An example file selection window is shown below. As with the other file selection windows, `Browse` can be pressed to change directory. In the file selection column, multiple files can be selected simultaneously using `Ctrl` or `Shift` + click. When `Add Data` is pressed, the data in the selected files will be added to the phase diagram.

![`binaryPhaseDiagramGUI` example `Add Data` window](/doc/images/phaseDiagramGUI-example-experimental.png)

An example of the Pd - Ru phase diagram with two (fictitious) settings of experimental data added is shown below.

![`binaryPhaseDiagramGUI` example diagram with experiment](/doc/images/pd-ru-phaseDiagram-experiment.png)

## Macros

While a phase diagram is being generated, refined, and labelled, the button and settings used are automatically logged to a macro, i.e. a collection of commands. By pressing `Export Macro`, the logged commands can be saved to a file. This file will be placed in the `python` directory, and is actually a full-fledged python script, which means it can be run outside of the `binaryPhaseDiagramGUI` environment. Macros can also be edited to add or remove commands as necessary. This makes them an effective way to 'save' your work when using the `binaryPhaseDiagramGUI`.

Pressing `Macro Settings` opens a window (shown below) with a file browser and a field for entering `Macro File Save Name`. The file browser and file selection column are used to select a previously-saved macro. This selection is confirmed by pressing `Select Macro`. The `Macro File Save Name` input box is used to set the file name for the macro saved when `Export Macro` is pressed in the main `binaryPhaseDiagramGUI` window. This setting is confirmed by pressing `Set Save Name`. Note that the file will always have extension `py`, and the default name is `macroPhaseDiagram`.

![`binaryPhaseDiagramGUI` example `Macro Settings` window](/doc/images/phaseDiagramGUI-example-macrosettings.png)

When `Run Macro` is pressed, the macro selected in the `Macro Settings` window is run. Further `Refine`, `Label`, etc., operations will be applied to the diagram generated.

`Clear Macro` can be used to erase the commands in the currect macro log.
