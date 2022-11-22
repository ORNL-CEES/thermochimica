# `phaseDiagramGUI`
This interface is used to generate **binary** phase diagrams. It includes the following functions:
- [set up](#calculation-setup), 
- [refinement](#diagram-refinement), 
- [labelling](#labels),
- [plot settings](#plot-settings),
- [figure export](#export-plot),
- [data inspection](#inspect),
- [diagram overlay](#overlaying-diagrams),
- [comparison with experimental data](#experimental-data), and
- [generation and use of macros](#macros).

The `phaseDiagramGUI` window is shown below.

![`phaseDiagramGUI` default window](/doc/images/phaseDiagramGUI-default.png)

## Calculation Setup
The `Element 1` and `Element 2` dropdown menus are used to select the endmembers of the binary phase diagram. In this GUI, these must be pure elements.

The input boxes labelled `Minimum Temperature` and `Maximum Temperature` are used to set the limits of the diagram. If left blank, the default values for are 300 and 1000. The dropdown menu below is used to select the `Temperature unit`. Similarly, the `Pressure` input box can be used to set the pressure for which the diagram is calculated, with a default value of 1, and unit selected from the `Pressure unit` dropdown.

The `Initial grid density` input box is used to select the number of points along each axis used to populate the phase diagram. As explained in the [Diagram Refinement](#diagram-refinement) section, the phase diagram can be progressively improved from the initial calculation, and therefore it is not important to modify this setting from its default value of 10.

Once `Run` is pressed, the settings described above will be loaded, the initial phase diagram produced, and most features will be unlocked. From this point forward, all further operations will be applied to this diagram, with the settings selected at this time. The only exception is the `Run` button itself, which will close any active diagrams and start fresh with the current settings.

An example configuration for the Pd - Ru phase diagram using the [`Kaye_NobleMetals.dat`](/data/Kaye_NobleMetals.dat) database is shown below, followed by the resulting initial phase diagram.

![`phaseDiagramGUI` example setup](/doc/images/phaseDiagramGUI-example-setup.png)

![`phaseDiagramGUI` example initial diagram](/doc/images/pd-ru-phaseDiagram-setup.png)

## Diagram Refinement
The initial phase diagram will likely be quite sparsely populated, with rough phase boundary edges and missing regions. Three functions are available to improve the phase diagram: (manual) `Refine`, `Auto Refine` and `Auto Smoothen`. Each performs new Thermochimica calculations with the aim of adding more points to the phase diagram, and thereby increasing the coverage of the phase diagram and smoothness of the phase boundaries.

### `Refine`
The first option is to manually refine the phase diagram by selecting regions in which to perform additional calculations. The window shown below is used for entry of the limits of this region. The input mirrors that used to set up the calculation, except different numbers of grid points along the two axes may be supplied.

![`phaseDiagramGUI` `Refine` window](/doc/images/phaseDiagramGUI-refine.png)

`Refine` has largely been superceded by the automated routines `Auto Refine` and `Auto Smoothen`. It is recommended for use when these routines fail to resolve some area or feature of a phase diagram.

An example refinement and the resulting improved phase diagram are shown below.

![`phaseDiagramGUI` example setup](/doc/images/phaseDiagramGUI-example-refine.png)

![`phaseDiagramGUI` example refined diagram](/doc/images/pd-ru-phaseDiagram-refine.png)

### `Auto Refine`
The suggested first step to improve a rough phase diagram is to use `Auto Refine`. This routine uses a geometric analysis of the existing regions of the phase diagram to determine where the diagram is missing phase information, and efficiently deploys calculations in those areas. The density of the refinement grid is set by an internal parameter, which is increased with each use of `Auto Refine`. Thus, successive calls to `Auto Refine` will continue to improve the resolution of the phase diagram.

Some phase regions may be skipped and the terminal output should indicate this.

The example phase diagram is shown after one call to `Auto Refine` below.

![`phaseDiagramGUI` example auto-refined diagram](/doc/images/pd-ru-phaseDiagram-autorefine.png)

### `Auto Smoothen`
`Auto Smoothen` is similar in operation to `Auto Refine`: it also automatically determines where calculations should be performed, and uses a ratcheting internal parameter to progressively increase the diagram resolution. However, whereas `Auto Refine` searches for regions of missing phase information to populate, `Auto Smoothen` increases the smoothness of phase boundaries by performing calculations within known two-phase regions.

The example phase diagram is shown after one call to `Auto Smoothen` below.

![`phaseDiagramGUI` example auto-smoothened diagram](/doc/images/pd-ru-phaseDiagram-autosmoothen.png)

`Auto Smoothen` also performs a second function, which is to detect overlapping phase regions. This is required in the somewhat uncommon case of a two-phase region existing at low temperature, disappears, and then reappears as temperature is increased. Consider the following Mo - Ru phase diagram, again using the [`Kaye_NobleMetals.dat`](/data/Kaye_NobleMetals.dat) database. After `Auto Refine`, it has multiple overlapping two-phase regions.

![`phaseDiagramGUI` example phase diagram with overlapping two-phase regions](/doc/images/mo-ru-phaseDiagram-overlap.png)

After calling `Auto Smoothen`, this issue is resolved as shown below.

![`phaseDiagramGUI` example phase diagram with resolved two-phase regions](/doc/images/mo-ru-phaseDiagram-resolved.png)

## Labels
Labels of phase regions can be added and removed manually, as well as automatically generated.

### `Add Label`
Pressing the `Add Label` button opens a very simple dialogue box, as shown below, which requests the location at which a label will be added. Pressing the `Add Label` button in this dialogue performs a new Thermochimica calculation at the specified point, adds phase boundary data corresponding to that point to the existing phase diagram, and adds text corresponding to the equilibrium phases determined at that point.

![`phaseDiagramGUI` `Add Label` window](/doc/images/phaseDiagramGUI-addlabel.png)

An example label point and the result on the phase diagram are shown below.

![`phaseDiagramGUI` example add label](/doc/images/phaseDiagramGUI-example-addlabel.png)

![`phaseDiagramGUI` example diagram with label](/doc/images/pd-ru-phaseDiagram-label.png)

### `Auto Label`
`Auto Label` attempts to create one label per phase region on the phase diagram. This routine is different from `Add Label` in that it does not run any new Thermochimica calculations, so no new phase region data is added to the diagram. Instead, the geometric analysis used for [`Auto Refine`](#auto-refine) is used to detect the regions. The algorithm also attempts to determine the midpoints of the regions, so that the labels are placed in convenient positions. 

As with `Auto Refine`, some phase regions may be skipped and the terminal output should indicate this.

The example phase diagram with automatically-generated labels is shown below.

![`phaseDiagramGUI` example diagram with automatic labels](/doc/images/pd-ru-phaseDiagram-autolabel.png)

### `Remove Label`
When there are labels on the current phase diagram, pressing `Remove Label` opens a dialogue box with a list of the labels. The label text as well as the concentration (of element 2) and temperature coordinates are displayed. Under the `Remove Label?` heading, a checkbox is present corresponding to each label. When `Remove Label(s)` is pressed, all the labels for which the box is checked are removed, and the dialogue window closed automatically.

An example dialogue window for the Pd - Ru phase diagram with automatically generated labels is shown below.

![`phaseDiagramGUI` example remove label](/doc/images/phaseDiagramGUI-example-removelabel.png)

## `Plot Settings`

## `Export Plot`

## `Inspect`

## Overlaying Diagrams

## Experimental Data

## Macros