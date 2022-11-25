# `plotGUI`
The `plotGUI` is used to make figures displaying data calculated by Thermochimica and stored in a JSON database. It can be launched by using the following command:

```bash
./scripts/plotGUI.sh
```

The `plotGUI` default window is shown below.

![`plotGUI` default window](/doc/images/plotGUI-default.png)

## Axis setup
Variables to be plotted on the x-axis, left-hand y-axis, and right-hand y-axis may be selected from the dropdown menus. For the x-axis, the following options are available:
- iteration
- temperature
- pressure

For either y-axis, the following options are available:
- temperature
- pressure
- moles
- mole fraction
- chemical potential
- driving force
- vapor pressure
- moles of element in phase
- mole fraction of phase by element
- mole fraction of element by phase
- mole fraction of endmembers
- moles of elements
- element potential
- integral Gibbs energy
- functional norm
- GEM iterations
- \# phases
- heat capacity
- enthalpy
- entropy

All axes allow for linear or log scale to be used.

Many of the y-axis options will open a new dialogue window when selected. An example of such a window is shown below. Items in the box on the left-hand side can be selected and added to the right-hand side. Items on the right-hand side will be included in the plot. Multiple selections can be made by using `Ctrl` or `Shift` + click. 

The dropdown menu is used to select the phase for which options are displayed in the left-hand box. The `Add All` button will add all options for **all phases**, not just the currently selected phase.

![`plotGUI` selection window](/doc/images/plotGUI-selection.png)

## Generate the plot
Once variables have been selected for each desired axis, pressing the `Plot` button generates a figure using [`matplotlib`](https://matplotlib.org/).

## Plot settings
The `Plot Settings` button may be pressed at any time to configure some options for the displayed plot as well as for exported figures. The settings window is shown below.

![`plotGUI` settings window](/doc/images/plotGUI-settings.png)

For both y-axes, the following settings are available:
- `Marker Style`: 
    - `Lines`: display only lines connecting data points
    - `Points`: display only closed circles (•) at data points
    - `Both`: display both lines and closed circles
- `Plot Colors`:
    - `Colorful`: the color for each series is determined by taking equal intervals in the `rainbow` colorspace
    - `Black`: all lines/points are set to black
- `Show Legend`: selects whether series on this axis are included in the plot legend

The following settings are available to configure exported figures:
- `Export Filename`: name of exported figure (default is `plot`), which will be saved to the `outputs` directory (note that the extension corresponding to the format selected in `Export Format` will be appended)
- `Export Format`: image file format may be selected from the following:
    - `png`: portable network graphic
    - `pdf`: portable document format
    - `ps`: postscript image
    - `eps`: encapsulated postscript
    - `svg`: scalable vector graphic
- `Export DPI`: set resolution in dots per inch (default is 300)

## Saving a plot
Press `Export Plot` to save the current figure to the location specified using `Export Filename` in the settings menu.

## Exporting a plot script
Pressing `Export Plot Script` will save a `python` script to `python/generatedPlotScript.py`. This script can be used to regenerate the current figure (including updates if the JSON database is updated), as well as to customize the setup of that figure.