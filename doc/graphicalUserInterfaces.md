# Graphical User Interfaces (GUIs)
While the primary purpose of Thermochimica is to be used within multi-physics applications, stand-alone calculations are also supported. Four GUIs are provided to assist with setting up, running, and plotting data from stand-alone Thermochimica calculations. These are: 

- [`thermoGUI`](/doc/thermoGUI.md), which performs simple calculations or loops over calculations, 
- [`plotGUI`](/doc/plotGUI), which can generate plots of many parameters calculated by Thermochimica
- [`phaseDiagramGUI`](/doc/phaseDiagramGUI), which creates binary phase diagrams, and
- [`pseudoBinaryPhaseDiagramGUI`](/doc/pseudoBinaryPhaseDiagramGUI), which creates phase diagrams in which the endpoints are compounds rather than individual elements.

Each GUI is documented below. Installation instructions are available in the main [Thermochimica readme](/README.md#method-3-guis). A shell script to start each GUI is available in the `scripts` directory. From the root Thermochimica directory, simply run (for example):
```bash
./scripts/thermoGUI.sh
```
These scripts (except for `plotGUI.sh`) build Thermochimica before launching the selected GUI, so may take a few seconds to start when used for the first time or after changes to the source code.

`thermoGUI`, `phaseDiagramGUI`, and `pseudoBinaryPhaseDiagramGUI` begin by opening a database selection window. In this interface, all files with extension `.dat` (case insensitive), are presented. `plotGUI` opens a similar window, but `.json` files are displayed and the default directory is `outputs`. The user can select one of the available databases by clicking on it, or use the `Browse` button to change directories. The default directory is `data`. Selecting a file will open a calculation window. The details of these windows will be described in the corresponding sections below. Multiple calculation windows of a given type can be active at a time.

![Database selection window](/doc/images/databaseSelection.png)
