# `pseudoBinaryPhaseDiagramGUI`
This interface is used to create phase diagrams in which the endmembers are compounds rather than pure elements. It can be launched by using the following command:

```bash
./scripts/pseudoBinaryPhaseDiagramGUI.sh
```

The calculation window, shown below, is similar to that used for the [`phaseDiagramGUI`](/doc/phaseDiagramGUI.md), except for the entry of the endmembers. The functions that have been implemented for this interface operate in the same way as for the `phaseDiagramGUI`, and users are directed to the documentation for that interface.

The `Use fuzzy stoichiometry` option is available. This setting is often helpful for pseudo-binary phase diagrams, which tend to involve under-determined systems. If the output phase diagram does not look reasonable without this setting (i.e. has overlapping phase regions), then try enabling it.

![`pseudoBinaryPhaseDiagramGUI` default window](/doc/images/pseudoBinaryPhaseDiagramGUI-default.png)

## Calculation Setup
Endmember compounds are configured using the `Composition 1` and `Composition 2` columns. These function in the same manner as composition entry in [`thermoGUI`](/doc/thermoGUI.md). In the `pseudoBinaryPhaseDiagramGUI`, the number of initial steps along the temperature and composition axes can be specified separately, and both have a default value of 10.

Setup otherwides proceeds as described in the [`phaseDiagramGUI`](/doc/phaseDiagramGUI.md#calculation_setup) documentation.

An example NaCl - AlCl<sub>3</sub> pseudo-binary phase diagram is shown below. The [NaCl-AlCl3.dat](/data/NaCl-AlCl3.dat) database was used.

![Example NaCl - AlCl<sub>3</sub> generated using `pseudoBinaryPhaseDiagramGUI`](/doc/images/nacl-alcl3-phaseDiagram.png)

## `Plot Settings`
The default exported figure name is `thermochimicaPseudoBinaryPhaseDiagram`.