#!/usr/bin/env bash

make test > make.out

source scripts/setPython.sh

$python_for_thermochimica python/phaseDiagramGui.py
