#!/usr/bin/env bash

make -j > make.out

source scripts/setPython.sh

$python_for_thermochimica python/thermoGui.py &
