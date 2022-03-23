#!/usr/bin/env bash

make > make.out

source scripts/setPython.sh

$python_for_thermochimica python/justPdGui.py
