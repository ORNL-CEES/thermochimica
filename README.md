THERMOCHIMICA [![Basic Status](https://github.com/ORNL-CEES/thermochimica/actions/workflows/quick.yml/badge.svg)](https://github.com/ORNL-CEES/thermochimica/actions/workflows/quick.yml) [![Extended Status](https://github.com/ORNL-CEES/thermochimica/actions/workflows/main.yml/badge.svg)](https://github.com/ORNL-CEES/thermochimica/actions/workflows/main.yml)
=============

THERMOCHIMICA is a computational library for chemical thermodynamics. It determines a unique combination of phases and their compositions for a prescribed chemical composition, temperature and pressure. The solver and the underlying thermodynamic models can be used to estimate the chemical state and various constitutive and transport properties necessary for modeling of materials and processes.

# Installation instructions
## General
Thermochimica requires a Fortran compiler and is regularly built using the GNU gfortran compiler. Also required are the LAPACK/BLAS libraries. If use of the Graphical User Interface (GUI) is desired, then python (tested with 3.8+) and some pip packages (listed in python/requirements.txt) are also required.

Suggested modes of installation and operation for Ubuntu, MacOS, and Windows Subsystem for Linux (WSL) are decribed below. Please contact us or open an issue on GitHub should you encounter any difficulties. Also, suggestions for improvements to these instructions are welcome, particularly for native Windows installation.

## Prerequisites
The prerequisites for cloning and building Thermochimica vary depending on the operating system. The instructions for some common OSs are provided below. 
### Ubuntu
Get packages required for basic Fortran compilation:
```bash
sudo apt install build-essential gfortran liblapack-dev git
```
Now follow the [steps for obtaining and building Thermochimica](#building-thermochimica).
### MacOS
If not already installed, install Xcode from the [Apple developer website](https://developer.apple.com/downloads/index.action) or using the [Mac App Store](https://apps.apple.com/us/app/xcode/id497799835).

With Xcode installed, you can install Apple's Command Line Developer Tools by running the following command in the terminal:
```bash
xcode-select --install
```
The other dependencies for Thermochimica can be installed and managed using either [Homebrew](https://brew.sh/) or [MacPorts](https://www.macports.org/index.php). Once either Brew or MacPorts is installed, they can be used to install the necessary packages. Using _Homebrew_, execute the following in terminal:
```bash
brew install make gcc
```
or using _MacPorts_, execute:
```bash
sudo port install gcc8
```
You can now proceed to [steps for obtaining and building Thermochimica](#building-thermochimica).
### Windows (WSL)
You can install the Windows Subsystem for Linux by following [the instructions from Microsoft](https://docs.microsoft.com/en-us/windows/wsl/install).

The default (Ubuntu) is good, so just run the one line there, restart, and try to start it. So in the windows shell (as administrator):
```bash
wsl --install
```
If WSL fails to start properly after reboot, you may need to create virtual disk space for WSL following [these instructions](https://utf9k.net/blog/wsl2-vhd-issue/).

When Ubuntu app is successfully installed, run the following.
```bash
sudo apt update && sudo apt upgrade
```
```bash
sudo apt install build-essential gfortran liblapack-dev git
```
Now you can follow the [Building Thermochimica](#building-thermochimica) instructions below.

## Building Thermochimica
Clone the repository and navigate to the root Thermochimica directory:
```bash
git clone https://github.com/ORNL-CEES/thermochimica.git
cd thermochimica
```
Build Thermochimica (with tests):
```bash
make test
```
Run tests:
```bash
./run_tests
```

# Operation
Thermochimica can be operated in three modes:
1. Writing and compiling Fortran driver files (like those in the `test` directory).
2. Using input scripts (like those in the `inputs` directory).
3. Using the GUIs (can be run via scripts in the `scripts` directory).

[Method 1](#method-1-fortran-driver-files) is the most involved, but also the most powerful, in the sense that there are no limitations on the series of calculations that can be run or output obtained. [Method 2](#method-2-input-scripts) and [Method 3](#method-3-guis) have features that allow for simple loops over temperature, pressure, or composition to be performed.

## Method 1: Fortran Driver Files
Create a new `.F90` file in the `test` directory, for example `demo.F90`. You may want to start by copying an existing test file as a template. For example, from the root Thermochimica directory:
```bash
cp test/Thermo.F90 test/demo.F90
```
Make whatever changes you like to `demo.F90`, and then when you are ready to try it, the next step is to compile.
```bash
make
```
Now you can run your calculation:
```bash
./bin/demo
```
If you wish to save the output as a JSON, add the following line of code after Thermochimica is called. 
```bash
call WriteJSON(.TRUE.)
```
The output will be saved as the `outputs\thermoout.json`.
## Method 2: Input Scripts
Again, we start by creating a new file, this time in the `inputs` directory, and can start by copying an existing file:
```bash
cp inputs/advanced-input.ti inputs/demo.ti
```
Take a look at `demo.ti`, and make some changes. Note the `pressure` and `temperature` lines use the format `start:stop:step` to set ranges of values to loop over. You may also want to add the following line to enable output of results in JSON format:
```bash
write json        = .TRUE.
```
The output will be saved as the `outputs\thermoout.json`.
When you are done, the script can be run:
```bash
./bin/InputScriptMode inputs/demo.ti
```

## Method 3: GUIs
The GUIs for Thermochimica depend on Python(3.8+) and some additional Python packages that can be installed via pip. For Ubuntu or WSL with Ubuntu, you can follow these instructions.

First get pip and tkinter for Python:
```bash
sudo apt install python3-pip python3-tk
```
Navigate to Thermochimica folder, now pip can be used to install the required packages listed in the `python/requirements.txt` file:
```bash
pip install -r python/requirements.txt
```

For WSL only, you will need X11 for the Ubuntu app to be able to open a window on your desktop. First, install [Xming](https://sourceforge.net/projects/xming/). Then (still in Windows), run XLaunch with "Display number" set to 0, and on the third page make sure to check "No Access Control". After this, there should be an instance of Xming in your system tray. Now open the Ubuntu app and add:
```bash
export DISPLAY=$(route.exe print | grep 0.0.0.0 | head -1 | awk '{print $4}'):0.0
```
to your .bashrc. This will be applied in any new terminal, so either restart your Ubuntu terminal, or just run the same command in your current one.

Finally, launch a GUI:
```bash
./scripts/thermoGUI.sh
```

If you have installed a python version other than your system default that you would like to use to run the Thermochimica GUIs, you can set the environment variable `python_for_thermochimica` to point to that python. Multiple interactive matplotlib windows seem to work better on python3.9 than with python3.8.

Further documentation on the use of the GUIs is available in the [GUI docs](/doc/graphicalUserInterfaces.md).

# License
Thermochimica has a [BSD 3-clause open-source license](LICENSE).
