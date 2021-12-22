{% include styles.html %}
# Installation instructions
## General

Thermochimica requires a Fortran compiler and is regularly built using the GNU gfortran compiler. Also required are the LAPACK/BLAS libraries. If use of the Graphical User Interface (GUI) is desired, then python (tested with 3.8+) and some pip packages (listed in python/requirements.txt) are also required.

Suggested modes of installation and operation for Ubuntu, MacOS, and Windows Subsystem for Linux (WSL) are decribed below. Please contact us or open an issue on GitHub should you encounter any difficulties. Also, suggestions for improvements to these instructions are welcome, particularly for native Windows installation.

## Ubuntu

Get packages required for basic Fortran compilation:
{% include codeHeader.html %}
```bash
sudo apt install build-essential gfortran liblapack-dev git
```

Clone the repository and navigate to the root Thermochimica directory:
{% include codeHeader.html %}
```bash
git clone https://github.com/ORNL-CEES/thermochimica.git
cd thermochimica
```

Build Thermochimica (with tests):
{% include codeHeader.html %}
```bash
make test
```

Run tests:
{% include codeHeader.html %}
```bash
./run_tests
```

## MacOS
You can install [Brew](https://brew.sh/) to manage and install dependencies.

Once Brew is installed, it can be used to install the necessary packages:
{% include codeHeader.html %}
```bash
brew install make gcc git 
```

## Windows (WSL)
You can install the Windows Subsystem for Linux by following [the instructions from Microsoft](https://docs.microsoft.com/en-us/windows/wsl/install).

The default (Ubuntu) is good, so just run the one line there, restart, and try to start it. So in the windows shell (as administrator):
{% include codeHeader.html %}
```bash
wsl --install
```

If WSL fails to start properly after reboot, you may need to create virtual disk space for WSL following [these instructions](https://utf9k.net/blog/wsl2-vhd-issue/).

Once the Ubuntu app is successfully installed, you can follow the Ubuntu instructions above.

# Operation
Thermochimica can be operated in three modes:
1. Writing and compiling Fortran driver files (like those in the `test` directory).
2. Using input scripts (like those in the `inputs` directory).
3. Using the GUIs (can be run via scripts in the `scripts` directory).

Method 1 is the most involved, but also the most powerful, in the sense that there are no limitations on the series of calculations that can be run or output obtained. Methods 2 and 3 have features that allow for simple loops over temperature, pressure, or composition to be performed.

## Method 1: Fortran Driver Files
Create a new `.F90` file in the `test` directory, for example `demo.F90`. You may want to start by copying an existing test file as a template. For example, from the root Thermochimica directory:
{% include codeHeader.html %}
```bash
cp test/Thermo.F90 test/demo.F90
```

Make whatever changes you like to `demo.F90`, and then when you are ready to try it, the next step is to compile.
{% include codeHeader.html %}
```bash
make
```

Now you can run your calculation:
{% include codeHeader.html %}
```bash
./bin/demo
```

## Method 2: Input Scripts
Again, we start by creating a new file, this time in the `inputs` directory, and can start by copying an existing file:
{% include codeHeader.html %}
```bash
cp inputs/advanced-input.ti inputs/demo.ti
```

Take a look at `demo.ti`, and make some changes. Note the `pressure` and `temperature` lines use the format start:stop:step to set ranges of values to loop over. You may also want to add the following line to enable output of results in JSON format:
{% include codeHeader.html %}
```bash
write json        = .TRUE.
```

When you are done, the script can be run:
{% include codeHeader.html %}
```bash
./bin/InputScriptMode inputs/demo.ti
```

## Method 3: GUIs
The GUIs for Thermochimica depend on Python(3.8+) and some additional Python packages that can be installed via pip. For Ubuntu or WSL with Ubuntu, you can follow these instructions.

First get pip and tkinter for Python:
{% include codeHeader.html %}
```bash
sudo apt install python3-pip python3-tk
```

Now pip can be used to install the required packages:
{% include codeHeader.html %}
```bash
pip install -r python/requirements.txt
```

For WSL only, you will need X11 for the Ubuntu app to be able to open a window on your desktop. First, install [Xming](https://sourceforge.net/projects/xming/). Then (still in Windows), run XLaunch with "Display number" set to 0, and on the third page make sure to check "No Access Control". After this, there should be an instance of Xming in your system tray. Now open the Ubuntu app and add:
{% include codeHeader.html %}
```bash
export DISPLAY=$(route.exe print | grep 0.0.0.0 | head -1 | awk '{print $4}'):0.0
```
to your .bashrc. This will be applied in any new terminal, so either restart your Ubuntu terminal, or just run the same command in your current one.

Finally, launch a GUI:
{% include codeHeader.html %}
```bash
python3 python/thermoGui.py
```

{% include copyCode.html %}
