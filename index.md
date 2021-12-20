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


## Windows (WSL)


{% include copyCode.html %}
