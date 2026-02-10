#####################################################################################
##										    ##
##	Makefile for Thermochimica						    ##
##										    ##
##	Base Makefile written by M. Rodriguez Rodriguez and M.H.A. Piro	            ##
##	Modified S. Simunovic
##	make		- compile shared libraries				    ##
##      make test 	- compile all tests					    ##
##	make dailytest 	- compile application and unit tests that run daily	    ##
##	make weeklytest - compile application and unit tests that run weekly	    ##
##	make doc	- compile dOxygen HTML and LaTex documents	            ##
##      make docHTML    - compile dOxygen HTML document only 			    ##
##	make doclatex	- compile dOxygen LaTeX document only	 	            ##
##	make clean	- clean object files					    ##
##      make cleandoc   - clean dOxygen files					    ##
##	make veryclean  - clean object, executable, module and document files.	    ##
##										    ##
######################################################################################

## ===================
## COMPILER VARIABLES:
## ===================
AR          = ar
# Check if mpif90 is available, otherwise use gfortran
MPIF90      := $(shell command -v mpif90 2> /dev/null)
ifdef MPIF90
    FC      = mpif90
else
    FC      = gfortran
endif
CC          = g++
FFPE_TRAPS  ?= zero
FCFLAGS     = -Wall -O2 -ffree-line-length-none -fbounds-check -ffpe-trap=$(FFPE_TRAPS) -cpp -D"DATA_DIRECTORY='$(DATA_DIR)'" -D"OUTPUT_DIRECTORY='$(OUTPUT_DIR)'"
CCFLAGS     = -std=gnu++17

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    # links to lapack and blas libraries:
    LDLOC   = -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas -lgfortran
    # link flags for linux users:
    LDFLAGS = -O2 -fbounds-check
endif
ifeq ($(UNAME_S),Darwin)
    # link flags for mac users:
    LDFLAGS = -O2 -framework Accelerate -fbounds-check
endif
ifneq (,$(findstring NT,$(UNAME_S)))
    LDLOC   =  -llapack -lblas -lgfortran
    # link flags for Windows users:
    LDFLAGS = -O2 -fbounds-check
endif

## ====================
## DIRECTORY VARIABLES:
## ====================
MKDIR_P     = mkdir -p
DOC_DIR     = doc
TEX_DIR     = $(DOC_DIR)/latex
BIN_DIR     = bin
OBJ_DIR     = obj
SRC_DIR     = src
SRC_SDR     = debug gem module parser postprocess reinit reset setup ctz api
EXE_DIR     = $(SRC_DIR)/exec
TST_DIR     = test
LIB_DIR     = lib
DTST_DIR    = $(TST_DIR)/daily
SHARED_DIR  = $(SRC_DIR)
SHARED_DIR += $(addprefix $(SRC_DIR)/,$(SRC_SDR))
CURR_DIR    = $(shell pwd)
DATA_DIR    = $(CURR_DIR)/data/
OUTPUT_DIR    = $(CURR_DIR)/outputs/
VPATH		= $(SHARED_DIR)

# Separate modules and non-modules
modfiles := $(shell find src -name "Module*.f90")
srcfiles := $(shell find src -iname "*.f90" -and -not -name "Module*")

##
OBJ_FILES			=  $(addprefix $(OBJ_DIR)/,$(patsubst %.f90, %.o, $(patsubst %.F90, %.o, $(notdir $(srcfiles)))))

## ========
## MODULES:
## ========
MODS_OBJ    = $(patsubst %.f90, %.o, $(notdir $(modfiles)))
MODS_LNK    = $(addprefix $(OBJ_DIR)/,$(MODS_OBJ))

## =================
## LIBRARIES:
## =================
TC_LIB      = libthermochimica.a
SHARED_SRC  = $(foreach dir,$(SHARED_DIR),$(notdir $(wildcard $(dir)/*.f90)))
SHARED_SRCF = $(foreach dir,$(SHARED_DIR),$(notdir $(wildcard $(dir)/*.F90)))
SHARED_OBJ  = $(SHARED_SRC:.f90=.o)
SHARED_OBJ += $(SHARED_SRCF:.F90=.o)
SHARED_LNK  = $(addprefix $(OBJ_DIR)/,$(SHARED_OBJ))
SHARED_LIB  = $(OBJ_DIR)/$(TC_LIB)

## =================
## C interface library:
## =================
C_SRC       = Thermochimica-c.C Thermochimica-cxx.C
C_OBJ       = $(C_SRC:.C=.o)
C_LNK       = $(addprefix $(OBJ_DIR)/,$(C_OBJ))
TC-C_LIB    = libthermoc.a
C_LIB  			= $(OBJ_DIR)/$(TC-C_LIB)

## ============
## OLD EXECUTABLES:
## ============
EXEC_SRC    = $(notdir $(wildcard $(TST_DIR)/*.F90))
EXEC_SRC   += $(notdir $(wildcard $(EXE_DIR)/*.F90))
EXEC_OBJ    = $(EXEC_SRC:.F90=.o)
EXEC_LNK    = $(addprefix $(OBJ_DIR)/,$(EXEC_OBJ))
EXE_OBJ     = $(basename $(EXEC_SRC))
EXE_BIN     = $(addprefix $(BIN_DIR)/,$(EXE_OBJ))

## ============
## DAILY TESTS:
## ============
DTEST_SRC   = $(notdir $(wildcard $(DTST_DIR)/*.F90))
DTEST_OBJ   = $(DTEST_SRC:.F90=.o)
DTEST_LNK   = $(addprefix $(OBJ_DIR)/,$(DTEST_OBJ))
DTST_OBJ    = $(basename $(DTEST_SRC))
DTST_BIN    = $(addprefix $(BIN_DIR)/,$(DTST_OBJ))

## ====================
## TEST-DRIVE TESTS:
## ====================
TESTDRIVE_SRC     = testdrive.F90
TESTDRIVE_OBJ     = $(OBJ_DIR)/$(TESTDRIVE_SRC:.F90=.o)
TESTSUITE_SRC     = test_error_handling.F90 test_systems.F90
TESTSUITE_OBJ     = $(addprefix $(OBJ_DIR)/,$(TESTSUITE_SRC:.F90=.o))
TESTMAIN_SRC      = test_main.F90
TESTMAIN_OBJ      = $(OBJ_DIR)/$(TESTMAIN_SRC:.F90=.o)
TESTMAIN_BIN      = $(BIN_DIR)/test_main

## =======
## COMPILE
## =======
all:  directories $(MODS_LNK) $(SHARED_LNK) $(SHARED_LIB) $(EXEC_LNK) $(EXE_BIN) $(C_LNK) $(C_LIB)
	@echo "Compilation complete using $(FC)"

directories: ${OBJ_DIR} ${BIN_DIR}
	@echo "Using Fortran compiler: $(FC)"

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

# Enforce module dependency rules
$(OBJ_FILES): $(MODS_LNK)
$(EXEC_LNK) $(DTST_LNK): $(MODS_LNK)

$(OBJ_DIR)/%.o: %.f90 | $(OBJ_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: %.F90 | $(OBJ_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(TST_DIR)/%.F90 | $(OBJ_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(EXE_DIR)/%.F90 | $(OBJ_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(SHARED_LIB): $(SHARED_LNK)
	$(AR) rcs $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.C | $(OBJ_DIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(C_LIB): $(C_LNK)
	$(AR) rcs $@ $^

$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(SHARED_LNK) | $(BIN_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) $(LDFLAGS) -o $(BIN_DIR)/$* $< $(SHARED_LNK) $(LDLOC)

.PHONY: clean veryclean test doc cleandoc directories testdrive runtests

## =====
## CLEAN
## =====
clean:
	rm -f $(OBJ_DIR)/*
	find bin -name \*.dSYM -exec rm -rf {} \; > /dev/null 2>&1 | :
	rm -f $(BIN_DIR)/*

cleanexternal:
	rm -f $(SRC_DIR)/*.lo
	rm -f $(SRC_DIR)/*.lo.d
	rm -f $(SRC_DIR)/.libs/*
	rm -f $(SRC_DIR)/*.mod
	rm -f $(LIB_DIR)/*
	rm -f $(LIB_DIR)/.libs/*

veryclean: clean cleandoc cleanexternal
	rm -fr $(OBJ_DIR)/*
	rm -fr $(BIN_DIR)/*
	rm -f *.mod

## =======
## INSTALL
## =======
ifeq ($(PREFIX),)
  PREFIX := /usr/local
endif

install: $(SHARED_LIB)
	install -d $(DESTDIR)$(PREFIX)/lib/
	install -m 644 $(SHARED_LIB) $(DESTDIR)$(PREFIX)/lib/
	install -d $(DESTDIR)$(PREFIX)/include/
	install -m 644 $(OBJ_DIR)/*.mod $(DESTDIR)$(PREFIX)/include/

c-thermo: $(C_LIB)
	install -d $(DESTDIR)$(PREFIX)/lib/
	install -m 644 $(C_LIB) $(DESTDIR)$(PREFIX)/lib/

libraries: install c-thermo

## =============
## DOCUMENTATION
## =============
doc: dochtml doclatex

dochtml:
	doxygen Doxyfile

doclatex: dochtml
	cd $(TEX_DIR); make

doctest:
	cd $(TST_DIR); doxygen Doxyfile; cd $(TEX_DIR); make; cd ../..; mv $(DOC_DIR) ../$(DOC_DIR)/$(TST_DIR)

cleandoc:
	rm -r -f $(DOC_DIR)/html; rm -r -f $(TEX_DIR); rm -r -f $(TST_DIR)/$(DOC_DIR)/html; rm -r -f $(TST_DIR)/$(TEX_DIR); rm -r -f $(DOC_DIR)/$(TST_DIR)

## ===========
## DAILY TESTS
## ===========
dailytest: $(DTEST_LNK) $(SHARED_LNK) $(MODS_LNK) $(DTST_BIN)

$(OBJ_DIR)/%.o: $(DTST_DIR)/%.F90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

## =====================
## TEST-DRIVE FRAMEWORK
## =====================
# Compile testdrive framework
$(TESTDRIVE_OBJ): $(TST_DIR)/$(TESTDRIVE_SRC) | $(OBJ_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

# Compile test suite modules (depend on testdrive and Thermochimica modules)
$(OBJ_DIR)/test_error_handling.o: $(TST_DIR)/test_error_handling.F90 $(TESTDRIVE_OBJ) $(MODS_LNK) | $(OBJ_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/test_systems.o: $(TST_DIR)/test_systems.F90 $(TESTDRIVE_OBJ) $(MODS_LNK) | $(OBJ_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

# Compile test main (depends on test suites and testdrive)
$(TESTMAIN_OBJ): $(TST_DIR)/$(TESTMAIN_SRC) $(TESTDRIVE_OBJ) $(TESTSUITE_OBJ) | $(OBJ_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

# Link test main executable
$(TESTMAIN_BIN): $(TESTMAIN_OBJ) $(TESTSUITE_OBJ) $(TESTDRIVE_OBJ) $(SHARED_LNK) | $(BIN_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) $(LDFLAGS) -o $@ $(TESTMAIN_OBJ) $(TESTSUITE_OBJ) $(TESTDRIVE_OBJ) $(SHARED_LNK) $(LDLOC)

# Target to build test-drive test suite
testdrive: $(MODS_LNK) $(SHARED_LNK) $(TESTMAIN_BIN)
	@echo "Test-drive test suite built successfully: $(TESTMAIN_BIN)"

# Target to build and run test-drive tests
runtests: testdrive
	@echo ""
	@echo "Running test-drive test suite..."
	@$(TESTMAIN_BIN)

## ===========
## ALL TESTS:
## ===========
test: all dailytest testdrive

## ===========
## DEBUG:
## ===========
setdebug:
	$(eval FCFLAGS = -Wall -O0 -g -fbounds-check -ffpe-trap=$(FFPE_TRAPS) -D"DATA_DIRECTORY='$(DATA_DIR)'")

debug: setdebug all dailytest
