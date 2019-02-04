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
FC          = gfortran
FCFLAGS     = -Wall -g -O0 -fno-automatic -fbounds-check -ffpe-trap=zero -D"DATA_DIRECTORY='$(DATA_DIR)'"
#FCFLAGS     = -Wall -g -fbounds-check
#FCFLAGS     = -Wall -g -O0 -fno-automatic -fbounds-check
#LDFLAGS     = -framework Accelerate -g -fbounds-check
#LDFLAGS     = -O0 -framework Accelerate -g -fno-automatic -fbounds-check
#LDFLAGS     =  -O0 -g -fno-automatic -fbounds-check

# links to lapack and blas libraries:
LDLOC     =  -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas -lgfortran

# link flags for linux users:
LDFLAGS     =  -O0 -g -fno-automatic -fbounds-check

# link flags for mac users:
#LDFLAGS     = -O0 -framework Accelerate -g -fno-automatic -fbounds-check


## ====================
## DIRECTORY VARIABLES:
## ====================

MKDIR_P     = mkdir -p

DOC_DIR     = doc
TEX_DIR     = $(DOC_DIR)/latex
BIN_DIR     = bin
OBJ_DIR     = obj
SRC_DIR     = src
TST_DIR     = test
DTST_DIR    = $(TST_DIR)/daily
WTST_DIR    = $(TST_DIR)/weekly
SHARED_DIR  = $(SRC_DIR)

CURR_DIR    = $(shell pwd)
DATA_DIR    = $(CURR_DIR)/data/

## ========
## MODULES:
## ========

MODS_SRC    = ModuleThermo.o ModuleThermoIO.o ModuleGEMSolver.o ModuleSubMin.o ModuleParseCS.o ModuleSS.o
MODS_LNK    = $(addprefix $(OBJ_DIR)/,$(MODS_SRC))


## =================
## SHARED LIBRARIES:
## =================

TC_LIB      = libthermochimica.a
SHARED_SRC  = $(foreach dir,$(SHARED_DIR),$(notdir $(wildcard $(dir)/*.f90)))
SHARED_OBJ  = $(SHARED_SRC:.f90=.o)
SHARED_LNK  = $(addprefix $(OBJ_DIR)/,$(SHARED_OBJ))
SHARED_LIB  = $(OBJ_DIR)/$(TC_LIB)

## ============
## OLD EXECUTABLES:
## ============

EXEC_SRC    = $(notdir $(wildcard $(TST_DIR)/*.F90))
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


## =============
## WEEKLY TESTS:
## =============

WTEST_SRC   = $(notdir $(wildcard $(WTST_DIR)/*.F90))
WTEST_OBJ   = $(WTEST_SRC:.F90=.o)
WTEST_LNK   = $(addprefix $(OBJ_DIR)/,$(WTEST_OBJ))

WTST_OBJ    = $(basename $(WTEST_SRC))
WTST_BIN    = $(addprefix $(BIN_DIR)/,$(WTST_OBJ))

## =======
## COMPILE
## =======

all:  directories $(MODS_LNK) $(SHARED_LNK) $(SHARED_LIB) $(EXEC_LNK) $(EXE_BIN)

directories: ${OBJ_DIR} ${BIN_DIR}

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(TST_DIR)/%.F90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(SHARED_LIB): $(SHARED_LNK)
	$(AR) rcs $@ $^

$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(SHARED_LNK)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) $(LDFLAGS) -o $(BIN_DIR)/$* $< $(SHARED_LNK) $(LDLOC)

.PHONY: clean veryclean test doc cleandoc directories


## =====
## CLEAN
## =====

clean:
	rm -f $(OBJ_DIR)/*
	rm -f $(BIN_DIR)/*

veryclean: clean cleandoc
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

## ============
## WEEKLY TESTS
## ============

weeklytest: $(WTEST_LNK) $(SHARED_LNK) $(MODS_LNK) $(WTST_BIN)

$(OBJ_DIR)/%.o: $(WTST_DIR)/%.F90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

## ===========
## BOTH TESTS:
## ===========

test: dailytest weeklytest
