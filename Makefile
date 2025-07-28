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
CC          = g++
FFPE_TRAPS  ?= zero
FCFLAGS = -Wall -O2 -ffree-line-length-none -fno-automatic -fbounds-check \
          -ffpe-trap=$(FFPE_TRAPS) -cpp \
          -D"DATA_DIRECTORY='$(DATA_DIR)'" \
          -D"DATA_DIRECTORY2='$(DATA_DIR2)'"
CCFLAGS     = -std=gnu++17

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    # links to lapack and blas libraries:
    LDLOC   = -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas -lgfortran
    # link flags for linux users:
    LDFLAGS = -O2 -fno-automatic -fbounds-check
endif
ifeq ($(UNAME_S),Darwin)
    # link flags for mac users:
    LDFLAGS = -O2 -framework Accelerate -fno-automatic -fbounds-check
endif
ifneq (,$(findstring NT,$(UNAME_S)))
    LDLOC   =  -llapack -lblas -lgfortran
    # link flags for Windows users:
    LDFLAGS = -O2 -fno-automatic -fbounds-check
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
DTSTT_DIR	= $(TST_DIR)/regression
DTSVAL_DIR	= $(TST_DIR)/validation
SHARED_DIR  = $(SRC_DIR)
SHARED_DIR += $(addprefix $(SRC_DIR)/,$(SRC_SDR))
CURR_DIR    = $(shell pwd)
DATA_DIR    = $(CURR_DIR)/data/
CSVT_DIR 	= $(TST_DIR)/csv/ 
VPATH		= $(SHARED_DIR)

# Separate modules and non-modules
modfiles := $(shell find src -name "Module*.f90")
srcfiles := $(shell find src -name "[^(Module)]*.f90")

## ========
## MODULES:
## ========
MODS_SRC    = $(patsubst %.f90, %.o, $(notdir $(modfiles)))
MODS_LNK    = $(addprefix $(OBJ_DIR)/,$(MODS_SRC))

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
C_LIB  		= $(OBJ_DIR)/$(TC-C_LIB)

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

## ===============
## REGRESSION TESTS
## ===============
RTST_SRC   = $(notdir $(wildcard $(DTSTT_DIR)/*.F90))
RTST_OBJ   = $(RTST_SRC:.F90=.o)
RTST_LNK   = $(addprefix $(OBJ_DIR)/,$(RTST_OBJ))
RTST_BIN   = $(addprefix $(BIN_DIR)/,$(basename $(RTST_SRC)))

## ===============
## VALIDATION TESTS
## ===============
VTST_SRC   = $(notdir $(wildcard $(DTSVAL_DIR)/*.F90))
VTST_OBJ   = $(VTST_SRC:.F90=.o)
VTST_LNK   = $(addprefix $(OBJ_DIR)/,$(VTST_OBJ))
VTST_BIN   = $(addprefix $(BIN_DIR)/,$(basename $(VTST_SRC)))

## =======
## COMPILE
## =======
all:  directories $(MODS_LNK) $(SHARED_LNK) $(SHARED_LIB) $(EXEC_LNK) $(EXE_BIN) $(C_LNK) $(C_LIB)

directories: ${OBJ_DIR} ${BIN_DIR}

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

# Enforce module dependency rules

$(srcfiles): $(MODS_LNK)

$(OBJ_DIR)/ModuleTesting.o: $(OBJ_DIR)/ModuleThermo.o
$(OBJ_DIR)/ModuleTesting.o: $(OBJ_DIR)/ModuleThermoIO.o

%.o: %.f90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: %.f90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: %.F90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(TST_DIR)/%.F90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(EXE_DIR)/%.F90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(SHARED_LIB): $(SHARED_LNK)
	$(AR) rcs $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.C
	$(CC) $(CCFLAGS) -c $< -o $@

$(C_LIB): $(C_LNK)
	$(AR) rcs $@ $^

$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(SHARED_LNK)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) $(LDFLAGS) -o $(BIN_DIR)/$* $< $(SHARED_LNK) $(LDLOC)

.PHONY: clean veryclean test doc cleandoc directories

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

## ===========
## REGRESSION
## ===========
regressiontest: $(RTST_LNK) $(SHARED_LNK) $(MODS_LNK) $(RTST_BIN)

$(RTST_LNK): $(OBJ_DIR)/%.o: $(DTSTT_DIR)/%.F90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(RTST_BIN): $(BIN_DIR)/%: $(OBJ_DIR)/%.o $(SHARED_LNK)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) $(LDFLAGS) -o $@ $< $(SHARED_LNK) $(LDLOC)

## ===========
## VALIDATION
## ===========
validationtest: $(VTST_LNK) $(SHARED_LNK) $(MODS_LNK) $(VTST_BIN)

$(VTST_LNK): $(OBJ_DIR)/%.o: $(DTSVAL_DIR)/%.F90
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@

$(VTST_BIN): $(BIN_DIR)/%: $(OBJ_DIR)/%.o $(SHARED_LNK)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) $(LDFLAGS) -o $@ $< $(SHARED_LNK) $(LDLOC)

## ===========
## ALL TESTS:
## ===========
test: all dailytest regressiontest validationtest

## ===========
## DEBUG:
## ===========
setdebug:
	$(eval FCFLAGS = -Wall -O0 -g -fno-automatic -fbounds-check -ffpe-trap=$(FFPE_TRAPS) -D"DATA_DIRECTORY='$(DATA_DIR)'")

debug: setdebug all dailytest
