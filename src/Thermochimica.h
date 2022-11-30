#include <string>

// Extern declarations for the bind(C) fortran TCAPI_* subroutines

extern "C"
{
  void TCAPI_setThermoFilename(const char *, std::size_t);
  void TCAPI_setUnitTemperature(const char *, std::size_t);
  void TCAPI_setUnitPressure(const char *, std::size_t);
  void TCAPI_setUnitMass(const char *, std::size_t);
  void TCAPI_setStandardUnits();
  void TCAPI_setModelicaUnits();

  void TCAPI_getMolFraction(int *, double *, int *);
  void TCAPI_getChemicalPotential(int *, double *, int *);
  void TCAPI_getElementPotential(int *, double *, int *);

  void TCAPI_getOutputChemPot(const char *, std::size_t, double *, int *);
  void TCAPI_getOutputSolnSpecies(const char *, std::size_t, const char *, std::size_t, double *, double *, int *);
  void TCAPI_getOutputMolSpecies(const char *, std::size_t, double *, double *, int *);
  void TCAPI_getOutputMolSpeciesPhase(const char *, std::size_t, const char *, std::size_t, double *, int *);
  void TCAPI_getElementMolesInPhase(const char *, std::size_t, const char *, std::size_t, double *, int *);
  void TCAPI_getElementMoleFractionInPhase(const char *, std::size_t, const char *, std::size_t, double *, int *);

  void TCAPI_getSolnPhaseMol(const char *, std::size_t, double *, int *);
  void TCAPI_getPureConPhaseMol(const char *, std::size_t, double *, int *);
  void TCAPI_getPhaseIndex(const char *, std::size_t, int *, int *);

  void TCAPI_getOutputSiteFraction(const char *, std::size_t, int *, int *, double *, int *);
  void TCAPI_getSublSiteMol(const char *, std::size_t, int *, int *, double *, int *);

  void TCAPI_setPrintResultsMode(int *);

  void TCAPI_setElementMass(int *, double *);
  void TCAPI_presetElementMass(int *, double *);
  void TCAPI_setTemperaturePressure(double *, double *);
  void TCAPI_checkInfoThermo(int *);

  void TCAPI_sSParseCSDataFile();
  void TCAPI_thermochimica();

  void TCAPI_solPhaseParse(int *, double *);

  void TCAPI_thermoDebug();

  void TCAPI_printResults();
  void TCAPI_printState();

  void TCAPI_resetInfoThermo();
  void TCAPI_resetThermo();
  void TCAPI_resetThermoAll();

  // re-initialization-related functions
  void TCAPI_saveReinitData();
  void TCAPI_getReinitDataSizes(int *, int *);
  void TCAPI_getMolesPhase(double *);
  void TCAPI_getAssemblage(int *);
  void TCAPI_setReinitRequested(int *);
  void TCAPI_resetReinit();
  void TCAPI_getAllElementPotential(double *);
  void TCAPI_getElementFraction(int *, double *);

  // MQMQA functions
  void TCAPI_getMqmqaMolesPairs(const char *, std::size_t, double *, int *);
  void TCAPI_getMqmqaPairMolFraction(const char *, std::size_t, const char *, std::size_t, double *, int *);
  void TCAPI_getMqmqaNumberPairsQuads(const char *, std::size_t, int *, int *, int *);
  void TCAPI_getMqmqaConstituentFraction(const char *, std::size_t, int *, const char *, std::size_t, double *, int *);

  // MOOSE reinit functions
  void TCAPI_getReinitData(int *, double *, double *, double *, double *, int *, int *);
  void TCAPI_setReinitData(const int *, const int *, const int *, const double *, const double *, const double *, const double *, const int *);

  // Heat capacity, enthalpy, and entropy
  void TCAPI_setHeatCapacityEnthalpyEntropyRequested(int *);
  void TCAPI_getHeatCapacityEnthalpyEntropy(double *, double *, double *);

}
