#include <string>
#include <tuple>
#pragma once

#include <utility>
#include <vector>

// #include "checkUnits.h"

namespace Thermochimica
{
  void thermochimica();
  void setup();
  void solve();
  void init();
  void checkSystem();
  void compThermoData();

  void setThermoFilename(const std::string &filename);
  /** Exclude the named phases from subsequent equilibrium calculations; return zero on success. */
  int setExcludedPhases(const std::vector<std::string> &phases);
  /** Exclude every phase except the named phases; return zero on success. */
  int setIncludedPhases(const std::vector<std::string> &phases);
  void setUnitTemperature(const std::string &tunit);
  void setUnitPressure(const std::string &punit);
  void setUnitMass(const std::string &munit);
  void setUnits(const std::string &tunit, const std::string &punit, const std::string &munit);
  void setStandardUnits();
  void setModelicaUnits();
  void setTemperaturePressure(double temperature, double pressure);
  void presetElementMass(int element, double mass);
  void setElementMass(int element, double mass);
  int checkInfoThermo();
  void parseThermoFile();
  void setPrintResultsMode(int mode);
  void printResults();
  void resetInfoThermo();
  void resetThermo();
  void resetThermoAll();

  // utilitiy functions for consistency check / database record

  std::size_t getNumberElementsDatabase();
  std::vector<std::string> getElementsDatabase();
  std::string getElementAtIndex(int element_index);
  std::pair<std::size_t, std::size_t> getNumberPhasesSystem();
  std::vector<std::string> getPhaseNamesSystem();
  std::string getPhaseNameAtIndex(int phase_index);
  std::vector<std::size_t> getNumberSpeciesSystem();
  std::vector<std::string> getSpeciesInPhase(int phase_index);
  std::vector<std::vector<std::string>> getSpeciesSystem();
  std::vector<std::string> getThermodynamicSpeciesInPhase(int phase_index);

  // re-initialization-related functions
  void saveReinitData();
  std::pair<int, int> getReinitDataSizes();
  std::vector<double> getMolesPhase();
  std::vector<int> getAssemblage();
  void setReinitRequested(bool requested);
  void resetReinit();
  std::vector<double> getAllElementPotential();
  double getElementFraction(int atomicNumber);

  std::pair<double, int> getOutputChemPot(const std::string &elementName);
  std::tuple<double, double, int> getOutputSolnSpecies(const std::string &phaseName, const std::string &speciesName);
  std::tuple<double, double, int> getOutputMolSpecies(const std::string &speciesName);
  std::pair<double, int> getOutputMolSpeciesPhase(const std::string &phaseName, const std::string &speciesName);
  std::pair<double, int> getOutputMolSpeciesPhase(int phaseSystemIndex, int speciesPhaseIndex);
  std::pair<double, int> getElementMolesInPhase(const std::string &elementName, const std::string &phaseName);
  std::pair<double, int> getElementMolesInPhase(int elementSystemIndex, int phaseSystemIndex);
  std::pair<double, int> getElementMoleFractionInPhase(const std::string &elementName, const std::string &phaseName);
  std::pair<double, int> getSolnPhaseMol(const std::string &phaseName);
  std::pair<double, int> getPureConPhaseMol(const std::string &phaseName);
  std::pair<int, int> getPhaseIndex(const std::string &phaseName);
  std::pair<int, int> getPhaseIndex(int phaseSystemIndex);
  std::pair<double, int> getPhaseMoles(int phaseSystemIndex);
  std::pair<double, int> getSpeciesChemicalPotential(int phaseSystemIndex, int speciesPhaseIndex);
  std::pair<double, int> getMqmqaEndmemberStoichiometricPotential(int phaseSystemIndex, int endmemberIndex);
  struct PhaseGibbsEnergy
  {
    double total;
    double molar;
    int status;
  };
  PhaseGibbsEnergy getPhaseGibbsEnergy(int phaseSystemIndex);
  std::pair<double, int> getPhaseDrivingForce(int phaseSystemIndex);
  std::pair<double, int> getSystemGibbsEnergy();
  /** Return constituent names by sublattice for a zero-based system phase index. */
  std::vector<std::vector<std::string>> getConstituentsInPhase(int phaseSystemIndex);
  /** Return a constituent fraction using zero-based phase, sublattice, and constituent indices. */
  std::pair<double, int> getConstituentFraction(int phaseSystemIndex,
                                                int sublatticeIndex,
                                                int constituentIndex);
  std::pair<int, int> getElementIndex(int atomicNumber);
  std::pair<double, int> getElementPotential(int elementSystemIndex);
  std::pair<double, int> getOutputSiteFraction(const std::string &phaseName, int sublattice, int constituent);
  std::pair<double, int> getSublSiteMol(const std::string &phaseName, int sublattice, int constituent);

  bool isPhaseGas(const int phaseIndex);

  bool isPhaseMQM(const int phaseIndex);
  std::pair<double, int> getMqmqaMolesPairs(const std::string &phaseName);
  std::pair<double, int> getMqmqaPairMolFraction(const std::string &phaseName, const std::string &pairName);
  std::tuple<int, int, int> getMqmqaNumberPairsQuads(const std::string &phaseName);
  std::pair<double, int> getMqmqaConstituentFraction(const std::string &phaseName, int sublattice, const std::string &constituent);

  // Reinitialization data
  struct ReinitializationData
  {
    std::vector<int> assemblage;
    std::vector<double> molesPhase;
    std::vector<double> elementPotential;
    std::vector<double> chemicalPotential;
    std::vector<double> moleFraction;
    // Set elements used size here to match Thermochimica
    int nPeriodicTable = 169;
    std::vector<int> elementsUsed = std::vector<int>(nPeriodicTable);
    bool reinitAvailable;
    int GEM_iterations;
  };

  ReinitializationData getReinitData();
  void setReinitData(const ReinitializationData &data);

  struct ReinitializationDataView
  {
    int *assemblage;
    double *molesPhase;
    double *elementPotential;
    double *chemicalPotential;
    double *moleFraction;
    int *elementsUsed;
  };

  std::pair<bool, int> getReinitData(ReinitializationDataView data);
  void setReinitData(ReinitializationDataView data);

  // Heat capacity, enthalpy, and entropy
  void setHeatCapacityEnthalpyEntropyRequested(bool requested);
  std::tuple<double, double, double> getHeatCapacityEnthalpyEntropy();

  // Fuzzy stoichiometry
  void setFuzzyStoich(bool requested);
  void setFuzzyMagnitude(double magnitude);
  void setGibbsMinCheck(bool requested);
  void setMinMoleFraction(double min_mole_fraction);
  void setMassBalanceTolerance(double tolerance);

  /** Compute fixed-assemblage sensitivities for the current converged equilibrium state. */
  std::pair<int, double> computeEquilibriumSensitivity();
  /** Return the internal Kelvin, atmosphere, and elemental-mole sensitivity inputs. */
  std::pair<std::vector<double>, int> getEquilibriumSensitivityInputs();
  /** Apply a physical input displacement using the current fixed-assemblage linearization. */
  int applyEquilibriumSensitivity(const std::vector<double> &delta);
  /** Restore the exact state from which the current sensitivity was constructed. */
  int restoreEquilibriumSensitivity();
  /** Return the one-based indices of active constituents in stable solution phases. */
  std::pair<std::vector<int>, int> getActiveSolutionSpecies();
  /** Report whether every solution model in the current system is sensitivity-qualified. */
  std::pair<bool, int> getEquilibriumSensitivityQualification();
  struct EquilibriumSensitivityDiagnostics
  {
    double scaled_rcond = 0.0;
    double raw_rcond = 0.0;
    double base_residual = 0.0;
    double linear_residual = 0.0;
    double backward_error = 0.0;
    double minimum_active_fraction = 0.0;
    double minimum_active_moles = 0.0;
    bool qualified = false;
    int status = 0;
  };
  EquilibriumSensitivityDiagnostics getEquilibriumSensitivityDiagnostics();
  /** Return a versioned token describing the complete fixed active set. */
  std::pair<std::vector<int>, int> getEquilibriumAssemblageToken();
  struct EquilibriumSensitivityState
  {
    std::vector<double> log_amounts;
    /** Row-major derivatives of log amounts with respect to physical sensitivity inputs. */
    std::vector<double> jacobian;
    std::size_t parameters = 0;
    int status = 0;
  };
  EquilibriumSensitivityState getEquilibriumSensitivityState();

} // namespace Thermochimica
