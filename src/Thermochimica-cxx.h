#include <utility>
#include <tuple>
#include <string>
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

  // re-initialization-related functions
  void saveReinitData();
  std::pair<int, int> getReinitDataSizes();
  std::vector<double> getMolesPhase();
  std::vector<int> getAssemblage();
  void setReinitRequested(bool requested);
  void resetReinit();
  std::vector<double> getAllElementPotential();
  double getElementFraction(int atomicNumber);

  std::pair<double, int>
  getOutputChemPot(const std::string &elementName);
  std::tuple<double, double, int>
  getOutputSolnSpecies(const std::string &phaseName, const std::string &speciesName);
  std::tuple<double, double, int>
  getOutputMolSpecies(const std::string &speciesName);
  std::pair<double, int>
  getOutputMolSpeciesPhase(const std::string &phaseName, const std::string &speciesName);
  std::pair<double, int>
  getElementMolesInPhase(const std::string &elementName, const std::string &phaseName);
  std::pair<double, int>
  getElementMoleFractionInPhase(const std::string &elementName, const std::string &phaseName);
  std::pair<double, int>
  getSolnPhaseMol(const std::string &phaseName);
  std::pair<double, int>
  getPureConPhaseMol(const std::string &phaseName);
  std::pair<int, int>
  getPhaseIndex(const std::string &phaseName);
  std::pair<double, int>
  getOutputSiteFraction(const std::string &phaseName, int sublattice, int constituent);
  std::pair<double, int>
  getSublSiteMol(const std::string &phaseName, int sublattice, int constituent);

  bool isPhaseGas(const int phaseIndex);

  bool isPhaseMQM(const int phaseIndex);
  std::pair<double, int>
  getMqmqaMolesPairs(const std::string &phaseName);
  std::pair<double, int>
  getMqmqaPairMolFraction(const std::string &phaseName, const std::string &pairName);
  std::tuple<int, int, int>
  getMqmqaNumberPairsQuads(const std::string &phaseName);
  std::pair<double, int>
  getMqmqaConstituentFraction(const std::string &phaseName, int sublattice, const std::string &constituent);

  // Reinitialization data
  struct reinitData
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
  };

  reinitData getReinitData();
  void setReinitData(const reinitData &data);

  // Heat capacity, enthalpy, and entropy
  void setHeatCapacityEnthalpyEntropyRequested(bool requested);
  std::tuple<double, double, double> getHeatCapacityEnthalpyEntropy();

  // Fuzzy stoichiometry
  void setFuzzyStoich(bool requested);
  void setFuzzyMagnitude(double magnitude);
  void setGibbsMinCheck(bool requested);

}
