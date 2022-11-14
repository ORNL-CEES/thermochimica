#include <utility>
#include <tuple>
#include <string>
#include <vector>

// #include "checkUnits.h"

namespace Thermochimica
{

  void thermochimica();
  void setThermoFilename(const std::string &filename);
  void setUnitTemperature(const std::string &tunit);
  void setUnitPressure(const std::string &punit);
  void setUnitMass(const std::string &munit);
  void setUnits(const std::string &tunit, const std::string &punit, const std::string &munit);
  void setStandardUnits();
  void setModelicaUnits();
  void setTemperaturePressure(double temperature, double pressure);
  void setElementMass(int element, double mass);
  int checkInfoThermo();
  void parseThermoFile();
  void setPrintResultsMode(int mode);
  void printResults();
  void resetInfoThermo();
  void resetThermo();
  void resetThermoAll();

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

  std::pair<double, int>
  getMqmqaMolesPairs(const std::string &phaseName);
  std::pair<double, int>
  getMqmqaPairMolFraction(const std::string &phaseName, const std::string &pairName);
  std::tuple<int, int, int>
  getMqmqaNumberPairsQuads(const std::string &phaseName);
  std::pair<double, int>
  getMqmqaConstituentFraction(const std::string &phaseName, int sublattice, const std::string &constituent);

  // Reinitialization data
  struct reinitData {
    std::vector<int> assemblage;
    std::vector<double> molesPhase;
    std::vector<double> elementPotential;
    std::vector<double> chemicalPotential;
    std::vector<double> moleFraction;
    // Set elements used size here to match Thermochimica
    int nPeriodicTable = 169;
    std::vector<int> elementsUsed;
    bool reinitAvailable;
  };

  reinitData getReinitData();
  void setReinitData(reinitData data);

}
