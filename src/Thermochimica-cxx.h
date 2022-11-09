#include <utility>
#include <tuple>
#include <string>

#include "checkUnits.h"

namespace Thermocimica
{

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
  GetMqmqaPairMolFraction(const char *phaseName, const char *pairName, double *molesPairs, int *info);
  std::tuple<int, int, int>
  GetMqmqaNumberPairsQuads(const std::string &phaseName);

}