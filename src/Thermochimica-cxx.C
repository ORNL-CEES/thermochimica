#include <stdio.h>
#include <math.h>
#include <cstring>
#include <map>
#include <vector>
#include "Thermochimica.h"

namespace Thermochimica
{

  void thermochimica()
  {
    TCAPI_thermochimica();
  }

  void setThermoFilename(const std::string &filename)
  {
    TCAPI_setThermoFilename(filename.c_str(), filename.length());
  }

  void setUnitTemperature(const std::string &tunit)
  {
    TCAPI_setUnitTemperature(tunit.c_str(), tunit.length());
  }

  void setUnitPressure(const std::string &punit)
  {
    TCAPI_setUnitPressure(punit.c_str(), punit.length());
  }

  void setUnitMass(const std::string &munit)
  {
    TCAPI_setUnitMass(munit.c_str(), munit.length());
  }

  void setUnits(const std::string &tunit, const std::string &punit, const std::string &munit)
  {
    setUnitTemperature(tunit);
    setUnitPressure(punit);
    setUnitMass(munit);
  }

  void setStandardUnits()
  {
    TCAPI_setStandardUnits();
  }

  void setModelicaUnits()
  {
    TCAPI_setModelicaUnits();
  }

  void setTemperaturePressure(double temperature, double pressure)
  {
    TCAPI_setTemperaturePressure(&temperature, &pressure);
  }

  void setElementMass(int element, double mass)
  {
    TCAPI_setElementMass(&element, &mass);
  }

  int checkInfoThermo()
  {
    int idbg;
    TCAPI_checkInfoThermo(&idbg);
    return idbg;
  }

  void parseThermoFile()
  {
    TCAPI_sSParseCSDataFile();
  }

  void setPrintResultsMode(int mode)
  {
    TCAPI_setPrintResultsMode(&mode);
  }

  void printResults()
  {
    TCAPI_printResults();
  }

  void resetInfoThermo()
  {
    TCAPI_resetInfoThermo();
  }

  void resetThermo()
  {
    TCAPI_resetThermo();
  }

  void resetThermoAll()
  {
    TCAPI_resetThermoAll();
  }

  // re-initialization-related functions
  void saveReinitData()
  {
    TCAPI_saveReinitData();
  }

  std::pair<int, int> getReinitDataSizes()
  {
    int elements, species;
    TCAPI_getReinitDataSizes(&elements, &species);
    return {elements, species};
  }

  double getMolesPhase()
  {
    double molesPhase;
    TCAPI_getMolesPhase(&molesPhase);
    return molesPhase;
  }

  std::vector<int> getAssemblage()
  {
    std::pair<int, int> p = getReinitDataSizes();
    int assemblageArray[p.first];

    TCAPI_getAssemblage(assemblageArray);

    std::vector<int> assemblage(assemblageArray, assemblageArray + sizeof assemblageArray / sizeof assemblageArray[0]);

    return assemblage;
  }

  void setReinitRequested(bool requested)
  {
    int req = (requested) ? req = 1 : 0;
    TCAPI_setReinitRequested(&req);
  }

  void resetReinit()
  {
    TCAPI_resetReinit();
  }

  std::vector<double> getAllElementPotential()
  {
    std::pair<int, int> p = getReinitDataSizes();
    double potentialArray[p.first];

    TCAPI_getAllElementPotential(potentialArray);

    std::vector<double> potential(potentialArray, potentialArray + sizeof potentialArray / sizeof potentialArray[0]);

    return potential;
  }

  double getElementFraction(int atomicNumber)
  {
    double fraction;
    TCAPI_getElementFraction(&atomicNumber, &fraction);

    return fraction;
  }

  // Data extraction APIs

  std::pair<double, int>
  getOutputChemPot(const std::string &elementName)
  {
    double chemPot;
    int info;
    TCAPI_getOutputChemPot(elementName.c_str(), elementName.length(), &chemPot, &info);
    return {chemPot, info};
  }

  std::tuple<double, double, int>
  getOutputSolnSpecies(const std::string &phaseName, const std::string &speciesName)
  {
    double moleFrac, chemPot;
    int info;
    TCAPI_getOutputSolnSpecies(phaseName.c_str(), phaseName.length(), speciesName.c_str(), speciesName.length(), &moleFrac, &chemPot, &info);
    return {moleFrac, chemPot, info};
  }

  std::tuple<double, double, int>
  getOutputMolSpecies(const std::string &speciesName)
  {
    double moleFrac, moles;
    int info;
    TCAPI_getOutputMolSpecies(speciesName.c_str(), speciesName.length(), &moleFrac, &moles, &info);
    return {moleFrac, moles, info};
  }

  std::pair<double, int>
  getOutputMolSpeciesPhase(const std::string &phaseName, const std::string &speciesName)
  {
    double moleFrac;
    int info;
    TCAPI_getOutputMolSpeciesPhase(phaseName.c_str(), phaseName.length(), speciesName.c_str(), speciesName.length(), &moleFrac, &info);
    return {moleFrac, info};
  }

  std::pair<double, int>
  getElementMolesInPhase(const std::string &elementName, const std::string &phaseName)
  {
    double molesElement;
    int info;
    TCAPI_getElementMolesInPhase(elementName.c_str(), elementName.length(), phaseName.c_str(), phaseName.length(), &molesElement, &info);
    return {molesElement, info};
  }

  std::pair<double, int>
  getElementMoleFractionInPhase(const std::string &elementName, const std::string &phaseName)
  {
    double moleFracElement;
    int info;
    TCAPI_getElementMoleFractionInPhase(elementName.c_str(), elementName.length(), phaseName.c_str(), phaseName.length(), &moleFracElement, &info);
    return {moleFracElement, info};
  }

  std::pair<double, int>
  getSolnPhaseMol(const std::string &phaseName)
  {
    double molesPhase;
    int info;
    TCAPI_getSolnPhaseMol(phaseName.c_str(), phaseName.length(), &molesPhase, &info);
    return {molesPhase, info};
  }

  std::pair<double, int>
  getPureConPhaseMol(const std::string &phaseName)
  {
    double molesPhase;
    int info;
    TCAPI_getPureConPhaseMol(phaseName.c_str(), phaseName.length(), &molesPhase, &info);
    return {molesPhase, info};
  }

  std::pair<int, int>
  getPhaseIndex(const std::string &phaseName)
  {
    int index, info;
    TCAPI_getPhaseIndex(phaseName.c_str(), phaseName.length(), &index, &info);
    return {index, info};
  }

  std::pair<double, int>
  getOutputSiteFraction(const std::string &phaseName, int sublattice, int constituent)
  {
    double siteFraction;
    int info;
    TCAPI_getOutputSiteFraction(phaseName.c_str(), phaseName.length(), &sublattice, &constituent, &siteFraction, &info);
    return {siteFraction, info};
  }

  std::pair<double, int>
  getSublSiteMol(const std::string &phaseName, int sublattice, int constituent)
  {
    double siteMoles;
    int info;
    TCAPI_getSublSiteMol(phaseName.c_str(), phaseName.length(), &sublattice, &constituent, &siteMoles, &info);
    return {siteMoles, info};
  }

  // MQMQA functions

  std::pair<double, int>
  getMqmqaMolesPairs(const std::string &phaseName)
  {
    double molesPairs;
    int info;
    TCAPI_getMqmqaMolesPairs(phaseName.c_str(), phaseName.length(), &molesPairs, &info);
    return {molesPairs, info};
  }

  std::pair<double, int>
  getMqmqaPairMolFraction(const std::string &phaseName, const std::string &pairName)
  {
    double molesPairs;
    int info;
    TCAPI_getMqmqaPairMolFraction(phaseName.c_str(), phaseName.length(), pairName.c_str(), pairName.length(), &molesPairs, &info);
    return {molesPairs, info};
  }

  std::tuple<int, int, int>
  getMqmqaNumberPairsQuads(const std::string &phaseName)
  {
    int nPairs;
    int nQuads;
    int info;
    TCAPI_getMqmqaNumberPairsQuads(phaseName.c_str(), phaseName.length(), &nPairs, &nQuads, &info);
    return {nPairs, nQuads, info};
  }

}
