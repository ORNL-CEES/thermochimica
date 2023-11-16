#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cstring>
#include <map>
#include <vector>
#include "Thermochimica.h"
#include "Thermochimica-cxx.h"

namespace Thermochimica
{

  void thermochimica()
  {
    TCAPI_thermochimica();
  }

  void setup()
  {
    TCAPI_setup();
  }

  void solve()
  {
    TCAPI_solve();
  }

  void init()
  {
    TCAPI_init();
  }

  void checkSystem()
  {
    TCAPI_checkSystem();
  }

  void compThermoData()
  {
    TCAPI_compThermoData();
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

  void presetElementMass(int element, double mass)
  {
    TCAPI_presetElementMass(&element, &mass);
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
    TCAPI_parseCSDataFile();
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

  // utilitiy functions for consistency check / database record

  std::size_t getNumberElementsDatabase()
  {
    int n_elements;
    TCAPI_getNumberElementsDatabase(&n_elements);
    return (std::size_t)n_elements;
  }

  std::vector<std::string> getElementsDatabase()
  {
    auto n_elements = getNumberElementsDatabase();

    std::vector<std::string> elements(n_elements);

    for (std::size_t i = 0; i < n_elements; ++i)
      elements[i] = getElementAtIndex(i);

    return elements;
  }

  std::string getElementAtIndex(int element_index)
  {
    int length;
    auto index = element_index + 1; // Fortran indexing adjustment

    char *buffer = TCAPI_getElementAtIndex(&index, &length);

    return std::string(buffer, buffer + length);
  }

  std::pair<std::size_t, std::size_t> getNumberPhasesSystem()
  {
    int n_solution, n_condensed;
    TCAPI_getNumberPhasesSystem(&n_solution, &n_condensed);
    return {(std::size_t)n_solution, (std::size_t)n_condensed};
  }

  std::vector<std::string> getPhaseNamesSystem()
  {
    auto [n_soln_phases, n_cond_phases] = getNumberPhasesSystem();

    auto n_phases = n_soln_phases + n_cond_phases;

    std::vector<std::string> phase_names(n_phases);

    for (std::size_t i = 0; i < n_phases; ++i)
      phase_names[i] = getPhaseNameAtIndex(i);

    return phase_names;
  }

  std::vector<std::size_t> getNumberSpeciesSystem()
  {
    auto [n_soln_phases, n_cond_phases] = getNumberPhasesSystem();
    (void)n_cond_phases;
    std::vector<int> n_sp(n_soln_phases);
    std::vector<std::size_t> n_species(n_soln_phases);
    TCAPI_getNumberSpeciesSystem(n_sp.data());

    for (std::size_t i = 0; i < n_soln_phases; ++i)
      n_species[i] = (std::size_t)n_sp[i];

    return n_species;
  }

  std::string getPhaseNameAtIndex(int phase_index)
  {
    int length;
    auto index = phase_index + 1; // Fortran indexing starts at 1 instead of 0

    char *buffer = TCAPI_getPhaseNameAtIndex(&index, &length);

    return std::string(buffer, buffer + length);
  }

  std::vector<std::vector<std::string>> getSpeciesSystem()
  {
    auto [n_soln_phases, n_cond_phases] = getNumberPhasesSystem();

    (void)n_cond_phases;
    std::vector<std::vector<std::string>> species(n_soln_phases);

    for (std::size_t i = 0; i < n_soln_phases; ++i)

      species[i] = getSpeciesInPhase(i);

    return species;
  }

  std::vector<std::string> getSpeciesInPhase(int phase_index)
  {
    int length, index;
    auto ph_index = phase_index + 1;
    auto n_species = getNumberSpeciesSystem();
    auto is_mqm = isPhaseMQM(phase_index);

    std::size_t n_species_phase;

    if (is_mqm)
    {
      auto [n_pairs, n_quads, idbg] = getMqmqaNumberPairsQuads(getPhaseNameAtIndex(phase_index));
      n_species_phase = n_pairs;
    }
    else
      n_species_phase = phase_index == 0 ? n_species[phase_index] : n_species[phase_index] - n_species[phase_index - 1];


    std::vector<std::string> species(n_species_phase);

    for (std::size_t i = 0; i < n_species_phase; ++i)
    {

      if (is_mqm)
      {
        index = i + 1;
        char *buffer = TCAPI_getMqmqaPairAtIndex(&ph_index, &index, &length);
        species[i] = std::string(buffer, buffer + length);
      }
      else
      {
        index = phase_index == 0 ? i + 1 : n_species[phase_index - 1] + i + 1;
        char *buffer = TCAPI_getSpeciesAtIndex(&index, &length);
        species[i] = std::string(buffer, buffer + length);
      }
    }

    return species;
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

  std::vector<double> getMolesPhase()
  {
    auto elements = getReinitDataSizes().first;
    std::vector<double> molesPhase(elements);
    TCAPI_getMolesPhase(molesPhase.data());
    return molesPhase;
  }

  std::vector<int> getAssemblage()
  {
    auto elements = getReinitDataSizes().first;
    std::vector<int> assemblage(elements);
    TCAPI_getAssemblage(assemblage.data());
    return assemblage;
  }

  void setReinitRequested(bool requested)
  {
    int req = (requested) ? 1 : 0;
    TCAPI_setReinitRequested(&req);
  }

  void resetReinit()
  {
    TCAPI_resetReinit();
  }

  std::vector<double> getAllElementPotential()
  {
    auto elements = getReinitDataSizes().first;
    std::vector<double> potential(elements);
    TCAPI_getAllElementPotential(potential.data());
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
    return std::make_tuple(moleFrac, chemPot, info);
  }

  std::tuple<double, double, int>
  getOutputMolSpecies(const std::string &speciesName)
  {
    double moleFrac, moles;
    int info;
    TCAPI_getOutputMolSpecies(speciesName.c_str(), speciesName.length(), &moleFrac, &moles, &info);
    return std::make_tuple(moleFrac, moles, info);
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

  // Gas phase functions
  bool isPhaseGas(const int phase_index)
  {
    bool isGas;
    auto index = phase_index + 1;
    TCAPI_isPhaseGas(&index, &isGas);

    return isGas;
  }

  // MQMQA functions
  bool isPhaseMQM(const int phase_index)
  {
    bool isMQM;
    auto index = phase_index + 1; // Adjust indexing for Fortran
    TCAPI_isPhaseMQM(&index, &isMQM);

    return isMQM;
  }

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
    return std::make_tuple(nPairs, nQuads, info);
  }

  std::pair<double, int>
  getMqmqaConstituentFraction(const std::string &phaseName, int sublattice, const std::string &constituent)
  {
    double moleFraction;
    int info;
    TCAPI_getMqmqaConstituentFraction(phaseName.c_str(), phaseName.length(), &sublattice, constituent.c_str(), constituent.length(), &moleFraction, &info);
    return {moleFraction, info};
  }

  ReinitializationData
  getReinitData()
  {
    ReinitializationData data;
    int available;
    auto [elements, species] = getReinitDataSizes();
    data.assemblage.resize(elements);
    data.molesPhase.resize(elements);
    data.elementPotential.resize(elements);
    data.chemicalPotential.resize(species);
    data.moleFraction.resize(species);

    TCAPI_getReinitData(data.assemblage.data(), data.molesPhase.data(), data.elementPotential.data(), data.chemicalPotential.data(), data.moleFraction.data(), data.elementsUsed.data(), &available);

    data.reinitAvailable = available;

    return data;
  }

  void
  setReinitData(const ReinitializationData & data)
  {
    auto [elements, species] = getReinitDataSizes();
    TCAPI_setReinitData(&elements, &species, data.assemblage.data(), data.molesPhase.data(), data.elementPotential.data(), data.chemicalPotential.data(), data.moleFraction.data(), data.elementsUsed.data());
  }

  void setHeatCapacityEnthalpyEntropyRequested(bool requested)
  {
    int req = (requested) ? 1 : 0;
    TCAPI_setHeatCapacityEnthalpyEntropyRequested(&req);
  }

  std::tuple<double, double, double> getHeatCapacityEnthalpyEntropy()
  {
    double heatCapacity, enthalpy, entropy;
    TCAPI_getHeatCapacityEnthalpyEntropy(&heatCapacity, &enthalpy, &entropy);
    return std::make_tuple(heatCapacity, enthalpy, entropy);
  }

  // Fuzzy stoichiometry
  void setFuzzyStoich(bool requested)
  {
    TCAPI_setFuzzyStoich(&requested);
  }

  void setFuzzyMagnitude(double magnitude)
  {
    TCAPI_setFuzzyMagnitude(&magnitude);
  }

  void setGibbsMinCheck(bool requested)
  {
    TCAPI_setGibbsMinCheck(&requested);
  }

}
