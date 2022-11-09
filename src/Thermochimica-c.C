#include <stdio.h>
#include <math.h>
#include <cstring>
#include <map>
#include "Thermochimica.h"

void SetThermoFilename(const char *filename)
{
  TCAPI_setThermoFilename(filename, strlen(filename));
}

void SetUnitTemperature(const char *unit)
{
  TCAPI_setUnitTemperature(unit, strlen(unit));
}

void SetUnitPressure(const char *unit)
{
  TCAPI_setUnitPressure(unit, strlen(unit));
}

void SetUnitMass(const char *unit)
{
  TCAPI_setUnitMass(unit, strlen(unit));
}

void SetUnits(const char *tunit, const char *punit, const char *munit)
{
  SetUnitTemperature(tunit);
  SetUnitPressure(punit);
  SetUnitMass(munit);
}

void GetOutputChemPot(char *elementName, double *chemPot, int *info)
{
  TCAPI_getOutputChemPot(elementName, strlen(elementName), chemPot, info);
}

void GetOutputSolnSpecies(const char *phaseName, const char *speciesName, double *moleFraction, double *chemPot, int *info)
{
  TCAPI_getOutputSolnSpecies(phaseName, strlen(phaseName), speciesName, strlen(speciesName), moleFraction, chemPot, info);
}

void GetOutputMolSpecies(const char *speciesName, double *moleFraction, double *moles, int *info)
{
  TCAPI_getOutputMolSpecies(speciesName, strlen(speciesName), moleFraction, moles, info);
}

void GetOutputMolSpeciesPhase(const char *phaseName, const char *speciesName, double *moleFraction, int *info)
{
  TCAPI_getOutputMolSpeciesPhase(phaseName, strlen(phaseName), speciesName, strlen(speciesName), moleFraction, info);
}

void GetElementMolesInPhase(const char *elementName, const char *phaseName, double *molesElement, int *info)
{
  TCAPI_getElementMolesInPhase(elementName, strlen(elementName), phaseName, strlen(phaseName), molesElement, info);
}

void GetElementMoleFractionInPhase(const char *elementName, const char *phaseName, double *molesElement, int *info)
{
  TCAPI_getElementMoleFractionInPhase(elementName, strlen(elementName), phaseName, strlen(phaseName), molesElement, info);
}

void GetSolnPhaseMol(const char *phaseName, double *molesPhase, int *info)
{
  TCAPI_getSolnPhaseMol(phaseName, strlen(phaseName), molesPhase, info);
}

void GetPureConPhaseMol(const char *phaseName, double *molesPhase, int *info)
{
  TCAPI_getPureConPhaseMol(phaseName, strlen(phaseName), molesPhase, info);
}

void GetPhaseIndex(const char *phaseName, int *index, int *info)
{
  TCAPI_getPhaseIndex(phaseName, strlen(phaseName), index, info);
}

void GetOutputSiteFraction(const char *phaseName, int *sublattice, int *constituent, double *siteFraction, int *info)
{
  TCAPI_getOutputSiteFraction(phaseName, strlen(phaseName), sublattice, constituent, siteFraction, info);
}

void GetSublSiteMol(const char *phaseName, int *sublattice, int *constituent, double *siteMoles, int *info)
{
  TCAPI_getSublSiteMol(phaseName, strlen(phaseName), sublattice, constituent, siteMoles, info);
}

// MQMQA functions

void GetMqmqaMolesPairs(const char *phaseName, double *molesPairs, int *info)
{
  TCAPI_getMqmqaMolesPairs(phaseName, strlen(phaseName), molesPairs, info);
}

void GetMqmqaPairMolFraction(const char *phaseName, const char *pairName, double *molesPairs, int *info)
{
  TCAPI_getMqmqaPairMolFraction(phaseName, strlen(phaseName), pairName, strlen(pairName), molesPairs, info);
}

void GetMqmqaNumberPairsQuads(const char *phaseName, int *nPairs, int *nQuads, int *info)
{
  TCAPI_getMqmqaNumberPairsQuads(phaseName, strlen(phaseName), nPairs, nQuads, info);
}

struct cmp_str
{
  bool operator()(const char *a, const char *b) const
  {
    return std::strcmp(a, b) < 0;
  }
};

unsigned int
atomicNumber(const char *element)
{
  static std::map<const char *, unsigned int, cmp_str> _atomic_number_map = {
      {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40}, {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80}, {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90}, {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100}, {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110}, {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114}, {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}};

  auto it = _atomic_number_map.find(element);
  if (it != _atomic_number_map.end())
    return it->second;
  else
    return -1;
}
