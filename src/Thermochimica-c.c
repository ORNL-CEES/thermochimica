#include <stdio.h>
#include <math.h>
#include <cstring>
#include <map>
#include "Thermochimica.h"

void ConvertToFortran(char* fstring, std::size_t fstring_len,
                      const char* cstring)
{
  std::size_t inlen = std::strlen(cstring);
  std::size_t cpylen = std::min(inlen, fstring_len);

  std::copy(cstring, cstring + cpylen, fstring);
  std::fill(fstring + cpylen, fstring + fstring_len, ' ');
}

void SetThermoFilename(const char* filename)
{
  char cfilename[120];
  int lcfilename = sizeof cfilename;
  ConvertToFortran(cfilename,lcfilename,filename);
  // printf("Datafile name set to %s\n",cfilename);

  FORTRAN_CALL(setthermofilename)(cfilename,&lcfilename);
}

void SetUnitTemperature(const char* unit)
{
  char cunit[15];
  int lcunit = sizeof cunit;
  ConvertToFortran(cunit,lcunit,unit);
  // printf("Temperature unit set to %s\n",unit);

  FORTRAN_CALL(setunittemperature)(cunit);
}

void SetUnitPressure(const char* unit)
{
  char cunit[15];
  int lcunit = sizeof cunit;
  ConvertToFortran(cunit,lcunit,unit);

  FORTRAN_CALL(setunitpressure)(cunit);
}

void SetUnitMass(const char* unit)
{
  char cunit[15];
  int lcunit = sizeof cunit;
  ConvertToFortran(cunit,lcunit,unit);

  FORTRAN_CALL(setunitmass)(cunit);
}

void SetUnits(const char* tunit, const char* punit, const char* munit)
{
  SetUnitTemperature(tunit);
  SetUnitPressure(punit);
  SetUnitMass(munit);
}

void GetOutputChemPot(char* elementName, double* chemPot, int* info)
{
  char celementName[25];
  int lcelementName = sizeof celementName;
  ConvertToFortran(celementName,lcelementName,elementName);

  FORTRAN_CALL(getoutputchempot)(celementName, chemPot, info);
}

void GetOutputSolnSpecies(const char* phaseName, const char* speciesName, double* moleFraction, double* chemPot, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  char cspeciesName[25];
  int lcspeciesName = sizeof cspeciesName;
  ConvertToFortran(cspeciesName,lcspeciesName,speciesName);

  FORTRAN_CALL(getoutputsolnspecies)(cphaseName, &lcphaseName, cspeciesName, &lcspeciesName, moleFraction, chemPot, info);
}

void GetOutputMolSpecies(const char* speciesName, double* moleFraction, double* moles, int* info)
{
  char cspeciesName[25];
  int lcspeciesName = sizeof cspeciesName;
  ConvertToFortran(cspeciesName,lcspeciesName,speciesName);

  FORTRAN_CALL(getoutputmolspecies)(cspeciesName, &lcspeciesName, moleFraction, moles, info);
}

void GetOutputMolSpeciesPhase(const char* phaseName, const char* speciesName, double* moleFraction, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  char cspeciesName[25];
  int lcspeciesName = sizeof cspeciesName;
  ConvertToFortran(cspeciesName,lcspeciesName,speciesName);

  FORTRAN_CALL(getoutputmolspeciesphase)(cphaseName, &lcphaseName, cspeciesName, &lcspeciesName, moleFraction, info);
}

void GetElementMolesInPhase(const char* elementName, const char* phaseName, double* molesElement, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  char celementName[3];
  int lcelementName = sizeof celementName;
  ConvertToFortran(celementName,lcelementName,elementName);

  FORTRAN_CALL(getelementmolesinphase)(celementName, &lcelementName, cphaseName, &lcphaseName, molesElement, info);
}

void GetElementMoleFractionInPhase(const char* elementName, const char* phaseName, double* molesElement, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  char celementName[3];
  int lcelementName = sizeof celementName;
  ConvertToFortran(celementName,lcelementName,elementName);

  FORTRAN_CALL(getelementmolefractioninphase)(celementName, &lcelementName, cphaseName, &lcphaseName, molesElement, info);
}

void GetSolnPhaseMol(const char* phaseName, double* molesPhase, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  FORTRAN_CALL(getsolnphasemol)(cphaseName, molesPhase, info);
}

void GetPureConPhaseMol(const char* phaseName, double* molesPhase, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  FORTRAN_CALL(getpureconphasemol)(cphaseName, molesPhase, info);
}

void GetPhaseIndex(const char* phaseName, int* index, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  FORTRAN_CALL(getphaseindex)(cphaseName, &lcphaseName, index, info);
}

void GetOutputSiteFraction(const char* phaseName, int* sublattice, int* constituent, double* siteFraction, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  FORTRAN_CALL(getoutputsitefraction)(cphaseName, &lcphaseName, sublattice, constituent, siteFraction, info);
}

void GetSublSiteMol(const char* phaseName, int* sublattice, int* constituent, double* siteMoles, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  FORTRAN_CALL(getsublsitemol)(cphaseName, &lcphaseName, sublattice, constituent, siteMoles, info);
}

// MQMQA functions

void GetMqmqaMolesPairs(const char* phaseName, double* molesPairs, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  FORTRAN_CALL(getmqmqamolespairs)(cphaseName, molesPairs, info);
}

void GetMqmqaPairMolFraction(const char* phaseName, const char* pairName, double* molesPairs, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  char cpairName[25];
  int lcpairName = sizeof cpairName;
  ConvertToFortran(cpairName,lcpairName,pairName);

  FORTRAN_CALL(getmqmqapairmolfraction)(cphaseName, &lcphaseName, cpairName, &lcpairName, molesPairs, info);
}

void GetMqmqaNumberPairsQuads(const char* phaseName, int* nPairs, int* nQuads, int* info)
{
  char cphaseName[25];
  int lcphaseName = sizeof cphaseName;
  ConvertToFortran(cphaseName,lcphaseName,phaseName);

  // GetMqmqaNumberPairsQuads was given a mismatched number of arguments between Fortran and C
  // The last argument in C is the string length.
  FORTRAN_CALL(getmqmqanumberpairsquads)(cphaseName, nPairs, nQuads, info, &lcphaseName);
}

struct cmp_str
{
   bool operator()(const char *a, const char *b) const
   {
      return std::strcmp(a, b) < 0;
   }
};

unsigned int
atomicNumber(const char* element)
{
  static std::map<const char*, unsigned int, cmp_str> _atomic_number_map = {
    {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5},
    {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
    {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15},
    {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20},
    {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25},
    {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35},
    {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40},
    {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45},
    {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50},
    {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55},
    {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60},
    {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65},
    {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70},
    {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75},
    {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
    {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85},
    {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
    {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95},
    {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100},
    {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105},
    {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110},
    {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114}, {"Mc", 115},
    {"Lv", 116}, {"Ts", 117}, {"Og", 118}
  };

  auto it = _atomic_number_map.find(element);
  if (it != _atomic_number_map.end())
    return it->second;
  else
    return -1;
}

