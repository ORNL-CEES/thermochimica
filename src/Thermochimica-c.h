#include <string>

void SetThermoFilename(const char *);
void SetUnitTemperature(const char *);
void SetUnitPressure(const char *);
void SetUnitMass(const char *);
void SetUnits(const char *, const char *, const char *);

void GetOutputChemPot(char *, double *, int *);
void GetOutputSolnSpecies(const char *, char *, double *, double *, int *);
void GetOutputMolSpecies(const char *, double *, double *, int *);
void GetOutputMolSpeciesPhase(const char *, const char *, double *, int *);
void GetElementMolesInPhase(const char *, const char *, double *, int *);
void GetElementMoleFractionInPhase(const char *, const char *, double *, int *);

void GetSolnPhaseMol(const char *, double *, int *);
void GetPureConPhaseMol(const char *, double *, int *);
void GetPhaseIndex(const char *, int *, int *);

void GetOutputSiteFraction(const char *, int *, int *, double *, int *);
void GetSublSiteMol(const char *, int *, int *, double *, int *);

// MQMQA functions
void GetMqmqaMolesPairs(const char *, double *, int *);
void GetMqmqaPairMolFraction(const char *, const char *, double *, int *);
void GetMqmqaNumberPairsQuads(const char *, int *, int *, int *);

unsigned int atomicNumber(const char *);

unsigned int checkTemperature(const std::string &temperature_unit);
unsigned int checkPressure(const std::string &pressure_unit);
unsigned int checkMass(const std::string &mass_unit);

std::string elementName(unsigned int atomic_number);
