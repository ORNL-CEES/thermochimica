#include <string>
#define FORTRAN_CALL(name) name##_

extern "C"
{
  void FORTRAN_CALL(setthermofilename)(const char *, int *);
  void SetThermoFilename(const char *);
  void FORTRAN_CALL(setunittemperature)(const char *);
  void SetUnitTemperature(const char *);
  void FORTRAN_CALL(setunitpressure)(const char *);
  void SetUnitPressure(const char *);
  void FORTRAN_CALL(setunitmass)(const char *);
  void SetUnitMass(const char *);
  void FORTRAN_CALL(setunits)(const char *, const char *, const char *);
  void SetUnits(const char *, const char *, const char *);
  void FORTRAN_CALL(setstandardunits)();

  void FORTRAN_CALL(getmolfraction)(int *, double *, int *);
  void FORTRAN_CALL(getchemicalpotential)(int *, double *, int *);
  void FORTRAN_CALL(getelementpotential)(int *, double *, int *);

  void FORTRAN_CALL(getoutputchempot)(char *, double *, int *);
  void GetOutputChemPot(char *, double *, int *);
  void FORTRAN_CALL(getoutputsolnspecies)(const char *, int *, const char *, int *, double *, double *, int *);
  void GetOutputSolnSpecies(const char *, char *, double *, double *, int *);
  void FORTRAN_CALL(getoutputmolspecies)(const char *, int *, double *, double *, int *);
  void GetOutputMolSpecies(const char *, double *, double *, int *);
  void FORTRAN_CALL(getoutputmolspeciesphase)(const char *, int *, const char *, int *, double *, int *);
  void GetOutputMolSpeciesPhase(const char *, const char *, double *, int *);
  void FORTRAN_CALL(getelementmolesinphase)(const char *, int *, const char *, int *, double *, int *);
  void GetElementMolesInPhase(const char *, const char *, double *, int *);
  void FORTRAN_CALL(getelementmolefractioninphase)(const char *, int *, const char *, int *, double *, int *);
  void GetElementMoleFractionInPhase(const char *, const char *, double *, int *);

  void FORTRAN_CALL(getsolnphasemol)(const char *, double *, int *);
  void GetSolnPhaseMol(const char *, double *, int *);
  void FORTRAN_CALL(getpureconphasemol)(const char *, double *, int *);
  void GetPureConPhaseMol(const char *, double *, int *);
  void FORTRAN_CALL(getphaseindex)(const char *, int *, int *, int *);
  void GetPhaseIndex(const char *, int *, int *);

  void FORTRAN_CALL(getoutputsitefraction)(const char *, int *, int *, int *, double *, int *);
  void GetOutputSiteFraction(const char *, int *, int *, double *, int *);
  void FORTRAN_CALL(getsublsitemol)(const char *, int *, int *, int *, double *, int *);
  void GetSublSiteMol(const char *, int *, int *, double *, int *);

  void FORTRAN_CALL(setprintresultsmode)(int *);

  void FORTRAN_CALL(setelementmass)(int *, double *);
  void FORTRAN_CALL(presetelementmass)(int *, double *);
  void FORTRAN_CALL(settemperaturepressure)(double *, double *);
  void FORTRAN_CALL(checkinfothermo)(int *);


  void FORTRAN_CALL(ssparsecsdatafile)();
  void FORTRAN_CALL(thermochimica)();

  void FORTRAN_CALL(solphaseparse)(int *, double *);

  void FORTRAN_CALL(thermodebug)();

  void FORTRAN_CALL(printresults)();
  void FORTRAN_CALL(printstate)();

  void FORTRAN_CALL(resetinfothermo)();
  void FORTRAN_CALL(resetthermo)();
  void FORTRAN_CALL(resetthermoall)();

  // re-initialization-related functions
  void FORTRAN_CALL(savereinitdata)();
  void FORTRAN_CALL(getreinitdatasizes)(int *, int *);
  void FORTRAN_CALL(getmolesphase)(double *);
  void FORTRAN_CALL(getassemblage)(int *);
  void FORTRAN_CALL(setreinitrequested)(int *);
  void FORTRAN_CALL(resetreinit)();
  void FORTRAN_CALL(getallelementpotential)(double *);
  void FORTRAN_CALL(getelementfraction)(int *, double *);

  unsigned int atomicNumber(const char *);
}

void ConvertToFortran(char *, std::size_t, const char *);

unsigned int checkTemperature(const std::string & temperature_unit);
unsigned int checkPressure(const std::string & pressure_unit);
unsigned int checkMass(const std::string & mass_unit);

std::string elementName(unsigned int atomic_number);
