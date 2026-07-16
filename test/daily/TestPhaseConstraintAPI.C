#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

#include "Thermochimica.h"
#include "Thermochimica-c.h"
#include "Thermochimica-cxx.h"

#ifndef THERMOCHIMICA_DATA_DIRECTORY
#define THERMOCHIMICA_DATA_DIRECTORY "data/"
#endif

namespace
{
bool closeEnough(double lhs, double rhs, double tolerance = 1.0e-8)
{
  return std::abs(lhs - rhs) <= tolerance;
}
}

int main()
{
  bool pass = true;
  int info = 0;

  Thermochimica::resetThermoAll();
  Thermochimica::setThermoFilename(std::string(THERMOCHIMICA_DATA_DIRECTORY) + "CO.dat");
  Thermochimica::setStandardUnits();
  Thermochimica::setTemperaturePressure(1000.0, 1.0);
  Thermochimica::setElementMass(0, 0.0);
  Thermochimica::setElementMass(6, 1.0);
  Thermochimica::setElementMass(8, 1.0);
  Thermochimica::parseThermoFile();

  pass = pass && (Thermochimica::addPhaseFractionConstraint("gas_ideal", 0.6) == 0);
  pass = pass && (Thermochimica::addPhaseFractionConstraint("C_Graphite(s)", 0.4) == 0);
  pass = pass && (Thermochimica::getNumberPhaseFractionConstraints() == 2);

  auto [beforeSolve, beforeSolveInfo] = Thermochimica::getPhaseFractionConstraintAtIndex(0);
  pass = pass && (beforeSolveInfo == 2);
  pass = pass && beforeSolve.phaseName.empty();
  pass = pass && closeEnough(beforeSolve.targetFraction, 0.0);

  Thermochimica::thermochimica();
  pass = pass && (Thermochimica::checkInfoThermo() == 0);

  auto [gas, gasInfo] = Thermochimica::getPhaseFractionConstraintAtIndex(0);
  pass = pass && (gasInfo == 0);
  pass = pass && (gas.phaseName == "gas_ideal");
  pass = pass && closeEnough(gas.targetFraction, 0.6);
  pass = pass && closeEnough(gas.achievedFraction, 0.6);
  pass = pass && closeEnough(gas.residual, 0.0);

  int fortranIndex = 1;
  int nameLength = 0;
  double target = 0.0;
  double achieved = 0.0;
  double residual = 0.0;
  double multiplier = 0.0;
  char *name = TCAPI_getPhaseFractionConstraintAtIndex(&fortranIndex,
                                                        &nameLength,
                                                        &target,
                                                        &achieved,
                                                        &residual,
                                                        &multiplier,
                                                        &info);
  pass = pass && (info == 0);
  pass = pass && (name != nullptr);
  pass = pass && (std::string(name, name + nameLength) == gas.phaseName);
  pass = pass && closeEnough(target, gas.targetFraction);
  pass = pass && closeEnough(achieved, gas.achievedFraction);
  pass = pass && closeEnough(residual, gas.residual);
  pass = pass && closeEnough(multiplier, gas.lagrangeMultiplier);

  int constraintCount = 0;
  GetNumberPhaseFractionConstraints(&constraintCount);
  pass = pass && (constraintCount == 2);

  fortranIndex = 2;
  name = GetPhaseFractionConstraintAtIndex(&fortranIndex,
                                           &nameLength,
                                           &target,
                                           &achieved,
                                           &residual,
                                           &multiplier,
                                           &info);
  pass = pass && (info == 0);
  pass = pass && (name != nullptr);
  pass = pass && (std::string(name, name + nameLength) == "C_Graphite(s)");
  pass = pass && closeEnough(target, 0.4);
  pass = pass && closeEnough(achieved, 0.4);
  pass = pass && closeEnough(residual, 0.0);

  pass = pass && (Thermochimica::addPhaseFractionConstraint("gas_ideal", 0.6) == 0);
  auto [afterUpdate, afterUpdateInfo] = Thermochimica::getPhaseFractionConstraintAtIndex(0);
  pass = pass && (afterUpdateInfo == 2);
  pass = pass && afterUpdate.phaseName.empty();

  fortranIndex = 3;
  name = TCAPI_getPhaseFractionConstraintAtIndex(&fortranIndex,
                                                  &nameLength,
                                                  &target,
                                                  &achieved,
                                                  &residual,
                                                  &multiplier,
                                                  &info);
  pass = pass && (info == 1);
  pass = pass && (name == nullptr);
  pass = pass && (nameLength == 0);
  pass = pass && closeEnough(target, 0.0);
  pass = pass && closeEnough(achieved, 0.0);
  pass = pass && closeEnough(residual, 0.0);
  pass = pass && closeEnough(multiplier, 0.0);

  Thermochimica::resetThermoAll();
  pass = pass && (Thermochimica::getNumberPhaseFractionConstraints() == 0);

  if (pass)
  {
    std::cout << " TestPhaseConstraintAPI: PASS\n";
    return 0;
  }

  std::cout << " TestPhaseConstraintAPI: FAIL <---\n";
  return 1;
}
