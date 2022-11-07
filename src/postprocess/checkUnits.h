#include <iostream>
#include <string>
#include <algorithm>
#include "MooseApp.h"

namespace Thermochimica
{

  unsigned int
  checkTemperature(const std::string &temperature_unit);

  unsigned int
  checkPressure(const std::string &pressure_unit);

  unsigned int
  checkMass(const std::string &mass_unit);

  unsigned int
  atomicNumber(const std::string &element);

  std::string
  elementName(unsigned int atomic_number);

}
