#include <iostream>
#include <string>
#include <algorithm>
#include <map>

namespace Thermochimica
{

  void
  checkTemperature(const std::string &temperature_unit)
  {
    static std::vector<std::string> _temperature_names = {"K", "C", "F", "R"};

    if (std::find(_temperature_names.begin(), _temperature_names.end(), temperature_unit) == _temperature_names.end())
      throw std::range_error("checkTemperature: unknown temperature unit '" + temperature_unit + "'.");
  }

  void
  checkPressure(const std::string &pressure_unit)
  {
    static std::vector<std::string> _pressure_names = {"atm", "psi", "bar", "Pa", "kPa"};

    if (std::find(_pressure_names.begin(), _pressure_names.end(), pressure_unit) == _pressure_names.end())
      throw std::range_error("checkPressure: unknown pressure unit '" + pressure_unit + "'.");
  }

  void
  checkMass(const std::string &mass_unit)
  {
    static std::vector<std::string> _mass_names = {
        "mole fraction", "atom fraction", "atoms", "moles",
        "gram-atoms", "mass fraction", "kilograms", "grams",
        "pounds"};

    if (std::find(_mass_names.begin(), _mass_names.end(), mass_unit) == _mass_names.end())
      throw std::range_error("checkMass: unknown mass unit '" + mass_unit + "'.");
  }

  unsigned int
  atomicNumber(const std::string &element)
  {
    static std::map<std::string, unsigned int> _atomic_number_map = {
        {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40}, {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80}, {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90}, {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100}, {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110}, {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114}, {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}};

    auto it = _atomic_number_map.find(element);

    if (it != _atomic_number_map.end())
      return it->second;
    else
      throw std::range_error("atomicNumber: unknown chemical element '" + element + "'.");
  }

  std::string
  elementName(unsigned int atomic_number)
  {
    static std::vector<std::string> _element_names = {
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

    if (atomic_number >= 1 && atomic_number <= _element_names.size())
      return _element_names[atomic_number - 1];
    else
      throw std::range_error("elementName: Requesting an element with atomic number " + std::to_string(atomic_number) + ", which is outside of the known range.");
  }
}
