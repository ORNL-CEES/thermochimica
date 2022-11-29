#include <iostream>
#include <string>
#include <algorithm>
#include "MooseApp.h"

namespace Thermochimica
{

void ConvertToFortran(char* fstring, std::size_t fstring_len,
                      const char* cstring)
{
  std::size_t inlen = std::strlen(cstring);
  std::size_t cpylen = std::min(inlen, fstring_len);

  if (inlen > fstring_len)
    {
      Moose::out << "ConvertToFortran " << inlen << " " << fstring_len << "\n";
    }

  std::copy(cstring, cstring + cpylen, fstring);
  std::fill(fstring + cpylen, fstring + fstring_len, ' ');
}
}
