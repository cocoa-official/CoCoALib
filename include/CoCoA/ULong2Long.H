#ifndef CoCoA_ULong2Long_H
#define CoCoA_ULong2Long_H

//   Copyright (c)  2011  John Abbott, and Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "CoCoA/PREPROCESSOR_DEFNS.H"

#include <limits>
// using std::numeric_limits

namespace CoCoA
{

  #ifndef CoCoA_ULONG2LONG
  // line below should trigger fatal compilation error
  "ERROR: -DCoCoA_ULONG2LONG=n setting missing";
  #endif
  
  // Convert back from unsigned long into signed long.
  // Definition used is determined by CPP flag COCOA_ULONG2LONG
  inline long ULong2Long(unsigned long n) noexcept
  {
#if CoCoA_ULONG2LONG == 1
    // Simplest defn; works on some platforms -- very fast.
    return n;
#elif CoCoA_ULONG2LONG == 2
    // Assumes same bit repn for long and ulong; probably fast but not clean.
    return reinterpret_cast<long&>(n);
#else
    // A clean and portable implementation (but not so fast).
    // Assumes MaxULong >= MaxLong - MinLong & abs(MinLong) <= 1+MaxLong
    constexpr unsigned long MaxLong = std::numeric_limits<long>::max();
    constexpr unsigned long MaxULong = std::numeric_limits<unsigned long>::max();
    if (n <= MaxLong)
      return n;
    const unsigned long negated = MaxULong - n;
    return -1-static_cast<long>(negated);
#endif
  }


} // end of namespace CoCoA

#endif
