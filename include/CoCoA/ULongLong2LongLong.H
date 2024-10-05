#ifndef CoCoA_ULongLong2LongLong_H
#define CoCoA_ULongLong2LongLong_H

//   Copyright (c)  2022  John Abbott, and Anna M. Bigatti

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

  #ifndef CoCoA_ULONGLONG2LONGLONG
  // line below should trigger fatal compilation error
  "ERROR: -DCoCoA_ULONGLONG2LONGLONG=n setting missing";
  #endif
  
  // Convert back from unsigned long into signed long.
  // Definition used is determined by CPP flag COCOA_ULONGLONG2LONGLONG
  inline long long ULongLong2LongLong(unsigned long long n) noexcept
  {
#if CoCoA_ULONGLONG2LONGLONG == 1
    // Simplest defn; works on some platforms -- very fast.
    return n;
#elif CoCoA_ULONGLONG2LONGLONG == 2
    // Assumes same bit repn for long and ulong; probably fast but not clean.
    return reinterpret_cast<long long&>(n);
#else
    // A clean and portable implementation (but not so fast).
    // Assumes MaxULongLong >= MaxLongLong - MinLongLong & abs(MinLongLong) <= 1+MaxLongLong
    constexpr unsigned long MaxLongLong = std::numeric_limits<long long>::max();
    constexpr unsigned long MaxULongLong = std::numeric_limits<unsigned long long>::max();
    if (n <= MaxLongLong)
      return n;
    const unsigned long long negated = MaxULongLong - n;
    return -1-static_cast<long long>(negated);
#endif
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ULongLong2LongLong.H,v 1.3 2024/03/18 20:57:28 abbott Exp $
// $Log: ULongLong2LongLong.H,v $
// Revision 1.3  2024/03/18 20:57:28  abbott
// Summary: Removed const from return value type (redmine 1606)
//
// Revision 1.2  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.1  2022/02/16 21:30:56  abbott
// Summary: Added long long handling to MachineInt (redmine 581)
//
//
//

#endif
