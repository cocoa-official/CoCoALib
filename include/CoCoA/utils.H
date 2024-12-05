#ifndef CoCoA_utils_H
#define CoCoA_utils_H

//   Copyright (c)  2006-2009,2011  John Abbott, and Anna M. Bigatti

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

#include "CoCoA/convert.H"

#include <istream>
// using std::istream
#include <ostream>
// using std::ostream
#include <limits>
// using std::numeric_limits
#include <cmath>
// using std::floor;
#include <string>
// using std::string;
#include <type_traits>
// using std::is_signed, std::make_unsigned
#include <vector>
// using std::vector;


namespace CoCoA
{

  // I think this fn should really be in the C++ standard library.
  // Round to nearest integer (rounding halves towards +infinity).
  inline double round(double x) noexcept // NOT constexpr (std::floor)
  {
    return std::floor(x+0.5);
  }


  inline int cmp(long n1, long n2) noexcept
  {
    if (n1 == n2) return 0;
    return (n1 < n2)? -1 : 1;
  }

  inline int cmp_ul(unsigned long ul1, unsigned long ul2) noexcept
  {
    if (ul1 == ul2) return 0;
    return (ul1 < ul2)? -1 : 1;
  }

// BUG BUG: changing this fn name to "cmp" causes several overload resolution problems :-(
  // This weird template stuff came from "cppreference" (dumbed down to be C++14 compatible)
  // SEE NEATER IMPL IMMEDIATELY BELOW (but needs C++20)
  template< class T, class U >
  constexpr int CMP( T t, U u ) noexcept
  {
    using UT = std::make_unsigned_t<T>;
    using UU = std::make_unsigned_t<U>;
    if /*constexpr*/ (std::is_signed<T>::value == std::is_signed<U>::value)
                   return (t == u)?0:(t < u)?-1:1;
    else if /*constexpr*/ (std::is_signed<T>::value)
                        return (t < 0) ? -1 : (UT(t) == u)? 0 : (UT(t) < u) ? -1 : 1;
    else
      return (u < 0) ? -1 : (t == UU(u))? 0 : (t < UU(u)) ? -1 : 1;
  }
 

  // // This might/should work with C++20
  // template<typename T1, typename T2> // only for builtin integral types
  // inline cmp(T1 a, T2 b)
  // {
  //   if (std::cmp_equal(a,b)) return 0;
  //   return (std:cmp_less(a,b))? -1 : 1;
  // }


  template< class T, class U >
  constexpr bool IsInRange_( T lo, U n, T hi ) noexcept
  {
    return (CMP(lo,n) <= 0) && (CMP(n,hi) <= 0);
  }

  // Computes hi - lo; assumes hi >= lo.
  inline unsigned long ULongDiff(long hi, long lo) noexcept
  {
    CoCoA_ASSERT(lo <= hi);
    return static_cast<unsigned long>(hi) - lo; // lo is silently converted to ulong
  }


  inline std::vector<long> LongRange(long lo, long hi) // lo & hi both included
  {
    if (lo > hi)  return std::vector<long>();
    const unsigned long n = ULongDiff(hi, lo);
    if (n+1 == 0)  CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    std::vector<long> ans(n+1);
    ans[n] = hi; // this way it works even when hi=LONG_MAX
    std::vector<long>::iterator ptr = ans.begin();
    for (long val=lo; val < hi; ++val)
      *ptr++ = val;
    return ans;
  }


  bool IsDecimal(const std::istream& in) noexcept;
  bool IsDecimal(const std::ostream& out) noexcept;


  std::string fold(const std::string& str, long MaxLineLen);


  // This template fn should be used only for integral types T.
  // MaxSquarableInteger<T>() gives largest value whose square fits in T.
  template<typename T>
  T MaxSquarableInteger()
  {
    CoCoA_ASSERT(std::numeric_limits<T>::is_integer);
    CoCoA_ASSERT(std::numeric_limits<T>::is_bounded);
    const static T SqrtMax = FloorSqrt(std::numeric_limits<T>::max());
    return SqrtMax;
  }


/////  template <typename T> class IntegerCast; // fwd decl -- defined in convert.H

  // len is the same as .size() except that result is long rather than size_t
  template<typename T>
  long len(const T& container)
  {
#ifdef CoCoA_DEBUG
    return IntegerCast<long>(container.size());  // IntegerCast throws in case of overflow
#else
    return static_cast<long>(container.size());  // cast to long without checking
#endif
  }


} // end of namespace CoCoA

#endif
