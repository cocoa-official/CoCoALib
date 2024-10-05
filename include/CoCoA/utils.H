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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/utils.H,v 1.37 2024/07/31 07:52:47 abbott Exp $
// $Log: utils.H,v $
// Revision 1.37  2024/07/31 07:52:47  abbott
// Summary: Switched to new CoCoA_THROW_ERROR macros (with suffix 1 or 2)
//
// Revision 1.36  2024/05/03 12:07:55  abbott
// Summary: changed spacing
//
// Revision 1.35  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.34  2021/09/01 08:15:56  abbott
// Summary: Added comment
//
// Revision 1.33  2021/03/12 16:50:12  abbott
// Summary: Minor change to comment
//
// Revision 1.32  2021/02/11 23:04:30  abbott
// Summary: Added temporary templates CMP and IsInRange_
//
// Revision 1.31  2021/02/10 19:32:51  abbott
// Summary: Added noexcept (redmine 1572), and added cmp fn for long and ulong
//
// Revision 1.30  2021/01/07 14:48:34  abbott
// Summary: Corrected copyright
//
// Revision 1.29  2020/12/05 13:05:10  abbott
// Summary: Added IsDecimal for istream and ostream
//
// Revision 1.28  2020/06/17 15:49:21  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.27  2020/02/03 17:04:47  abbott
// Summary: Added new fn fold
//
// Revision 1.26  2017/12/12 14:18:02  abbott
// Summary: Corrected use of CoCoA_DEBUG CPP symbol
//
// Revision 1.25  2017/02/16 11:56:02  abbott
// Summary: Improved comment
//
// Revision 1.24  2017/02/15 12:23:11  abbott
// Summary: Removed square fn
//
// Revision 1.23  2015/11/24 12:47:04  abbott
// Summary: Renamed "isqrt" --> "FloorSqrt"
//
// Revision 1.22  2014/07/14 15:05:35  abbott
// Summary: Removed LexCmp3 (moved into UtilsTemplate.H)
// Author: JAA
//
// Revision 1.21  2014/04/30 16:03:08  abbott
// Summary: Modified len to enable range check only when CoCoA_DEBUG is defined
// Author: JAA
//
// Revision 1.20  2014/04/11 13:36:19  abbott
// Summary: Added missing std:: prefix in two assertions
// Author: JAA
//
// Revision 1.19  2014/04/04 10:14:44  abbott
// Summary: New shorter impl of MaxSquarableInteger
// Author: JAA
//
// Revision 1.18  2014/04/03 15:52:15  abbott
// Summary: Improved commented out impl of BiggestSqBelow2Power
// Author: JAA
//
// Revision 1.17  2014/04/03 15:32:35  abbott
// Summary: New impls for LongRange, ULongDiff, MaxSquarableInteger
// Author: JAA
//
// Revision 1.16  2013/02/14 17:35:14  bigatti
// -- added LongRange -- useful for submat calls
//
// Revision 1.15  2011/12/05 16:54:40  bigatti
// -- changed: MachineInteger --> MachineInt (just in comment)
//
// Revision 1.14  2011/09/06 13:37:50  abbott
// Moved ULong2Long into separate file (ULong2Long.H)
//
// Revision 1.13  2011/08/27 21:48:11  abbott
// Commented out template fn "cmp" -- because it preempts cmp fn for MachineInt (qv.)
//
// Revision 1.12  2011/08/27 20:49:20  abbott
// Improved defn of ULong2Long -- needs more work.
//
// Revision 1.11  2011/03/14 10:25:52  abbott
// Added ULong2Long: an inline fn to convert unsigned long into long
// (i.e. invert the automatic conversion long to unsigned long).
//
// Revision 1.10  2011/03/09 17:21:21  bigatti
// -- added #include "CoCoA/convert.H"
//
// Revision 1.9  2011/03/09 15:52:52  abbott
// Added fwd decl for NumericCast.
//
// Revision 1.8  2011/03/09 14:56:19  abbott
// Added template fn "len" (v. similar to mem fn size).
//
// Revision 1.7  2009/10/29 18:30:40  abbott
// Corrected copyright years.
//
// Revision 1.6  2009/06/04 14:39:29  abbott
// Added template class NumericLimits (similar to std::numeric_limits) so that
// one can easily find the largest integer whose square fits into a given type.
//
// Revision 1.5  2008/11/06 12:50:44  abbott
// Moved definitions of square and round to utils.H from ApproxPts.H
//
// Revision 1.4  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/05/21 12:53:08  abbott
// No real change in the end.
//
// Revision 1.2  2007/03/23 18:38:43  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/08 18:42:05  cocoa
// Cleaned up whitespace.
//
// Revision 1.4  2006/12/06 17:18:40  cocoa
// -- removed #include "config.H"
//
// Revision 1.3  2006/11/27 13:06:23  cocoa
// Anna and Michael made me check without writing a proper message.
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//

#endif
