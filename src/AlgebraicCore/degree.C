//   Copyright (c)  2002-2009,2011  John Abbott and Anna M. Bigatti

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

#include "CoCoA/degree.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
#include<algorithm>
using std::max;
//#include <vector>
using std::vector;

namespace CoCoA
{

  // inline void degree::CheckCompatible(const degree& d1, const degree& d2, const char* fn)
  // {
  //   if (len(d1.myCoords) != len(d2.myCoords))
  //     CoCoA_THROW_ERROR(ERR::MixedDegrees, fn);
  // }

  inline void degree::CheckCompatible(const degree& d1, const degree& d2, const ErrorContext& ErrCtx)
  {
    if (len(d1.myCoords) != len(d2.myCoords))
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::MixedDegrees, ErrCtx);
  }


  const BigInt& degree::operator[](long index) const
  {
    if (index < 0 || index >= len(myCoords))
      CoCoA_THROW_ERROR(ERR::BadIndex, "degree::operator[]");
    return myCoords[index];
  }


  void degree::mySetComponent(long index, const BigInt& N)
  {
    CoCoA_ASSERT(0 <= index && index < len(myCoords));
    myCoords[index] = N;
  }


  void degree::mySetComponent(long index, const MachineInt& n)
  {
    CoCoA_ASSERT(0 <= index && index < len(myCoords));
    myCoords[index] = BigInt(n);
  }


  degree& degree::operator+=(const degree& d)
  {
    degree::CheckCompatible(*this, d, CoCoA_ERROR_CONTEXT);
    const long dim = len(myCoords);
    for (long i=0; i < dim; ++i)
      myCoords[i] += d[i];
    return *this;
  }


  degree& degree::operator-=(const degree& d)
  {
    degree::CheckCompatible(*this, d, CoCoA_ERROR_CONTEXT);
    const long dim = len(myCoords);
    for (long i=0; i < dim; ++i)
      myCoords[i] -= d[i];
    return *this;
  }


  bool IsZero(const degree& d) noexcept
  {
    const long dim = GradingDim(d);
    //??? return find_if(coords.begin(), coords.end(), NonZero) != coords.end();
    for (long i=0; i < dim; ++i)
      if (d[i] != 0) return false;
    return true;
  }


  degree operator+(const degree& d1, const degree& d2)
  {
    degree::CheckCompatible(d1, d2, CoCoA_ERROR_CONTEXT);
    const long dim = GradingDim(d1);
    degree ans(dim);
    for (long i=0; i < dim; ++i)
      ans.mySetComponent(i, d1[i] + d2[i]);
    return ans;
  }


  degree operator-(const degree& d1, const degree& d2)
  {
    degree::CheckCompatible(d1, d2, CoCoA_ERROR_CONTEXT);
    const long dim = GradingDim(d1);
    degree ans(dim);
    for (long i=0; i < dim; ++i)
      ans.mySetComponent(i, d1[i] - d2[i]);
    return ans;
  }


  degree top(const degree& d1, const degree& d2)
  {
    degree::CheckCompatible(d1, d2, CoCoA_ERROR_CONTEXT);
    const long dim = GradingDim(d1);
    degree ans(dim);
    for (long i=0; i < dim; ++i)
      ans.mySetComponent(i, max(d1[i], d2[i]));
    return ans;
  }


  int cmp(const degree& d1, const degree& d2)
  {
    degree::CheckCompatible(d1, d2, CoCoA_ERROR_CONTEXT);
    // The rest is the same as the body of FastCmp; writing it explicitly
    // here avoids compiler foibles (i.e. choosing not to make FastCmp inline).
//     const long dim = GradingDim(d1);
    return FastCmp(d1, d2);
//     return LexCmp3(&d1.myCoords[0], &d1.myCoords[dim],
//                    &d2.myCoords[0], &d2.myCoords[dim]);
//     for (long i=0; i < dim; ++i)
//       if (d1[i] != d2[i])
//         return (d1[i] > d2[i] ? 1 : -1);
//     return 0;
  }


  ostream& operator<<(ostream& out, const degree& d)
  {
    if (!out) return out;  // short-cut for bad ostreams
    const long dim = GradingDim(d);
    if (dim == 0) return out << "()";     // no grading -- graded over N^0
    //  if (dim == 1) return out << d[0]; // omit parens if grading is over N
    // General case
    out << "(";
    for (long i=0; i < dim-1; ++i)
      out << d[i] << ", ";
    out << d[dim-1] << ")";
    return out;
  }

  void SetComponent(degree& d, long index, const BigInt& N)
  {
    if (index < 0 || index >= GradingDim(d))
      CoCoA_THROW_ERROR(ERR::BadIndex, "SetComponent(degree, index, N)");
    d.mySetComponent(index, N);
  }

  void SetComponent(degree& d, long index, const MachineInt& n)
  {
    if (index < 0 || index >= GradingDim(d))
      CoCoA_THROW_ERROR(ERR::BadIndex, "SetComponent(degree, index, n)");
    d.mySetComponent(index, n);
  }


} // end of namespace CoCoA

// RCS header/log
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/degree.C,v 1.20 2024/07/22 15:43:58 abbott Exp $
// $Log: degree.C,v $
// Revision 1.20  2024/07/22 15:43:58  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.19  2024/07/03 09:44:11  bigatti
// Summary: error codes (redmine 92):
// removed BadIntedIndex, BadDegIndex, BadComptIndex
//
// Revision 1.18  2024/03/16 10:15:12  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.17  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.16  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.15  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.14  2020/06/17 15:49:29  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.13  2018/05/18 12:25:54  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.12  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.11  2014/07/14 15:09:34  abbott
// Summary: Added include of utils.H
// Author: JAA
//
// Revision 1.10  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.9  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.8  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.7  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.6  2009/02/20 09:55:43  bigatti
// -- changed: now also degrees with GradingDim=1 are printed with parentheses
//
// Revision 1.5  2008/12/17 11:53:08  abbott
// Indexes into degree objects are now MachineInts rather than size_t.
// Hope this proves not to be a duff idea.
//
// Revision 1.4  2008/12/16 21:36:41  abbott
// Changed long into MachineInt for specifying values when setting components.
//
// Revision 1.3  2008/04/21 12:32:54  abbott
// Corrected size_t into std::size_t in several header files; in some cases,
// replaced size_t with MachineInt (with consequent changes to impl files).
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.2  2005/06/27 16:23:04  cocoa
// -- Added: GradingDim, operator+=, operator-=
// -- Removed: IsSubtractable (and calls to it in operator-)
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.5  2005/04/29 14:34:14  cocoa
// -- changed: in mySetComponent  myCoords[pos] = ZZ(value);  (was "= value")
//
// Revision 1.4  2005/04/29 14:21:17  cocoa
// -- changed: ElementType is now ZZ (was int)
//
// Revision 1.3  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.2  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.5  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.4  2004/11/11 14:08:19  cocoa
// -- minor changes for doxygen
// -- moved CVS log to the bottom
//
// Revision 1.3  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.2  2004/01/28 15:29:42  cocoa
// Added IsZero function, and changed name of coordinate type used
// in the representation of degree objects.
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.3  2003/06/23 16:23:03  abbott
// Minor cleaning prior to public release.
//
// Revision 1.2  2003/05/14 17:09:08  abbott
// Added definitions of operator[] (read-only access),
// top and SetComponent.
//
// Revision 1.1  2002/12/18 18:19:20  abbott
// Initial revision
//
