#ifndef CoCoA_degree_H
#define CoCoA_degree_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti

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


// Header file for the class degree


#include "CoCoA/BigInt.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/utils.H"
// WARNING: next line must be last include (see redmine 591)
#include "CoCoA/VectorOps.H"  // for LexCmp3


#include <vector>
//using std::vector;
#include <iosfwd>
// using std::ostream;


namespace CoCoA
{

  // NOTE: most fns here are NOT noexcept because of arg checks

  namespace OrdvArith { class base; } // forward declaration for friend declaration

  class degree
  {
  public:
    explicit degree(long dim);
//???    friend long dim(const degree& deg);  // unnecessary because of GradingDim???
    const BigInt& operator[](long index) const;                      ///< deg[index] (READ ONLY), assumes 0 <= index < dim
    long myGradingDim() const noexcept;
    void mySetComponent(long index, const BigInt& N);                ///< deg[index] = N, assumes 0 <= index < dim
    void mySetComponent(long index, const MachineInt& n);            ///< deg[index] = n, assumes 0 <= index < dim
    degree& operator+=(const degree& d);                             ///< Computes deg1+=deg2
    degree& operator-=(const degree& d);                             ///< Computes deg1-=deg2
    friend int FastCmp(const degree& d1, const degree& d2) noexcept; ///< like cmp but inline and no checks
    friend inline bool IsZero(const degree& d) noexcept;
//    static void CheckCompatible(const degree& d1, const degree& d2, const char* fn); ///< checks dims are equal, throws if not.
    static void CheckCompatible(const degree& d1, const degree& d2, const ErrorContext& ErCtx); ///< checks dims are equal, throws if not.
  private: // data members
    std::vector<BigInt> myCoords;
  };


  long GradingDim(const degree& d) noexcept;                         ///< Computes "dimension" of the degree
  degree operator+(const degree& d1, const degree& d2);              ///< Computes d1+d2
  degree operator-(const degree& d1, const degree& d2);              ///< Computes d1-d2
  degree top(const degree& d1, const degree& d2);                    ///< i-th compt is max(d1[i], d2[i])
  int cmp(const degree& d1, const degree& d2);                       ///< lex comparison, result is <0,=0,>0
  void SetComponent(degree& d, long index, const BigInt& N);
  void SetComponent(degree& d, long index, const MachineInt& n);
  std::ostream& operator<<(std::ostream& out, const degree& d);      ///< prints out the degree, handles dim=1 specially.

  //---------------------------------------------------------------------------
  // inline functions

  inline degree::degree(long dim):
      myCoords(dim)
  {}


  inline long degree::myGradingDim() const noexcept
  {
    return len(myCoords);
  }


  // FastCmp is the same as cmp except it is inline and does not check compatibility
  inline int FastCmp(const degree& d1, const degree& d2) noexcept
  {
    const long dim = d1.myGradingDim();
    return LexCmp3(&d1.myCoords[0], &d1.myCoords[dim],
                   &d2.myCoords[0], &d2.myCoords[dim]);
  }


  inline bool operator<(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) < 0;
  }


  inline bool operator<=(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) <= 0;
  }


  inline bool operator>(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) > 0;
  }


  inline bool operator>=(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) >= 0;
  }


  inline bool operator==(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) == 0;
  }


  inline bool operator!=(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) != 0;
  }


  inline long GradingDim(const degree& d) noexcept
  {
    return d.myGradingDim();
  }

} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/degree.H,v 1.19 2024/03/16 10:15:11 abbott Exp $
// $Log: degree.H,v $
// Revision 1.19  2024/03/16 10:15:11  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.18  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.17  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.16  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.15  2018/05/18 12:13:37  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.14  2018/05/17 16:06:19  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.13  2017/01/13 16:05:42  abbott
// Summary: Added comment (see redmine 591)
//
// Revision 1.12  2014/07/31 14:45:16  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.11  2014/07/15 16:50:05  abbott
// Summary: Reordered include directives (o/w get weird bug with some g++ versions)
// Author: JAA
//
// Revision 1.10  2014/07/14 15:05:07  abbott
// Summary: Changed include after moving LexCmp3
// Author: JAA
//
// Revision 1.9  2013/06/03 15:27:50  bigatti
// -- added include IntOperations.H
// -- removed doxygen include
//
// Revision 1.8  2011/11/09 13:51:29  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.7  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.6  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
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
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/08 18:42:05  cocoa
// Cleaned up whitespace.
//
// Revision 1.4  2006/12/06 17:17:44  cocoa
// -- removed #include "config.H"
//
// Revision 1.3  2006/10/06 14:04:15  cocoa
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
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.5  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.4  2005/07/01 16:09:14  cocoa
// Degrees may now have negative components.
//
// Revision 1.3  2005/06/27 16:23:04  cocoa
// -- Added: GradingDim, operator+=, operator-=
// -- Removed: IsSubtractable (and calls to it in operator-)
//
// Revision 1.2  2005/05/04 16:44:57  cocoa
// -- changed: in FastCmp  ">" --> cmp
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.4  2005/04/29 14:21:17  cocoa
// -- changed: ElementType is now ZZ (was int)
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/28 15:58:56  cocoa
// Resynch after some minor changes.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.6  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.5  2004/11/11 13:02:58  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.4  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
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
// Revision 1.3  2003/06/23 16:22:46  abbott
// Minor cleaning prior to public release.
//
// Revision 1.2  2003/05/14 17:07:47  abbott
// Changed name of PPOrderingBase (following PPOrdering.H).
// Added (read-only) accessor function (via []).
// Added FastCmp (like cmp but no checks) with inline definition.
// Added top function (need a better name?).
// Added SetComponent private member function.
//
// Revision 1.1  2002/12/18 18:19:06  abbott
// Initial revision
//

#endif
