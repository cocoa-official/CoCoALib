#ifndef CoCoA_PPMonoidEvOv_H
#define CoCoA_PPMonoidEvOv_H

//   Copyright (c)  2005,2007,2010  John Abbott,  Anna M. Bigatti

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


#include <vector>
// using vector;

namespace CoCoA
{

  class symbol;         // forward declaration -- defined in symbol.H
  class PPMonoid;       // forward declaration -- defined in PPMonoid.H
  class PPOrdering;     // forward declaration -- defined in PPOrdering.H
  class PPOrderingCtor; // forward declaration -- defined in PPOrdering.H

  // Pseudo-ctor
  PPMonoid NewPPMonoidEvOv(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
  PPMonoid NewPPMonoidEvOv(const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord);

} // end of namespace CoCoA

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/PPMonoidEvOv.H,v 1.6 2022/02/18 14:11:47 abbott Exp $
// $Log: PPMonoidEvOv.H,v $
// Revision 1.6  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2021/01/07 14:20:19  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.3  2010/02/03 16:13:52  abbott
// Added new single word tags for specifying the ordering in PPMonoid
// pseudo-ctors.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.4  2006/11/27 13:41:32  cocoa
// -- added #include <cstddef> //using std::size_t
// -- removed #include "config.H"
//
// Revision 1.3  2006/11/24 17:31:33  cocoa
// -- reorganized includes of header files
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.5  2006/03/14 17:21:18  cocoa
// Moved concrete PPMonoid impls entirely into their respective .C files.
// Now the corresponding .H files are very compact.
//
// Revision 1.4  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.3  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.2  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.8  2005/10/11 16:37:31  cocoa
// Added new small prime finite field class (see RingFpDouble).
//
// Cleaned makefiles and configuration script.
//
// Tidied PPMonoid code (to eliminate compiler warnings).
//
// Fixed bug in RingFloat::myIsInteger.
//
// Revision 1.7  2005/08/08 16:36:33  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.6  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.5  2005/07/08 15:09:29  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.4  2005/07/01 16:08:16  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.3  2005/06/23 15:42:41  cocoa
// Fixed typo in GNU fdl -- all doc/*.txt files affected.
// Minor corrections to PPMonoid (discovered while writing doc).
//
// Revision 1.2  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from ZZ).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.5  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.4  2004/11/11 13:42:15  cocoa
// -- minor changes for doxygen
//
// Revision 1.3  2004/11/11 11:56:09  cocoa
// (1) Tidied makefiles, and introduced common.mki
// (2) Improved several tests, and cleaned them so that they
//     handle sanely any otherwise unhandled exceptions.
//
// Revision 1.2  2004/11/02 14:56:33  cocoa
// -- changed *Print* into *Output* (myPrint --> myOutput)
// -- changed *Var* into *Indet* (myPrintVarName --> myOutputIndetName)
// -- removed suffix "IgnoreDivMask"
// -- added myComputeDivMask
// -- improved storing of IndetNames
// -- changed ExpvElem into SmallExponent_t
//
// Revision 1.1  2004/10/29 15:31:25  cocoa
// -- new PPMonoid for compatibility with OrdvArith (without DivMask)
//

#endif
