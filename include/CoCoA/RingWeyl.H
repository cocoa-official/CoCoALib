#ifndef CoCoA_RingWeyl_H
#define CoCoA_RingWeyl_H

//   Copyright (c)  2003-2007, 2011  John Abbott, and Anna M. Bigatti

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


#include "CoCoA/SparsePolyRing.H"

#include <vector>
//using std::vector;

namespace CoCoA
{

  /*-----------------------------------------------------------------*/
  /** \include RingWeyl.txt  */
  /*-----------------------------------------------------------------*/

  // Pseudo-ctors for (sparse) polynomial rings.

  std::vector<symbol> WANAMES(const std::vector<symbol>& names); ///< {x,y,z,..} --> {x,y,z,..,dx,dy,dz,..}

  SparsePolyRing NewWeylAlgebra(const ring& CoeffRing, long NumTrueIndets, const std::vector<long>& ElimIndets);
  SparsePolyRing NewWeylAlgebra(const ring& CoeffRing, const std::vector<symbol>& names, const std::vector<long>& ElimIndets);
//   const RingElem& indet(const RingWeyl& RW, long var);
//   const RingElem& derivation(const RingWeyl& RW, long var);


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RingWeyl.H,v 1.6 2022/02/18 14:11:48 abbott Exp $
// $Log: RingWeyl.H,v $
// Revision 1.6  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.3  2011/03/08 17:57:21  bigatti
// -- changed: ElimIndets list is now of  long  instead of size_t
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2006/12/06 17:21:55  cocoa
// -- removed #include "config.H"
//
// Revision 1.7  2006/11/24 17:25:00  cocoa
// -- reorganized includes of header files
//
// Revision 1.6  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.5  2006/08/17 09:41:33  cocoa
// -- changed:  I think it works again
//
// Revision 1.4  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.3  2006/07/20 14:11:27  cocoa
// -- more stable version (more similar to RingDistrMPoly)
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.6  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.5  2006/03/21 09:43:14  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
//
// Revision 1.4  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.3  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.2  2006/02/13 14:46:45  cocoa
// -- changes by Max
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.2  2005/05/04 16:47:45  cocoa
// -- updated with new function/class names
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.13  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.12  2004/11/11 14:27:17  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.11  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.10  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//
// Revision 1.9  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.8  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.7  2004/06/17 12:44:58  cocoa
// -- applied new coding conventions: compiles and seems to work
//
// Revision 1.6  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.5  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.4  2004/01/30 14:07:10  cocoa
// Tidied RingRawValue union: now it contains just two fields,
// and has no need of forward declarations of types used internally
// by the concrete rings -- it uses explicitly a void* instead.
//
// I have tidied the "import" functions used by most concrete rings.
//
// I have moved the choice of representation type for RingFp and RingFpLog
// into a typedef in config.H -- this is to recognise that different
// choices may work best on different platforms.
//
// Revision 1.3  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.3  2003/06/23 16:58:25  abbott
// Minor cleaning prior to public release.
// Just consequential changes.
//
// Revision 1.2  2003/05/30 15:23:28  abbott
// Changed beyond all recognition -- completely new version.
//
// Revision 1.1  2003/05/02 13:07:48  abbott
// Initial revision
//

#endif
