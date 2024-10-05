#ifndef CoCoA_RingQQ_H
#define CoCoA_RingQQ_H

//   Copyright (c)  2003--2008  John Abbott, and Anna M. Bigatti

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


#include "CoCoA/FractionField.H"

// RingQQ() is the CoCoA ring representing the rationals.
namespace CoCoA
{

  class RingHom; // fwd decl, defined in RingHom.H

  const FractionField& RingQQ(); ///< The unique field of rationals.
  bool IsQQ(const ring& R);
  RingHom QQEmbeddingHom(const ring& codomain); ///< NB this can produce a *partial* homomorphism.

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RingQQ.H,v 1.3 2022/02/18 14:11:48 abbott Exp $
// $Log: RingQQ.H,v $
// Revision 1.3  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.1  2012/02/10 10:22:53  bigatti
// -- was RingQ
//
// Revision 1.5  2012/02/08 13:37:35  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.4  2009/09/24 15:45:13  abbott
// Improved two comments.
//
// Revision 1.3  2008/11/23 18:32:17  abbott
// Now includes FractionField.H, otherwise you cannot do much with a RingQ
// (without explicitly including FractionField.H yourself).
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/05 21:20:10  cocoa
// Better comments, renamed the embedding hom.
//
// Revision 1.6  2007/03/03 14:07:23  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.4  2007/02/26 15:00:54  bigatti
// -- added placeholders for new syntax based on unique Z implementation
//
// Revision 1.3  2006/12/06 17:21:55  cocoa
// -- removed #include "config.H"
//
// Revision 1.2  2006/10/06 10:15:53  cocoa
// In response to Susan's bug: a fiasco when compiling with CoCoA_MEMPOOL_DEBUG
// set wrongly.  Moved several implementation classes out of their header files
// into the implementation files.  Several functions had to be uninlined.
// Also corrected position of #include, etc.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.2  2005/06/22 14:42:16  cocoa
// Renamed MemPool data member to myMemMgr
// (seems more sensible than myMemory).
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.14  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.13  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.12  2004/11/11 13:04:30  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.11  2004/11/05 15:34:33  cocoa
// Consequential change following from the renaming of
// FieldIdealImpl and the introduction of the new pseudo-ctor.
//
// Revision 1.10  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.9  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.8  2004/07/16 15:45:12  cocoa
// First stage of new RingElem implementation completed.
//
// Revision 1.7  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.6  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.5  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.4  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.3  2004/01/30 14:07:10  cocoa
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
// Revision 1.2  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.1  2003/10/09 12:15:15  cocoa
// New coding convention for rings.
//

#endif
