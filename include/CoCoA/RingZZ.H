#ifndef CoCoA_RingZZ_H
#define CoCoA_RingZZ_H

//   Copyright (c)  2005-2007,2010  John Abbott, and Anna M. Bigatti

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


//---------------------------------------------------------------------------
// RingZZ() is the CoCoA ring representing the integers.


namespace CoCoA
{

  class ring; // fwd decl, defined in ring.H.
  class RingHom; // fwd decl, defined in RingHom.H.

  const ring& RingZZ();    ///< The unique ring of integers.
  bool IsZZ(const ring& R);
  RingHom ZZEmbeddingHom(const ring& codomain);
  RingHom ZZEmbeddingHom(const ring& Z, const ring& codomain);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RingZZ.H,v 1.3 2022/02/18 14:11:48 abbott Exp $
// $Log: RingZZ.H,v $
// Revision 1.3  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.1  2012/02/10 10:23:12  bigatti
// -- was RingZ
//
// Revision 1.4  2012/02/08 13:37:35  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.3  2010/09/30 14:35:34  abbott
// Added new pseudo-ctor for ZEmbeddingHom which requires the domain as first arg:
// this form is essential for the ctor for RingQImpl!
//
// Added (hidden) new fn which checks whether RingZ has ref count greater than 3:
// this fn tells us whether any object other than the GlobalManager and RingQ
// refers to RingZ.  The fn is called only from the dtor of GlobalManager when
// checking that no CoCoALib values still survive.
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/05 21:22:26  cocoa
// Better comments, renamed the embedding hom pseudo-ctor.
//
// Revision 1.7  2007/03/03 14:07:23  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2006/12/06 17:21:55  cocoa
// -- removed #include "config.H"
//
// Revision 1.4  2006/11/20 15:55:03  cocoa
// ring is now a class again.  Improved definitions of operator-> in derived classes.
//
// Revision 1.3  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
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
// Revision 1.8  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.7  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.6  2006/03/21 09:43:14  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
//
// Revision 1.5  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.4  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.3  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.2  2005/12/31 12:22:18  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/09/30 15:03:39  cocoa
// Minor cleaning and tidying.
// DistrMPolyInlPP: use of summands now rather cleaner.
//
// Revision 1.2  2005/06/22 14:42:16  cocoa
// Renamed MemPool data member to myMemMgr
// (seems more sensible than myMemory).
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
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
// Revision 1.15  2004/11/18 18:33:40  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.14  2004/11/11 13:20:47  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.13  2004/11/04 18:47:42  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.12  2004/07/27 16:03:38  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.11  2004/07/16 15:45:12  cocoa
// First stage of new RingElem implementation completed.
//
// Revision 1.10  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.9  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.8  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.7  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.6  2004/01/30 14:07:10  cocoa
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
// Revision 1.5  2003/11/14 13:06:05  cocoa
// -- New function "myIsPrintAtom" for printing polynomials and fractions
//
// Revision 1.4  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.3  2003/10/09 14:55:19  cocoa
// - minor debugging after merge
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.15  2003/06/23 16:40:44  abbott
// Prior to public release:
//   changed class name to follow coding conventions
//   updated characteristic following RingBase
//   added IsIntegralDomain following RingBase
//
// Revision 1.14  2003/05/14 17:15:17  abbott
// Several consequential changes stemming from revisions to ring.H;
// most notably added support for homomorphisms and ideals.
//
// Revision 1.13  2002/11/15 11:23:27  abbott
// Updated following the renaming in ring.H.
// Decided to use the default implementation of IsZeroAddMul.
//
// Revision 1.12  2002/07/05 15:19:53  abbott
// Added IsDivisible member function.
// The zero element is now held through an auto_ptr (for exception cleanliness).
//
// Revision 1.11  2002/06/22 17:10:41  abbott
// Changed name of "equal" member function to "IsEqual" (as per new ring.H).
//
// Revision 1.10  2002/06/03 13:32:26  abbott
// Added compatibility with the new OrderedRing abstract class.
//
// Revision 1.9  2002/05/30 13:51:43  abbott
// Added IsGCDDomain, IsField, and zero functions as required by the new ring.H.
// Constructor and destructor are no longer inline (no real need to be, I hope).
//
// Revision 1.8  2002/05/15 15:05:01  abbott
// General cleaning of the source.
// Added negate, characteristic and power functions (req'd by new ring.H).
// Made myMemory data member private -- previously it was public?!?
//
// Revision 1.7  2002/03/21 15:12:21  bigatti
// - added IsZeroAddMul
// - corrected pool name
//
// Revision 1.6  2002/02/15 11:59:21  bigatti
// - added: gcd
// - removed (virtual ... const) = 0
//
// Revision 1.5  2001/12/07 18:21:34  abbott
// Changed names in accordance with the new coding conventions.
//
// Revision 1.4  2001/11/23 20:58:34  abbott
// Added assignment from a long.
//
// Revision 1.3  2001/11/19 20:09:10  abbott
// Added std:: prefix to identifiers from the standard library
// (required by g++ version 3.0.2, and also properly correct).
//
// Revision 1.2  2001/11/06 16:46:47  abbott
// Two main changes:
//
//  * alignment with the change to ring::raw_elem from void* to a union
//
//  * the elements of a ring_Z object are now pointers to mpz_t structs
//    (rather than objects of class ZZ); the memory for the head of the
//    mpz_t is managed by a MemPool belonging to the ring_Z.
//
// Revision 1.1  2001/10/31 20:37:12  abbott
// Initial revision
//

#endif
