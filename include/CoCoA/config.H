#ifndef CoCoA_config_H
#define CoCoA_config_H

//   Copyright (c)  2005,2013  John Abbott,  Anna M. Bigatti

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


#include <climits>

namespace CoCoA
{

  // >>> EXPONENTS <<<

  // This is the type used internally for "small" exponents in some PPMonoids
  // (PPMonoidEv, PPMonoidEvOv).  Also implies the number of bits used for a
  // "compressed entry" in an order vector (see OrdvArith).
  typedef unsigned int SmallExponent_t;


  // -------------------------------------------------------
  // >>> ROUNDING HALVES <<<

  // The way CoCoALib rounds halves is determined by the fns RoundDiv in
  // IntOperations.C: one for machine ints, and one for BigInts.


  // -------------------------------------------------------
  // >>> INTEGER FACTORIZATION <<<

  // Trial division limit for factor(BigInt) defined in NumTheory-factor.C
  constexpr long FactorBigIntTrialLimit = 33554432; // arbitrary limit (happens to be 2^25)


  // -------------------------------------------------------
  // >>> INTEGER OVERFLOW <<<

  // Some functions check for BigInt overflow: power & factorial.
  // This constant sets the max bit-size allowed in result value.
  // Need to check whether platform is 32-bit (hence #if construct).
#if LONG_MAX/(1024L*1024L*1024L) >= 4
  constexpr long OVERFLOW_BITS = 4*1024L*1024L*1024L; // more than 32 bits
#else
  constexpr long OVERFLOW_BITS = 1024L*1024L*1024L; // 32-bitter
#endif

} // end of namespace CoCoA



// RCS header/log in the next few lines.
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/config.H,v 1.14 2022/02/18 14:11:51 abbott Exp $
// $Log: config.H,v $
// Revision 1.14  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.13  2021/09/29 16:39:38  abbott
// Summary: Modified to work on 32-bit platforms; improved comments
//
// Revision 1.12  2021/09/22 13:58:36  abbott
// Summary: Increased OVERFLOW_BITS, improved comment
//
// Revision 1.11  2021/09/13 15:41:02  abbott
// Summary: Added constexpr OVERFLOW_BITS
//
// Revision 1.10  2021/02/10 19:35:18  abbott
// Summary: Changed FactorBigIntTrialLimit to constexpr (redmine 1511)
//
// Revision 1.9  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.8  2018/05/05 15:21:42  abbott
// Summary: Increased FactorBigIntTrialLimit
//
// Revision 1.7  2014/09/01 16:23:04  abbott
// Summary: Added new const FactorBigIntTrialLimit
// Author: JAA
//
// Revision 1.6  2014/05/16 12:31:25  abbott
// Summary: Added comment about rounding convention (just refers to fn RoundDiv)
// Author: JAA
//
// Revision 1.5  2013/05/27 14:07:00  abbott
// Removed typedefs SmallFpElem_t & SmallFpLogElem_t; they are
// now nested typedefs in SmallFpImpl & SmallFpLogImpl.
// Removed PPExpSize enum; it now resides in PPMonoidEv.H.
//
// Revision 1.4  2013/02/26 16:49:05  abbott
// Changed SmallExponent_t into "unsigned int" (as suggested on redmine).
// Added useful comment about OrdvArith too.
//
// Revision 1.3  2013/02/21 14:10:53  abbott
// Added enum for choosing between "small" and "big" exponents in PPMonoidElems
// (to be used by various PPMonoid pseudo-ctors).
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
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
// Revision 1.2  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.4  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.3  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.2  2005/07/08 15:09:29  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.5  2004/11/08 11:59:33  cocoa
// Messages produced by a failed CoCoA_ASSERT are now printed on
// GlobalErrput() -- previously they were sent to cerr.
//
// Revision 1.4  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//
// Revision 1.3  2004/07/13 16:32:26  cocoa
// First stage of major revamp of ring elements.
// Implementation of RingFp has been split into "ring interface"
// and "algorithms plus data structures".
//
// Revision 1.2  2004/01/30 14:07:10  cocoa
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
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.6  2002/11/13 14:39:27  abbott
// Added a non-empty definition of the CoCoA_ASSERT macro for the
// non-debugging case.  Minor simplification to the debugging definition.
//
// Revision 1.5  2002/11/12 17:07:24  abbott
// Changed name of ASSERT macro into CoCoA_ASSERT (in accordance with
// the coding conventions).
//
// Revision 1.4  2002/03/28 15:59:02  abbott
// Revised ASSERT macro to call AssertionFailed function.
// Added AssertionFailed function to simplify debugging when an
// assertion fails.
//
// Revision 1.3  2001/12/07 18:08:13  abbott
// Removed IFDEBUG macro (and friends) because the GNU CPP rejects it.
// Added ASSERT macro.
//
// Revision 1.2  2001/12/07 15:52:50  abbott
// Added IFDEBUG macro (via the COMMENT and SLASH macros).
//
// Revision 1.1  2001/11/26 19:14:45  abbott
// Initial revision
//

#endif
