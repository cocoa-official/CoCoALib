#ifndef CoCoA_QuotientRing_H
#define CoCoA_QuotientRing_H

//   Copyright (c)  2003-2009  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/ideal.H"
#include "CoCoA/ring.H"

#include <string>
// #include <vector>  // included by ring.H


// Header file for the abstract class QuotientRingBase
// and associated class QuotientRing
namespace CoCoA
{

  class RingHom;        // fwd decl -- defined in RingHom.H


  class QuotientRing; // fwd declaration needed for friend declarations

  class QuotientRingBase: public RingBase
  {
  protected:
    QuotientRingBase(const ring& R, const ideal& I);

  public:
    const ring& myBaseRing() const override; // inline defn
    // functions special to a QuotientRing -- all concrete classes must impl these
    virtual RingElem myCanonicalRepr(ConstRawPtr rawx) const = 0; ///< result is element of myReprRing
    virtual void myReduction(RawPtr rawimage, ConstRawPtr rawarg) const = 0;
    virtual RingHom myQuotientingHomCtor() const;
    virtual RingHom myInducedHomCtor(const RingHom& InducingHom) const = 0;
    RingElem mySymbolValue(const symbol& sym) const override;     ///< returns the RingElem corresponding to sym
  protected: // data members common to all implementations of quotient rings
    ring myBaseRingValue;
    const ideal myDefiningIdeal;
    const ring myReprRing;
    const ideal myReducingIdeal;
    friend const ideal& DefiningIdeal(const QuotientRing& RmodI);
    friend const ring& ReprRing(const QuotientRing& RmodI);
    friend const ideal& ReducingIdeal(const QuotientRing& RmodI);
  };


  class QuotientRing: public ring
  {
  public:
    QuotientRing(const ring& R);
    explicit QuotientRing(const QuotientRingBase* RingPtr);
    // Default copy ctor & assignment work fine.
    // Default dtor works fine.
  public:
    const QuotientRingBase* operator->() const; ///< covariant return type -- allow member fns to be called
  };


  QuotientRing NewQuotientRing(const ring& R, const ideal& I);
  QuotientRing NewQuotientRing(const ring& R,const std::string& str);
  QuotientRing NewQuotientRing(const ring& R,const std::vector<std::string>& L);
  QuotientRing NewZZmod(const MachineInt& n);
  QuotientRing NewZZmod(const BigInt& N);
  bool IsQuotientRing(const ring& R);
  const ideal& DefiningIdeal(const QuotientRing& RmodI);
  const ring& ReprRing(const QuotientRing& RmodI);
  const ideal& ReducingIdeal(const QuotientRing& RmodI);
  RingHom QuotientingHom(const QuotientRing& RmodI);
  RingHom InducedHom(const QuotientRing& RmodI, const RingHom& phi);
  RingElem CanonicalRepr(ConstRefRingElem r); ///< given x in R/I produces X in R s.t. QuotientingHom(RmodI)(X)=x

  //---------------------------------------------------------------------------

  inline const QuotientRingBase* QuotientRing::operator->() const
  {
    return static_cast<const QuotientRingBase*>(ring::operator->());
  }

  inline const QuotientRingBase* QuotientRingPtr(const ring& R)
  {
    return dynamic_cast<const QuotientRingBase*>(R.myRawPtr());
  }

  inline bool IsQuotientRing(const ring& R)
  {
    return QuotientRingPtr(R) != nullptr;
  }





  /////////////////////////////////////////////////////////////////////////////
  // Inline functions (the order is important).

  inline const ring& QuotientRingBase::myBaseRing() const
  {
    return myBaseRingValue;
  }


  inline const ideal& DefiningIdeal(const QuotientRing& RmodI)
  {
    return RmodI->myDefiningIdeal;
  }


  inline const ring& ReprRing(const QuotientRing& RmodI)
  {
    return RmodI->myReprRing;
  }


  inline const ideal& ReducingIdeal(const QuotientRing& RmodI)
  {
    return RmodI->myReducingIdeal;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/QuotientRing.H,v 1.22 2022/02/18 14:11:48 abbott Exp $
// $Log: QuotientRing.H,v $
// Revision 1.22  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.21  2021/10/30 19:03:08  abbott
// Summary: Added some more "override" (clang said they were mising; redmine 1625)
//
// Revision 1.20  2021/10/30 17:15:19  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.19  2021/01/07 14:22:33  abbott
// Summary: Corrected copyright
//
// Revision 1.18  2019/10/08 13:00:53  bigatti
// -- added NewQuotientRing with string
//
// Revision 1.17  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.16  2016/04/07 14:35:45  abbott
// Summary: Removed code which disabled assignment of rings
//
// Revision 1.15  2015/06/25 14:21:02  abbott
// Summary: Added minor comment (fwd ref to inline defn)
// Author: JAA
//
// Revision 1.14  2014/07/28 16:03:33  abbott
// Summary: Renamed myQuotientingHom to myQuotientingHomCtor
// Author: JAA
//
// Revision 1.13  2014/07/28 15:10:52  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble)
// Author: JAA
//
// Revision 1.12  2014/07/08 15:21:39  abbott
// Summary: Removed old defn of BaseRing
// Author: JAA
//
// Revision 1.11  2014/07/08 13:14:40  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.10  2012/02/10 12:01:07  bigatti
// -- removed useless fwd decl
// -- removed old doxygen comment
//
// Revision 1.9  2012/02/08 13:37:35  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.8  2011/11/09 13:47:56  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.7  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.6  2010/10/01 15:20:33  bigatti
// -- added mySymbolValue
// -- added RingElem(R, sym)
//
// Revision 1.5  2009/09/24 15:35:17  abbott
// Removed pointless include directives.
//
// Revision 1.4  2009/07/24 14:19:23  abbott
// Cleaned up include directives, and added fwd decls to the header file.
//
// Revision 1.3  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.15  2007/03/08 10:23:29  bigatti
// -- CanonHom --> CanonicalHom
//
// Revision 1.14  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.13  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (think of a better name)
//
// Revision 1.12  2007/02/26 15:00:54  bigatti
// -- added placeholders for new syntax based on unique Z implementation
//
// Revision 1.11  2007/01/11 14:06:15  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.10  2006/12/18 15:16:05  cocoa
// Added definition of CanonRepr.
//
// Revision 1.9  2006/12/06 17:15:26  cocoa
// -- removed #include "config.H"
//
// Revision 1.8  2006/11/24 17:34:49  cocoa
// -- reorganized includes of header files
//
// Revision 1.7  2006/11/20 15:55:03  cocoa
// ring is now a class again.  Improved definitions of operator-> in derived classes.
//
// Revision 1.6  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
//
// Revision 1.5  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.4  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
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
// Revision 1.6  2006/04/21 15:01:36  cocoa
// Changed default implementation of RingBase::myGcd -- it now gives a SERIOUS
// error.  All fields must now handle a call to gcd explicitly: they can use
// the new myGcdInField function.  It's now cleaner than it was.
//
// Revision 1.5  2006/03/21 09:43:14  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
//
// Revision 1.4  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.3  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
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
// Revision 1.5  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.3  2005/04/01 16:18:20  cocoa
// Friday check-in.  Fixed a bug in the ctor for GeneralQuotientRingImpl.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.17  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.16  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.15  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.14  2004/11/11 14:26:51  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.13  2004/11/09 15:50:54  cocoa
// -- changed GeneralQuotientRing --> GeneralQuotientRingImpl
// -- added: NewZmod(Z, N)  with N in ZZ
//
// Revision 1.12  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.11  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.10  2004/06/30 16:46:06  cocoa
// End of day check-in: mostly tidying up, and correcting some
// overly lax access permissions.
//
// Revision 1.9  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.8  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.7  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.6  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.5  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.4  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.3  2003/10/09 14:55:20  cocoa
// - minor debugging after merge
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/06/23 16:47:33  abbott
// Minor cleaning prior to public release.
// Updated characteristic and IsIntegralDomain following RingBase.
//
// Revision 1.1  2003/05/14 17:11:54  abbott
// Initial revision
//

#endif
