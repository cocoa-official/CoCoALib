#ifndef CoCoA_FractionField_H
#define CoCoA_FractionField_H

//   Copyright (c)  2005-2007  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/ring.H"
#include "CoCoA/error.H"

namespace CoCoA
{

  class FractionFieldBase: public RingBase
  {
  public:
    FractionFieldBase(const ring& BaseRing): myBaseRingValue(BaseRing) {}
    // functions special to a FractionField
    virtual const ring& myBaseRing() const { return myBaseRingValue; }
    virtual ConstRawPtr myRawNum(ConstRawPtr rawq) const = 0;  ///< result belongs to BaseRing!!
    virtual ConstRawPtr myRawDen(ConstRawPtr rawq) const = 0;  ///< result belongs to BaseRing!!
    virtual RingElem mySymbolValue(const symbol& sym) const;   ///< returns the RingElem corresponding to sym
    virtual BigInt myCharacteristic() const;
    virtual bool IamCommutative() const;
    virtual bool3 IamIntegralDomain3(bool) const;
    virtual bool IamOrderedDomain() const;
    virtual bool IamField() const;
    virtual bool IamFiniteField() const;
    virtual bool IamExact() const;

    virtual RingHom myEmbeddingHomCtor() const = 0;
    virtual RingHom myInducedHomCtor(const RingHom& phi) const = 0;
    void myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const; ///< lhs = deriv(f, x)), only if BaseRing is PolyRing!!
  protected: // data member
    ring myBaseRingValue;
  };


  /*-----------------------------------------------------------------*/
  /** \include FractionField.txt  */
  /*-----------------------------------------------------------------*/

  class FractionField: public ring
  {
  public:
    FractionField(const ring& R);
    explicit FractionField(const FractionFieldBase* RingPtr);
    // Default copy ctor & assignment work fine.
    // Default dtor works fine.
  public:
    const FractionFieldBase* operator->() const; ///< allow const member fns to be called
  };


  // Pseudo-ctor and functions on FractionFields.
  FractionField NewFractionField(const ring& R);
  const FractionFieldBase* FractionFieldPtr(const ring& R);
///  const FractionFieldBase* FractionFieldPtr(const ring& R, const char* const FnName);
  const FractionFieldBase* FractionFieldPtr(const ring& R, const ErrorContext& ErrCtx);
  bool IsFractionField(const ring& Q);
  bool IsFractionFieldOfGCDDomain(const ring& Q);
  RingHom EmbeddingHom(const FractionField& FrF);
  RingHom InducedHom(const FractionField& FrF, const RingHom& InducingHom);

  // Functions on RingElems applicable only to elems of a FractionField.
  RingElem num(ConstRefRingElem q);
  RingElem den(ConstRefRingElem q);


  /////////////////////////////////////////////////////////////////////////////
  // Inline function

  inline bool IsFractionField(const ring& R)
  {
    return FractionFieldPtr(R) != nullptr;
  }

  inline const FractionFieldBase* FractionField::operator->() const
  {
    return static_cast<const FractionFieldBase*>(ring::operator->());
  }


  inline FractionField::FractionField(const ring& R):
      ring(FractionFieldPtr(R,CoCoA_ERROR_CONTEXT))
  {}

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/FractionField.H,v 1.18 2024/03/15 19:04:45 abbott Exp $
// $Log: FractionField.H,v $
// Revision 1.18  2024/03/15 19:04:45  abbott
// Summary: Arg checking fns now expect a ErrorContext object (instead of string; related to redmine 1726)
//
// Revision 1.17  2023/03/08 19:43:56  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.16  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.15  2021/01/07 14:13:40  abbott
// Summary: Corrected copyright
//
// Revision 1.14  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.13  2016/04/07 14:35:45  abbott
// Summary: Removed code which disabled assignment of rings
//
// Revision 1.12  2014/07/28 16:02:58  abbott
// Summary: Renamed myEmbeddingHom to myEmbeddingHomCtor
// Author: JAA
//
// Revision 1.11  2014/07/28 15:09:23  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble)
// Author: JAA
//
// Revision 1.10  2014/07/09 11:36:28  abbott
// Summary: Tidying: moved several simple fns to FractionalFieldBase
// Author: JAA
//
// Revision 1.9  2014/07/08 15:20:39  abbott
// Summary: Changed myBaseRingValue to protected; gave common defn of myBaseRing
// Author: JAA
//
// Revision 1.8  2014/07/08 08:42:25  abbott
// Summary: Removed AsFractionField, added FractionFieldPtr
// Author: JAA
//
// Revision 1.7  2013/05/21 11:57:25  abbott
// Replaced mem fns myGetNum (&Den) by myRawNum (&Den).
//
// Revision 1.6  2013/05/14 14:21:31  abbott
// Revised/improved impl of derivative of ratfns.
//
// Revision 1.5  2012/05/22 10:02:38  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.4  2010/10/01 15:20:33  bigatti
// -- added mySymbolValue
// -- added RingElem(R, sym)
//
// Revision 1.3  2009/09/24 14:36:03  abbott
// Corrected a comment.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.7  2007/01/11 14:06:15  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.6  2006/12/06 17:11:16  cocoa
// -- removed #include "config.H"
//
// Revision 1.5  2006/11/20 15:55:03  cocoa
// ring is now a class again.  Improved definitions of operator-> in derived classes.
//
// Revision 1.4  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
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
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.3  2006/05/12 16:10:59  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/07/08 15:09:29  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.2  2005/06/22 14:42:16  cocoa
// Renamed MemPool data member to myMemMgr
// (seems more sensible than myMemory).
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.4  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
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
// Revision 1.16  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.15  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.14  2004/11/08 14:03:32  cocoa
// -- small changes for doxygen support
//
// Revision 1.13  2004/11/05 15:34:33  cocoa
// Consequential change following from the renaming of
// FieldIdealImpl and the introduction of the new pseudo-ctor.
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
// Revision 1.10  2004/07/16 15:45:12  cocoa
// First stage of new RingElem implementation completed.
//
// Revision 1.9  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.8  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.7  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.6  2004/02/03 16:16:21  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.5  2004/01/30 14:07:10  cocoa
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
// Revision 1.4  2003/11/14 13:06:05  cocoa
// -- New function "myIsPrintAtom" for printing polynomials and fractions
//
// Revision 1.3  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.6  2003/06/23 16:45:30  abbott
// Prior to public release:
//   cleaned up num/den member functions
//   revised GeneralFractionField following the new RingBase class
//
// Revision 1.5  2003/05/27 15:53:27  abbott
// Changed name from FractionField to GeneralFractionField.
// Added EmbeddingHom and EmbeddingHom_t (to map from R into FrF(R)).
//
// Revision 1.4  2003/05/15 10:15:58  abbott
// Consequential changes from ring.H: old changes to names and some
// member functions.  Added homomorphisms and ideals (reqd by ring.H).
//
// Revision 1.3  2002/11/15 10:41:25  abbott
// Revised according the renaming in ring.H.
// A small amount of tidying up too.
//
// Revision 1.2  2002/07/05 15:31:13  abbott
// Added member function IsDivisible.
// Zero element of ring is now in an auto_ptr (possibly gives better exception safety).
//
// Revision 1.1  2002/06/22 16:47:02  abbott
// Initial revision
//
//

#endif
