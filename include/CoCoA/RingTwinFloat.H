#ifndef CoCoA_RingTwinFloat_H
#define CoCoA_RingTwinFloat_H

//   Copyright (c)  2004-2008  John Abbott, and Anna M. Bigatti

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


// Header file for the classes RingTwinFloat and RingTwinFloatImpl


#include "CoCoA/ring.H"
#include "CoCoA/error.H"
//#include "CoCoA/FloatApprox.H"


namespace CoCoA
{

  class MantExp2; // fwd decl -- defined in FloatApprox.H

  class RingTwinFloatBase: public RingBase
  {
    RingElem mySymbolValue(const symbol& sym) const override;
  public: // member functions specific to RingTwinFloat implementations
    virtual long myPrecisionBits() const = 0;
    virtual void myRecvTwinFloat(RingElemRawPtr rawlhs, RingElemConstRawPtr rawx) const override = 0;
    virtual MantExp2 myExport(ConstRawPtr rawx) const = 0;
  };


  class RingTwinFloat: public ring
  {
  public:
    RingTwinFloat(const ring& R);
    explicit RingTwinFloat(const RingTwinFloatBase* RingPtr);
    // Default copy ctor & assignment work fine.
    // Default dtor works fine.
  public:
    const RingTwinFloatBase* operator->() const; // allow member fns to be called

  public:
    // Exception class to indicate when the error growth impedes further computation
    class InsufficientPrecision: public ErrorInfo
    {
    public:
      InsufficientPrecision();
    };
  };


  RingTwinFloat NewRingTwinFloat(const MachineInt& AccuracyBits);
  RingTwinFloat NewRingTwinFloat(const MachineInt& AccuracyBits, const MachineInt& BufferBits, const MachineInt& NoiseBits);
  const RingTwinFloatBase* RingTwinFloatPtr(const ring& R);
///  const RingTwinFloatBase* RingTwinFloatPtr(const ring& R, const char* const FnName);
  const RingTwinFloatBase* RingTwinFloatPtr(const ring& R, const ErrorContext& ErrCtx);
  bool IsRingTwinFloat(const ring& RR);
  long PrecisionBits(const RingTwinFloat& RR);
  bool IsPracticallyEqual(ConstRefRingElem x, ConstRefRingElem y); // applicable only if x & y are elems of a RingTwinFloat
  RingHom NewApproxHom(const ring& TwinFloat, const ring& R);
  void DebugPrint(std::ostream& out, ConstRefRingElem x); // print out both components of a Twin float

  MantExp2 MantissaAndExponent2(const RingElem& x); // only for elements of a RingTwinFloat!

  //---------------------------------------------------------------------------
  // Inline functions

  inline const RingTwinFloatBase* RingTwinFloat::operator->() const
  {
    return static_cast<const RingTwinFloatBase*>(ring::operator->());
  }


  inline bool IsRingTwinFloat(const ring& R)
  {
    return RingTwinFloatPtr(R) != nullptr;
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RingTwinFloat.H,v 1.20 2024/03/16 10:15:11 abbott Exp $
// $Log: RingTwinFloat.H,v $
// Revision 1.20  2024/03/16 10:15:11  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.19  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.18  2021/10/30 17:16:40  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.17  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.16  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.15  2018/05/18 14:39:25  bigatti
// -- fixed typo
//
// Revision 1.14  2018/05/18 14:35:59  bigatti
// -- removed include FloatApprox
//
// Revision 1.13  2016/04/07 14:35:45  abbott
// Summary: Removed code which disabled assignment of rings
//
// Revision 1.12  2016/03/25 20:21:33  abbott
// Summary: Added new fn MantissaAndExponent2 for elements of RingTwinFloat
//
// Revision 1.11  2014/07/30 14:28:02  abbott
// Summary: Added myExport for MantExp2
// Author: JAA
//
// Revision 1.10  2014/07/09 11:45:07  abbott
// Summary: Removed AsRingTwinFloat
// Author: JAA
//
// Revision 1.9  2011/11/09 13:50:01  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.8  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.7  2011/03/14 10:31:16  abbott
// Changed size_t into long (in fn interfaces).
//
// Revision 1.6  2010/10/08 15:48:33  abbott
// Commented out unused formal parameter (to avoid warning when compiling with pedantic).
//
// Revision 1.5  2010/10/01 15:20:33  bigatti
// -- added mySymbolValue
// -- added RingElem(R, sym)
//
// Revision 1.4  2010/01/19 17:37:27  abbott
// Changed arg type of pseudo-ctors NewRingTwinFloat into MachineInt
// (instead of size_t).
//
// Revision 1.3  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.2  2007/03/21 14:52:49  bigatti
// -- added AsRingFloat(R)
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2006/12/06 17:21:56  cocoa
// -- removed #include "config.H"
//
// Revision 1.6  2006/11/20 15:55:03  cocoa
// ring is now a class again.  Improved definitions of operator-> in derived classes.
//
// Revision 1.5  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.4  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.3  2006/10/06 10:12:07  cocoa
// Corrected position of #include.  Fixed incorrect #include.
//
// Revision 1.2  2006/07/17 19:14:53  cocoa
// First stage of big change.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.6  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.5  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
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
// Revision 1.11  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.10  2004/11/11 13:04:30  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.9  2004/11/05 15:34:33  cocoa
// Consequential change following from the renaming of
// FieldIdealImpl and the introduction of the new pseudo-ctor.
//
// Revision 1.8  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.7  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.6  2004/07/16 15:45:12  cocoa
// First stage of new RingElem implementation completed.
//
// Revision 1.5  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.4  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.3  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.2  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.1  2004/01/28 15:54:09  cocoa
// Sundry additions.
//

#endif
