//   Copyright (c)  2005,2008  John Abbott and Anna M. Bigatti

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


// Source code for abstract class module and friends

#include "CoCoA/module.H"

#include "CoCoA/FGModule.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/ring-AutomaticConversion.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/submodule.H"

//#include <iostream>
using std::ostream;

namespace CoCoA
{


  // C++ needs this function to be defined
  ModuleBase::~ModuleBase()
  {}

  bool module::operator==(const module& M) const
  {
    if (IsFreeModule(*this) && IsFreeModule(M))
      return mySmartPtr==M.mySmartPtr;
    return IsContained(M,*this) && IsContained(*this,M);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Operations on ModuleElems

  ModuleElem::ModuleElem(const module& M):
      myM(M)
  {
    myM->myNew(myValue);
  }


  ModuleElem::ModuleElem(const ModuleElem& copy):
      myM(copy.myM)
  {
    myM->myNew(myValue, raw(copy));
  }


  ModuleElem::~ModuleElem()
  {
    myM->myDelete(myValue);
  }


  ModuleElem& ModuleElem::operator=(const ModuleElem& rhs)
  {
    if (this == &rhs) return *this;
    const module& Mlhs = owner(*this);
    const module& Mrhs = owner(rhs);
    if (Mlhs == Mrhs)
    {
      Mlhs->myAssign(raw(*this), raw(rhs));
      return *this;
    }
    CoCoA_THROW_ERROR(ERR::MixedModules, "Mvector = Mvector");
    return *this; // just to keep compiler quiet
  }


  ConstRefRingElem ModuleElem::operator[](long pos) const
  {
    const module& M = owner(*this);
    if (!IsFGModule(M))  CoCoA_THROW_ERROR(ERR::ReqFGModule, "ModuleElem[pos]");
    if (pos < 0 || pos >= NumCompts(M))
      CoCoA_THROW_ERROR(ERR::BadIndex, "ModuleElem[pos]");
    return FGModulePtr(M)->myCompt(raw(*this), pos);
  }


  ModuleElem operator-(const ModuleElem& v)
  {
    const module& M = owner(v);
    ModuleElem ans(M);
    M->myNegate(raw(ans), raw(v));
    return ans;
  }


  ModuleElem operator+(const ModuleElem& v, const ModuleElem& w)
  {
    const module& Mv = owner(v);
    const module& Mw = owner(w);
    if (Mv == Mw)
    {
      ModuleElem ans(Mv);
      Mv->myAdd(raw(ans), raw(v), raw(w));
      return ans;
    }
    CoCoA_THROW_ERROR(ERR::MixedModules, "Mvector + Mvector");
    return v; // just to keep compiler quiet
  }


  ModuleElem operator-(const ModuleElem& v, const ModuleElem& w)
  {
    const module& Mv = owner(v);
    const module& Mw = owner(w);
    if (Mv == Mw)
    {
      ModuleElem ans(Mv);
      Mv->mySub(raw(ans), raw(v), raw(w));
      return ans;
    }
    CoCoA_THROW_ERROR(ERR::MixedModules, "Mvector - Mvector");
    return v; // just to keep compiler quiet
  }



  ModuleElem operator*(ConstRefRingElem r, const ModuleElem& v)
  {
    CoCoA_STATIC_ERROR_MESG(ErrCannotConvert, ERR::MixedRings,"RingElem*ModuleElem");
    const ring& Rr = owner(r);
    const module& M = owner(v);
    const ring& R = RingOf(M);
    if (Rr != R)
    {
      if (RingID(Rr) > RingID(R)) ThrowException(ErrCannotConvert);
      const RingHom promote = AutomaticConversionHom(Rr,R,ErrCannotConvert); // throws ErrCannotConvert if auto-conv not possible
      return promote(r) * v;
    }
    // r is in correct ring
    ModuleElem ans(M);
    M->myMul(raw(ans), raw(r), raw(v));
    return ans;
  }


  ModuleElem operator*(const ModuleElem& v, ConstRefRingElem r)
  {
    CoCoA_STATIC_ERROR_MESG(ErrCannotConvert, ERR::MixedRings,"ModuleElem*RingElem");
    const ring& Rr = owner(r);
    const module& M = owner(v);
    const ring& R = RingOf(M);
    if (Rr != R)
    {
      if (RingID(Rr) > RingID(R)) ThrowException(ErrCannotConvert);
      const RingHom promote = AutomaticConversionHom(Rr,R,ErrCannotConvert); // throws ErrCannotConvert if auto-conv not possible
      return v * promote(r);
    }

    // r is in the correct ring, just do the multoplication
    if (IsCommutative(R)) return r*v; // placeholder until below is implemented
    ModuleElem ans(M);
    CoCoA_THROW_ERROR(ERR::NYI,"v*r");  // BUG BUG
    M->myMul(raw(ans), raw(r), raw(v));  // BUG should mult on right!!!
    return ans;
  }


  ModuleElem operator/(const ModuleElem& v, ConstRefRingElem r)
  {
    const ring& R = owner(r);
    const module& M = owner(v);
    if (R == RingOf(M))
    {
      ModuleElem ans(M);
      M->myDiv(raw(ans), raw(r), raw(v));  /// BUGLY;  mult on L or on R ????
      return ans;
    }
    CoCoA_THROW_ERROR(ERR::MixedRings, "Mvector / RingElem");
    return v; // just to keep compiler quiet
  }


  ModuleElem& operator+=(ModuleElem& v, const ModuleElem& w)
  {
    const module& Mv = owner(v);
    const module& Mw = owner(w);
    if (Mv != Mw) CoCoA_THROW_ERROR(ERR::MixedModules, "Mvector += Mvector");

    Mv->myAdd(raw(v), raw(v), raw(w));
    return v;
  }


  ModuleElem& operator-=(ModuleElem& v, const ModuleElem& w)
  {
    const module& Mv = owner(v);
    const module& Mw = owner(w);
    if (Mv == Mw)
    {
      Mv->mySub(raw(v), raw(v), raw(w));
      return v;
    }
    CoCoA_THROW_ERROR(ERR::MixedModules, "Mvector -= Mvector");
    return v; // just to keep compiler quiet
  }



  ModuleElem& operator*=(ModuleElem& v, ConstRefRingElem r)
  {
    const ring& R = owner(r);
    const module& M = owner(v);
    if (R == RingOf(M))
    {
      M->myMul(raw(v), raw(r), raw(v));
      return v;
    }
    CoCoA_THROW_ERROR(ERR::MixedRings, "Mvector *= RingElem");
    return v; // just to keep compiler quiet
  }


  ModuleElem& operator/=(ModuleElem& v, ConstRefRingElem r)
  {
    const ring& R = owner(r);
    const module& M = owner(v);
    if (R == RingOf(M))
    {
      M->myDiv(raw(v), raw(r), raw(v));
      return v;
    }
    CoCoA_THROW_ERROR(ERR::MixedRings, "Mvector /= RingElem");
    return v; // just to keep compiler quiet
  }


  ModuleElem operator*(const MachineInt& n, const ModuleElem& v)
  {
    return RingElem(RingOf(owner(v)), n)*v;
  }


  ModuleElem operator*(const ModuleElem& v, const MachineInt& n)
  {
    return RingElem(RingOf(owner(v)), n)*v;
  }


  ModuleElem operator/(const ModuleElem& v, const MachineInt& n)
  {
    return v/RingElem(RingOf(owner(v)), n);
  }


  ModuleElem& operator*=(ModuleElem& v, const MachineInt& n)
  {
    return v *= RingElem(RingOf(owner(v)), n);
  }


  ModuleElem& operator/=(ModuleElem& v, const MachineInt& n)
  {
    return v /= RingElem(RingOf(owner(v)), n);
  }


  // Arith between ModuleElems and BigInts
  ModuleElem operator*(const BigInt& N, const ModuleElem& v)
  {
    return RingElem(RingOf(owner(v)), N)*v;
  }

  ModuleElem operator*(const ModuleElem& v, const BigInt& N)
  {
    return RingElem(RingOf(owner(v)), N)*v;
  }

  ModuleElem operator/(const ModuleElem& v, const BigInt& N)
  {
    return v/RingElem(RingOf(owner(v)), N);
  }


  ModuleElem& operator*=(ModuleElem& v, const BigInt& N)
  {
    return v *= RingElem(RingOf(owner(v)), N);
  }

  ModuleElem& operator/=(ModuleElem& v, const BigInt& N)
  {
    return v /= RingElem(RingOf(owner(v)), N);
  }

  // Arith between ModuleElems and BigRats
  ModuleElem operator*(const BigRat& q, const ModuleElem& v)
  {
    return RingElem(RingOf(owner(v)), q)*v;
  }

  ModuleElem operator*(const ModuleElem& v, const BigRat& q)
  {
    return RingElem(RingOf(owner(v)), q)*v;
  }

  ModuleElem operator/(const ModuleElem& v, const BigRat& q)
  {
    return v/RingElem(RingOf(owner(v)), q);
  }


  ModuleElem& operator*=(ModuleElem& v, const BigRat& q)
  {
    return v *= RingElem(RingOf(owner(v)), q);
  }

  ModuleElem& operator/=(ModuleElem& v, const BigRat& q)
  {
    return v /= RingElem(RingOf(owner(v)), q);
  }



  std::ostream& operator<<(std::ostream& out, const ModuleElem& v)
  {
    if (!out) return out;  // short-cut for bad ostreams
    owner(v)->myOutput(out, raw(v));
    return out;
  }


  std::ostream& operator<<(std::ostream& out, const module& M)
  {
    if (!out) return out;  // short-cut for bad ostreams
    M->myOutputSelf(out);
    return out;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const ModuleElem& v)
  {
    owner(v)->myOutput_OM(OMOut, raw(v));
    return OMOut;
  }

  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const module& M)
  {
    M->myOutputSelf_OM(OMOut);
    return OMOut;
  }


  bool IsZero(const ModuleElem& v)
  {
    return owner(v)->myIsZero(raw(v));
  }


  bool operator==(const ModuleElem& x, const ModuleElem& y)
  {
    return owner(x)->myIsEqual(raw(x), raw(y));
  }


  bool operator!=(const ModuleElem& x, const ModuleElem& y)
  {
    //    return !owner(x)->myIsEqual(raw(x), raw(y));
    return !(x==y);
  }



} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/module.C,v 1.23 2024/07/03 14:46:49 bigatti Exp $
// $Log: module.C,v $
// Revision 1.23  2024/07/03 14:46:49  bigatti
// Summary: error codes: Not... changed into ReqFGModule, ReqFreeModule, ReqModuleSpPR
//
// Revision 1.22  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.21  2022/02/08 20:18:55  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.20  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.19  2020/06/22 15:44:33  abbott
// Summary: Some fns now do automatic ring promotion
//
// Revision 1.18  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.17  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.16  2014/07/30 14:13:24  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.15  2014/07/09 14:27:53  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.14  2013/06/06 05:54:41  bigatti
// -- resorted includes
//
// Revision 1.13  2013/06/03 15:32:27  bigatti
// -- better fix for operator==
//
// Revision 1.12  2013/06/03 14:02:33  bigatti
// -- changed: operator== in now defined in .C
//
// Revision 1.11  2013/01/23 14:09:48  bigatti
// -- removed "inline" for "operator=="
// -- redefined "operator!="
//
// Revision 1.10  2012/10/11 14:37:02  abbott
// Removed non-const operator[], thereby removing the need for RefRingElem.
//
// Revision 1.9  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.8  2011/08/24 10:32:04  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.7  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.6  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.5  2009/12/03 17:26:34  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.4  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
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
// Revision 1.4  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.3  2006/11/02 13:25:43  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/27 19:09:45  cocoa
// Replaced some member functions of CoCoA::symbol by friend functions.
// Removed some include dependency on symbol.H
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.6  2006/04/21 14:56:33  cocoa
// Changed return type of myCompt member function: now it returns a
// ConstRefRingElem instead of a RingElem (i.e. a copy).
//
// Revision 1.5  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.4  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.3  2005/11/29 13:04:47  cocoa
// -- added "const" to myCompt argument
//
// Revision 1.2  2005/11/24 16:09:38  cocoa
// -- added operator[] for ModuleElem
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.4  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.3  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.6  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.5  2004/11/11 13:56:25  cocoa
// -- moved CVS log to the bottom
//
// Revision 1.4  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.3  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.2  2004/01/28 15:35:37  cocoa
// Major update -- was very much "old style" code, and didn't
// compiler under the new organization.
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.1  2003/05/30 15:10:54  abbott
// Initial revision
//
//

