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

#include "CoCoA/FinGenModule.H"
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
    CoCoA_THROW_ERROR1(ERR::MixedModules);
    return *this; // just to keep compiler quiet
  }


  ConstRefRingElem ModuleElem::operator[](long pos) const
  {
    const module& M = owner(*this);
    if (!IsFinGenModule(M))  CoCoA_THROW_ERROR1(ERR::ReqFinGenModule);
    if (pos < 0 || pos >= NumCompts(M))
      CoCoA_THROW_ERROR1(ERR::BadIndex);
    return FinGenModulePtr(M)->myCompt(raw(*this), pos);
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
    CoCoA_THROW_ERROR1(ERR::MixedModules);
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
    CoCoA_THROW_ERROR1(ERR::MixedModules);
    return v; // just to keep compiler quiet
  }



  ModuleElem operator*(ConstRefRingElem r, const ModuleElem& v)
  {
    const ring& Rr = owner(r);
    const module& M = owner(v);
    const ring& R = RingOf(M);
    if (Rr != R)
    {
      if (RingID(Rr) > RingID(R))  CoCoA_THROW_ERROR1(ERR::MixedRings);
      const RingHom promote = AutomaticConversionHom(Rr,R,CoCoA_ERROR_CONTEXT); // throws ErrCannotConvert if auto-conv not possible
      return promote(r) * v;
    }
    // r is in correct ring
    ModuleElem ans(M);
    M->myMul(raw(ans), raw(r), raw(v));
    return ans;
  }


  ModuleElem operator*(const ModuleElem& v, ConstRefRingElem r)
  {
    const ring& Rr = owner(r);
    const module& M = owner(v);
    const ring& R = RingOf(M);
    if (Rr != R)
    {
      if (RingID(Rr) > RingID(R))  CoCoA_THROW_ERROR1(ERR::MixedRings);
      const RingHom promote = AutomaticConversionHom(Rr,R,CoCoA_ERROR_CONTEXT); // throws ErrCannotConvert if auto-conv not possible
      return v * promote(r);
    }

    // r is in the correct ring, just do the multoplication
    if (IsCommutative(R))  return r*v; // placeholder until below is implemented
    ModuleElem ans(M);
    CoCoA_THROW_ERROR1(ERR::NYI);  // BUG BUG
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
    CoCoA_THROW_ERROR1(ERR::MixedRings);
    return v; // just to keep compiler quiet
  }


  ModuleElem& operator+=(ModuleElem& v, const ModuleElem& w)
  {
    const module& Mv = owner(v);
    const module& Mw = owner(w);
    if (Mv != Mw)  CoCoA_THROW_ERROR1(ERR::MixedModules);

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
    CoCoA_THROW_ERROR1(ERR::MixedModules);
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
    CoCoA_THROW_ERROR1(ERR::MixedRings);
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
    CoCoA_THROW_ERROR1(ERR::MixedRings);
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
