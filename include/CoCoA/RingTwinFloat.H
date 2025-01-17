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

#endif
