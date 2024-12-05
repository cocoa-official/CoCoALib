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

#endif
