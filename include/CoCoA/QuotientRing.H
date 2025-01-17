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

#endif
