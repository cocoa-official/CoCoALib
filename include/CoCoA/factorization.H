#ifndef CoCoA_factorization_H
#define CoCoA_factorization_H

//   Copyright (c)  2009,2014  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/BigInt.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"
// clang says the following include directive is needed...
#include "CoCoA/VectorOps.H"

#include <vector>
// using std::vector;
#include <iostream>
//using std::ostream;  -- needed for the templated operator<<

namespace CoCoA
{

  // This is really just a struct: it merely contains the result of a factorization.
  template <typename T>
  class factorization
  {
  public:
    explicit factorization(const T& RemainingFactor):
        myFactorVec(),
        myMultiplicityVec(),
        myRemainingFactorValue(RemainingFactor)
      {
        ourCheckNotField(RemainingFactor);
        ourCheckNotZeroDiv(RemainingFactor);
      }

    factorization(const std::vector<T>& factors, const std::vector<long>& multiplicities, const T& RemainingFactor):
        myFactorVec(factors),
        myMultiplicityVec(multiplicities),
        myRemainingFactorValue(RemainingFactor)
      {
        ourCheckNotField(RemainingFactor);
        ourCheckNotZeroDiv(RemainingFactor);
        // Check vecs have same len; all mults are positive; all factors not invertible nor zero-divs:
        const long n = len(factors);
        if (n != len(multiplicities))
          CoCoA_THROW_ERROR1(ERR::IncompatDims);
        for (long i=0; i < n; ++i)
        {
          if (multiplicities[i] <= 0)
            CoCoA_THROW_ERROR2(ERR::ReqPositive, "multiplicities");
          ourCheckCompatible(factors[i], RemainingFactor);
          ourCheckNotZeroDiv(factors[i]);
          ourCheckNotUnit(factors[i]);
        }
      }
    // default copy ctor is OK
    // default assignment is OK
    // default dtor is OK

    const std::vector<T>& myFactors() const { return myFactorVec; }
    const std::vector<long>& myMultiplicities() const { return myMultiplicityVec; }
    const T& myRemainingFactor() const { return myRemainingFactorValue; }

    void myAppend(const T& fac, long mult);
    void myNewRemainingFactor(const T& RemFac);

    friend void swap(factorization<T>& A, factorization<T>& B) /*?noexcept?*/ { swap(A.myRemainingFactorValue, B.myRemainingFactorValue); std::swap(A.myFactorVec, B.myFactorVec); std::swap(A.myMultiplicityVec, B.myMultiplicityVec);}

  private: // Data members
    std::vector<T> myFactorVec;
    std::vector<long> myMultiplicityVec;
    T myRemainingFactorValue;
  private:
    static void ourCheckNotField(const T& fac);   // Each new type must impl these fns (throw if check fails, o/w do nothing)
    static void ourCheckNotZeroDiv(const T& fac); // Each new type must impl these fns (throw if check fails, o/w do nothing)
    static void ourCheckNotUnit(const T& fac);    // Each new type must impl these fns (throw if check fails, o/w do nothing)
    static void ourCheckCompatible(const T& fac, const T& RemFac); // give error if fac*RemFac cannot be computed
    static void ourConsistencyCheck(const std::vector<T>& factors, const std::vector<long>& multiplicities, const T& RemainingFactor);
  };


  template <typename T>
  std::ostream& operator<<(std::ostream& out, const factorization<T>& facs)
  {
    return out << "factorization(myFactors=" << facs.myFactors()
               << ", myMultiplicities=" << facs.myMultiplicities()
               << ", myRemainingFactor=" << facs.myRemainingFactor() << ")";
  }


  // Deliberately no default impls of ourCheckNonZero and ourCheckNotUnit!

  // default impl of ourCheckCompatible  (does nothing)
  template <typename T>
  void factorization<T>::ourCheckCompatible(const T& /*fac*/, const T& /*RemFac*/)
  {}


  // Specialization of ourCheckNotField for RingElem
  template <>
  inline void factorization<RingElem>::ourCheckNotField(const RingElem& fac)
  {
    if (IsField(owner(fac)))
      CoCoA_THROW_ERROR2(ERR::BadRing, "must be non-field");
  }


  // Specialization of ourCheckNotZeroDiv for RingElem
  template <>
  inline void factorization<RingElem>::ourCheckNotZeroDiv(const RingElem& fac)
  {
    if (IsZeroDivisor(fac))
      CoCoA_THROW_ERROR2(ERR::BadArg, "must be non-zero-divisor");
  }


  // Specialization ourCheckNotUnit for RingElem
  template <>
  inline void factorization<RingElem>::ourCheckNotUnit(const RingElem& fac)
  {
    if (IsInvertible(fac))
      CoCoA_THROW_ERROR2(ERR::BadArg, "must be non-unit");
  }


  // Specialization of ourCheckCompatible for RingElem
  template <>
  inline void factorization<RingElem>::ourCheckCompatible(const RingElem& fac, const RingElem& RemFac)
  { if (owner(fac) != owner(RemFac))  CoCoA_THROW_ERROR1(ERR::MixedRings); }


  //-------------------------------------------------------
  // Specialization of ourCheckNotField for BigInt
  template <>
  inline void factorization<BigInt>::ourCheckNotField(const BigInt& /*fac*/)
  { /*nothing to do*/ }


  // Specialization of ourCheckNotZeroDiv for BigInt
  template <>
  inline void factorization<BigInt>::ourCheckNotZeroDiv(const BigInt& fac)
  { if (IsZero(fac))  CoCoA_THROW_ERROR1(ERR::ReqNonZero); }


  // Specialization of ourCheckNotUnit for BigInt
  template <>
  inline void factorization<BigInt>::ourCheckNotUnit(const BigInt& fac)
  {
    if (IsOne(fac) || IsMinusOne(fac))
      CoCoA_THROW_ERROR2(ERR::BadArg, "must be non-unit");
  }


  //-------------------------------------------------------
  // Specialization of ourCheckNotField for long
  template <>
  inline void factorization<long>::ourCheckNotField(const long& /*fac*/)
  { /*nohting to do*/ }


  // Specialization of ourCheckNotZeroDiv for long
  template <>
  inline void factorization<long>::ourCheckNotZeroDiv(const long& fac)
  { if (fac == 0)  CoCoA_THROW_ERROR1(ERR::ReqNonZero); }


  // Specialization of ourCheckNotUnit for long
  template <>
  inline void factorization<long>::ourCheckNotUnit(const long& fac)
  {
    if (fac == 1 || fac == -1)
      CoCoA_THROW_ERROR2(ERR::BadArg, "must be non-unit");
  }

  // ***TEMPLATE SPECIALIZATIONS*** of ourCheckNotZeroDiv & ourCheckNotUnit also in DUPFp.H


  template <typename T>
  void factorization<T>::myAppend(const T& fac, long mult)
  {
    if (mult <= 0)
      CoCoA_THROW_ERROR2(ERR::ReqPositive, "multiplicity");
    ourCheckCompatible(fac, myRemainingFactorValue);
    ourCheckNotZeroDiv(fac);
    ourCheckNotUnit(fac);
    myMultiplicityVec.reserve(1+myMultiplicityVec.size()); // for exception safety
    myFactorVec.push_back(fac);
    myMultiplicityVec.push_back(mult);
  }


  template <typename T>
  void factorization<T>::myNewRemainingFactor(const T& RemFac)
  {
    ourCheckCompatible(RemFac, myRemainingFactorValue);
    ourCheckNotZeroDiv(RemFac);
    myRemainingFactorValue = RemFac;
  }


} // end of namespace CoCoA

#endif
