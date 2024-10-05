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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/factorization.H,v 1.21 2024/07/31 07:52:47 abbott Exp $
// $Log: factorization.H,v $
// Revision 1.21  2024/07/31 07:52:47  abbott
// Summary: Switched to new CoCoA_THROW_ERROR macros (with suffix 1 or 2)
//
// Revision 1.20  2022/06/08 15:15:41  abbott
// Summary: Added "explicit"; added new fn swap
//
// Revision 1.19  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.18  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.17  2020/10/27 09:47:38  abbott
// Summary: Added check that factorization is not over a field
//
// Revision 1.16  2020/10/14 20:04:00  abbott
// Summary: Some tidying; improved error mesg (redmine 1499)
//
// Revision 1.15  2020/10/05 19:30:06  abbott
// Summary: Changed err mesg if multiplicity not positive (redmine 1499)
//
// Revision 1.14  2020/06/17 15:49:21  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.13  2018/05/18 12:13:37  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.12  2018/05/17 16:06:04  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.11  2014/07/31 14:45:16  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.10  2014/04/02 14:41:01  abbott
// Summary: Added missing arg check in ctor; removed 2 lines of commented out cruft
// Author: JAA
//
// Revision 1.9  2014/03/25 15:46:25  abbott
// Summary: Updated copyright date
// Author: JAA
//
// Revision 1.8  2014/03/25 11:36:22  abbott
// Summary: Better design: replaced template impls by template specializations for RingElem
// Author: JAA
//
// Revision 1.7  2014/03/24 12:09:20  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.6  2013/10/22 14:00:48  abbott
// Added ourConsistencyCheck.
//
// Revision 1.5  2013/04/11 15:06:25  abbott
// Added new include (because clang needs it).
//
// Revision 1.4  2012/10/05 09:30:00  abbott
// Changed myExponents into myMultiplicities.
//
// Revision 1.3  2012/01/30 11:05:19  abbott
// Changed include "CoCoA/io.H" into include <iostream>
//
// Revision 1.2  2011/12/05 16:55:43  bigatti
// -- changed: MachineInteger --> MachineInt (just in comment)
//
// Revision 1.1  2009/06/05 12:14:55  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from ZZ.H/C
//   removed gcd from MachineInt.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
//

#endif
