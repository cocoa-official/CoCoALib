//   Copyright (c)  2023  John Abbott and Anna M. Bigatti
//   Original author Nico Mexis (2023)

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


#include "CoCoA/SparsePolyOps-IsPrimitivePoly.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"

#include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {

    // This fn should really be made public
    RingElem PowerMod(const RingElem& f, const BigInt& N, const RingElem& modulus)
    {
      // assume in polynomial ring
      // assume f = NR(f, {modulus})
      // assume N > 0
      if (IsOne(N)) return f;
      CheckForInterrupt("PowerMod");
      const RingElem EvenPwr = PowerMod(NR(f*f,{modulus}), N/2, modulus);
      if (IsEven(N)) return EvenPwr;
      return NR(f*EvenPwr, {modulus});
    }
    
  } // end of namespace anonymous
  

  bool IsPrimitivePoly_NoArgChecks(ConstRefRingElem f)
  {
    // Assume coeffs are in finite field, & f is univariate
    if (IsZero(ConstantCoeff(f))) return false;
    if (!IsOne(LC(f))) return false;
    if (!IsIrred(f)) return false;
    if (deg(f) == 1) return true; // or error???

    const ring &Px = owner(f);
    const long IndetIndex = UnivariateIndetIndex(f);
    const RingElem& x = indet(Px, IndetIndex);

    // We check if f is an n-th primitive polynomial in ZZ/(p)
    const long n = deg(f);
    const BigInt p = characteristic(CoeffRing(Px));
    const BigInt M = power(p, n) - 1;
    const vector<BigInt> fac = factor(M).myFactors(); // primes dividing M
    // Should check that the factors found are (prob)prime  -- all facs should be prime for M < 2^80 (I guess)
//?? SLOWER???   sort(fac.begin(), fac.end(), [](const BigInt& a, const BigInt& b) { return (a>b); } );  // reverse sorted!
    const auto CheckExponent = [&x, &M, &f](const BigInt& m) { return IsOne(PowerMod(x, M/m, f)); };
    return none_of(fac.begin(), fac.end(), CheckExponent);
//    return none_of(fac.begin(), fac.end(), [f, Px, IndetIndex, M](const BigInt &m) {
//      return IsOne(NR(IndetPower(Px, IndetIndex, M / m), {f}));
//    });
  }

  bool IsPrimitivePoly(ConstRefRingElem f)
  {
    if (IsZero(f)) return false;
    const ring &Px = owner(f);
    if (!IsSparsePolyRing(Px))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (UnivariateIndetIndex(f) < 0)  CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    const ring &P = CoeffRing(Px);
    if (!IsFiniteField(P) || !IsZZ(BaseRing(P)) /*|| !IsQuotientRing(P)*/) // TODO: Does that only allow ZZ/(p)?
      CoCoA_THROW_ERROR2(ERR::BadRing, "must be finite field or based on ZZ");

    return IsPrimitivePoly_NoArgChecks(f);
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-IsPrimitivePoly.C,v 1.4 2024/08/01 16:53:02 bigatti Exp $
// $Log: SparsePolyOps-IsPrimitivePoly.C,v $
// Revision 1.4  2024/08/01 16:53:02  bigatti
// Summary: error codes: Req(Elem)SparsePolyRing;
// now using CoCoA_THROW_ERROR1/2
//
// Revision 1.3  2024/07/02 09:57:58  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.2  2023/03/17 19:26:49  abbott
// Summary: Added comment about sorting the factors
//
// Revision 1.1  2023/03/15 21:33:24  abbott
// Summary: New fn IsPrimitivePoly
//
//
//
