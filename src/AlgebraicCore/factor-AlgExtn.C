//   Copyright (c)  2018  John Abbott and Anna M. Bigatti
//   Authors:  2018  Anna M. Bigatti
//             2017  Elisa Palezzato (original code and tests in CoCoA-5 pkg)

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


// Source code for ideals in SparsePolyRing (functions and member functions)

#include "CoCoA/factor.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyOps-ideal.H" // for GBasis
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/factorization.H"
#include "CoCoA/ideal.H"
#include "CoCoA/random.H" // for RandomLongStream
#include "CoCoA/time.H"
#include "CoCoA/verbose.H"

#include <algorithm>
#include <iostream>
using std::endl;
#include <list>
using std::vector;

namespace CoCoA
{

  namespace{ // anonymous


    template <typename T>
      std::vector<T> first(const std::vector<T>& v, long n)
    {
      std::vector<T> w;
      if (n < 0 || n>len(v))  CoCoA_THROW_ERROR1(ERR::OutOfRange);
      for (long i=0; i<n; ++i)
        w.push_back(v[i]);
      return w;
    }

    template <typename T>
      std::vector<T> last(const std::vector<T>& v, long n)
    {
      std::vector<T> w;
      if (n < 0 || n>len(v))  CoCoA_THROW_ERROR1(ERR::OutOfRange);
      for (long i=len(v)-n; i<len(v); ++i)
        w.push_back(v[i]);
      return w;
    }

    template <typename T>
      std::vector<T> concat(const std::vector<T>& v1, const std::vector<T>& v2)
    {
      std::vector<T> w;
      for (long i=0; i<len(v1); ++i) w.push_back(v1[i]);
      for (long i=0; i<len(v2); ++i) w.push_back(v2[i]);
      return w;
    }

    RingHom Hom_P_Px(const ring& P)
    {
      /// use anonymous symbols!
      const vector<symbol> NewIndets = NewSymbols(NumIndets(P)+1);
      SparsePolyRing Px = NewPolyRing(CoeffRing(P), NewIndets);
      return PolyAlgebraHom(P, Px, first(indets(Px), NumIndets(P)));
    }


    RingElem LC0(ConstRefRingElem m)
    {
      if (IsZero(m)) return zero(CoeffRing(owner(m)));
      return LC(m);
    }
    

  }  // end of anonymous namespace



  factorization<RingElem> factor_AlgExtn(ConstRefRingElem f)
  {
    VerboseLog VERBOSE("factor_AlgExtn");
    const ring& Lx = owner(f);  // L is Pa/I
    if (!IsQuotientRing(CoeffRing(Lx)))
      CoCoA_THROW_ERROR1(ERR::NotQuotientRing);
    long idx = UnivariateIndetIndex(f);
    if (idx < 0)
      CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    const ring& Pa = BaseRing(CoeffRing(Lx));
    const RingHom phi = Hom_P_Px(Pa);
    const std::vector<RingElem> coeffs = CoeffVecWRT(f, indet(Lx, idx));
    const ring Pax = codomain(phi);
    const RingElem X = indets(Pax).back();
    RingElem fp(Pax);
    for (long i=0; i<len(coeffs); ++i)
      fp += phi(CanonicalRepr(LC0(coeffs[i])))*power(X,i);
    const vector<RingElem> mP = phi(gens(DefiningIdeal(CoeffRing(Lx))));
    ideal J = ideal(mP) + ideal(fp);
    double T = CpuTime();
    const vector<ideal> PD = PrimaryDecomposition(J);
    VERBOSE(40) << "num factors " << len(PD) << endl;
    VERBOSE(90) << "PrimaryDec time = " << CpuTime()-T << endl;
    vector<RingElem> classX = ChainCanonicalHom(Pa,Lx)( indets(Pa) );
    classX.push_back(indet(Lx,idx));
    RingHom psi = PolyRingHom(Pax,Lx,ChainCanonicalHom(CoeffRing(Pax),Lx),classX);
    vector<RingElem> components;
    vector<RingElem> GB;
    for (long i=0; i<len(PD); ++i) // ring is a PID
    {
      if (HasGBasis(PD[i]))
        GB = GBasis(ideal(psi(GBasis(PD[i]))));
      else
        GB = GBasis(ideal(psi(gens(PD[i]))));
      CoCoA_ASSERT(len(GB) == 1);
      components.push_back(GB[0]);
    }
    VERBOSE(99) << components << endl;
    vector<RingElem> facs;
    vector<long> mult;
    RingElem rem=one(Lx);
    for (long i=0; i<len(components); ++i)
    {
      const factorization<RingElem> F = SqFreeFactor(components[i]);
      CoCoA_ASSERT(len(F.myFactors())==1);
      facs.push_back(F.myFactors()[0]);
      mult.push_back(F.myMultiplicities()[0]);
      rem *= F.myRemainingFactor();
    }
    VERBOSE(99) << facs << endl;
    VERBOSE(99) << mult << endl;
    return factorization<RingElem>(facs, mult, LC(f)*rem);
  }


} // end of namespace CoCoA
