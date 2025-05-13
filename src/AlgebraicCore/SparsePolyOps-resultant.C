//   Copyright (c)  2017,2022  John Abbott and Anna M. Bigatti

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


#include "CoCoA/SparsePolyOps-resultant.H"

#include "CoCoA/interrupt.H"
#include "CoCoA/ring.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/verbose.H"

//#include <vector>
using std::vector;
#include <iostream>
using std::endl;

namespace CoCoA
{

  namespace // anonymous
  {

    // Leading coeff viewed as poly in x
    RingElem lcx(ConstRefRingElem f, long x)
    {
      CoCoA_ASSERT(!IsZero(f));
      const vector<RingElem> c = CoeffVecWRT(f,indet(owner(f),x));
      return c.back();
    }

    // pseudo-remainder 
    RingElem prem(RingElem f, const RingElem& g, long x)
    {
      CoCoA_ASSERT(owner(f) == owner(g));
      CoCoA_ASSERT(!IsZero(g));
      if (IsZero(f)) return f;
      const ring& P = owner(f);
      const long degg = deg(g,x);
      if (degg == 0) return zero(P);
      long degf = deg(f,x);
      if (degf < degg) return f;
      const RingElem lcg = lcx(g,x);
      f *= power(lcg, degf-degg+1);
      while (degf >= degg)
      {
        const RingElem lcf = lcx(f,x);
        f -= (lcf/lcg)*g*power(indet(P,x), degf-degg); // exact division!
        if (IsZero(f)) return f;
        degf = deg(f,x);
      }
      return f;
    }

  } // end of namespace anonymous


  vector<RingElem> SubresultantSeq(ConstRefRingElem f, ConstRefRingElem g, long x)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (owner(g) != P)
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (x < 0 || x >= NumIndets(P))
      CoCoA_THROW_ERROR1(ERR::BadIndex);
    if (IsZero(f) || IsZero(g))
      CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    if (deg(f,x) < deg(g,x))  return SubresultantSeq(g, f, x);
    VerboseLog VERBOSE("SubresultantSeq");
    VERBOSE(90) << "Inputs:" << endl;
    VERBOSE(90) << "LPP(f) = " << LPP(f) << endl;
    VERBOSE(90) << "LPP(g) = " << LPP(g) << endl;
    vector<RingElem> S;
    RingElem s = power(lcx(g,x), deg(f,x)-deg(g,x));
    RingElem A = g;
    RingElem B = prem(f,-g,x);
    if (IsZero(B))
    { // g is a factor of f: two cases: g const or g not const
      if (!IsConstant(g))  { S.push_back(zero(P)); return S; }
      S.push_back(s*one(P));
      return S;
    }
    VERBOSE(90) << "Before loop LPP(prem) = " << LPP(B) << endl; // prem is sometimes huge...
    while (true)
    {
      CheckForInterrupt("Subresultant");
      if (IsZero(B))  return S;
      const long dA = deg(A,x);
      const long dB = deg(B,x);
      S.push_back(B);
      const long delta = dA-dB;
      RingElem C;
      if (delta > 1)
      {
        VERBOSE(90) << "case delta > 1; in fact delta = " << delta << endl;
        C = power(lcx(B,x), delta-1)*B/power(s, delta-1);
        VERBOSE(90) << "LPP(C) = " << LPP(C) << endl;
        S.push_back(C);
      }
      else
      {
        C = B;
      }
      if (dB == 0) return S;
      B = prem(A, -B, x);
      const RingElem D = power(s,delta)*lcx(A,x);
      B /= D; // exact div!
      if (!IsZero(B))  VERBOSE(90) << "LPP(B) = " << LPP(B) << endl;
      A = C;
      s = lcx(A,x);
    }
    // never get here
  }


  RingElem resultant(ConstRefRingElem f, ConstRefRingElem g) // f,g univariate; result is in CoeffRing!
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (owner(g) != P)
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (IsZero(f) || IsZero(g))
      CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    const long x = UnivariateIndetIndex(f);
    if (x < 0 || (!IsConstant(g) && x != UnivariateIndetIndex(g)))
      CoCoA_THROW_ERROR2(ERR::ReqUnivariate, "use 3 arg version for multivariate args");
    return ConstantCoeff(resultant(f,g,x));
  }

  RingElem resultant(ConstRefRingElem f, ConstRefRingElem g, long x)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (owner(g) != P)
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (IsZero(f) || IsZero(g))
      CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    if (x < 0 || x >= NumIndets(P))
      CoCoA_THROW_ERROR1(ERR::BadIndex);

    RingElem R = SubresultantSeq(f,g,x).back();
    if (IsZero(R) || deg(R,x) == 0)  return R; // redmine 1725... seems to be a hack, SubresultantSeq must be wrong
    return zero(P);
  }


  RingElem discriminant(ConstRefRingElem f)    ///< discr of univariate polynomial; result is in CoeffRing!
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsConstant(f))
      CoCoA_THROW_ERROR2(ERR::BadArg, "Poly must be non-constant");
    const long x = UnivariateIndetIndex(f);
    if (x < 0)
      CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    return ConstantCoeff(discriminant(f,x));
  }
  
  RingElem discriminant(ConstRefRingElem f, long x) ///< discr of multivariate polynomial
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsConstant(f))
      CoCoA_THROW_ERROR2(ERR::BadArg, "Poly must be non-constant");
    if (x < 0 || x >= NumIndets(P))
      CoCoA_THROW_ERROR1(ERR::BadIndex);

    const RingElem R = resultant(f,deriv(f,x),x);
    const long degf4 = deg(f)%4;
    const int sign = (degf4/2 == 0)?1:-1;
    return sign*(R/lcx(f,x));
  }

} // end of namespace CoCoA
  
