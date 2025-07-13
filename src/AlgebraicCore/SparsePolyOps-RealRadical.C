 //   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/SparsePolyOps-RealRadical.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/SparsePolyOps-SturmSeq.H"
#include "CoCoA/factorization.H"
#include "CoCoA/factor.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/utils.H"

#include <algorithm>
using std::max;

#include <vector>
using std::vector;


namespace CoCoA
{

  // export ????
  vector<long> degs(const RingElem& f)
  {
    const ring& P = owner(f);
    const int nvars = NumIndets(P);
    vector<long> ans(nvars);
    vector<long> expv(nvars);
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      exponents(expv, PP(it));
      for (int i=0; i < nvars; ++i)
        ans[i] = max(ans[i], expv[i]);
    }
    return ans;
  }


  namespace // anonymous
  {
    // Quick(ish): may return uncertain3
    bool3 HasRealRoot3_univariate(const RingElem& f)
    {
      const PolyRing P(owner(f), CoCoA_ERROR_CONTEXT);
////      CoCoA_ASSERT(IsPolyRing(P));
      const int var = UnivariateIndetIndex(f);
      CoCoA_ASSERT(var >= 0); // poly is univariate, ring may have several indets!
      if (IsZero(f)) return true3;
      if (deg(f) == 0) return false3;  // non-zero const
      if (IsOdd(deg(f))) return true3; // odd deg
      if (IsMonomial(f)) return true3; // non-const monomial

      const RingHom coeff = CoeffEmbeddingHom(P);
      RingElem g = f;
      const RingElem& x = indet(P,var);
      const int s = sign(LC(g));
      while (true)
      {
        const long d = deg(g);
        if (IsOdd(d)) return uncertain3;
        const RingElem lcg = coeff(LC(g));
        if (sign(LC(g)) != s) return uncertain3;
        if (d == 0) return false3; // it was sum of squares (with non-zero const part)
        RingElem term = power(x,d/2);
        RingElem rem = g - lcg*power(term,2);
        if (IsZero(rem)) return true3;  // sum of sqs (with zero const part)
        const long d2 = deg(rem);
        if (IsEven(deg(rem)) && sign(LC(rem)) == s) { g = rem; continue; }
        if (d2 < d/2) return uncertain3;
        RingElem term2 = coeff(LC(rem))/(2*lcg) * power(x, d2 - d/2);
        g -= lcg*power(term+term2, 2);
      }
      // never get here
    }

  
    bool3 HasRealRoot3_multivariate(const RingElem& f)
    {
      const ring& P = owner(f);
      CoCoA_ASSERT(IsPolyRing(P));
      const int nvars = NumIndets(P);
      CoCoA_ASSERT(nvars > 1);

      if (IsZero(f)) return true3;       // zero
      if (IsOne(LPP(f))) return false3;  // non-zero const
      if (IsMonomial(f)) return true3;   // just a non-const monomial
      
      // if deg is odd wrt some variable then surely there is a real root.
      const vector<long> d = degs(f);
      for (int i = 0; i < nvars; ++i)
        if (IsOdd(d[i])) return true3;

      // Poly is multivariate, has >= 2 terms, even deg in all indets.
      // Check if repr as sum of terms is explicit repr as sum of sqs.
      vector<long> expv(nvars);
      const int s = sign(LC(f));
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
      {
        if (IsOne(PP(it)))
          return bool3(sign(coeff(it)) != s);  // sum of sqs with non-zero const term
        if (sign(coeff(it)) != s) return uncertain3;
        exponents(expv, PP(it));
        for (int i=0; i < nvars; ++i)
          if (IsOdd(expv[i])) return uncertain3;
      }
      return true3; // because setting all vars to 0 makes f vanish
    }

  } // end of namespace anonymous

  
  // Return value:
  //   false3      then f has no real roots
  //   true3       then f surely has real roots
  //   uncertain3  then we don't know  (only for multivariate inputs)
  bool3 HasRealRoot3(const RingElem& f)
  {
    const int var = UnivariateIndetIndex(f);
    if (var < 0)
      return HasRealRoot3_multivariate(f); // cop out for multivariate
    const bool3 QuickCheck = HasRealRoot3_univariate(f);
    if (!IsUncertain3(QuickCheck))
      return QuickCheck;
    return bool3(NumRealRoots(f) > 0);
  }


  RingElem RealRadical(const RingElem& f)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P) || !IsQQ(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::BadRing);
    const factorization<RingElem> FacInfo = factor(radical(f));
    const vector<RingElem>& fac = FacInfo.myFactors();
    RingElem ans = one(P);
    for (int i=0; i < len(fac); ++i)
    {
      if (!IsFalse3(HasRealRoot3(fac[i])))
        ans *= fac[i];
    }
    return ans;
  }


} // end of namespace CoCoA
