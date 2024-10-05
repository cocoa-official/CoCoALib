//   Copyright (c)  2006,2009  John Abbott and Anna M. Bigatti

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


// These are from directory TmpFactorDir/
#include "mpz_alias.h"
#include "primes.h"
#include "FF.h"
#include "DUPFF.h"
#include "DUPFFfactor.h"
//#include "DUPFFlist.h"
#include "WGD.h"

// These are from directory TmpFactorDir/multivariate/
///#include "DMPFF.h"
///#include "DMPFFgcd.h"
#include "DMPZgcd.h"
#include "DMPZfactor.h"
#include "DMPZfactor_modp.h"

/***************************************************************************/
/* CoCoALib includes */

#include "CoCoA/factor.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/convert.H"
#include "CoCoA/utils.H"

//#include <iostream>
#include <cstdlib>
using std::malloc;// someone actually uses malloc!!!!
#include <vector>
using std::vector;

/***************************************************************************/

namespace CoCoA
{


  DUPFF ConvertToDUPFF(ConstRefRingElem f)
  {
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing Rx = owner(f);
    const ring R = CoeffRing(Rx);
    CoCoA_ASSERT(IsField(R) && IsQuotientRing(R) && IsZZ(BaseRing(R)));
    const BigInt P = characteristic(R);
    const long p = ConvertTo<long>(P);
    CoCoA_ASSERT(p < 32768); /// *** BUG BUG BUG ***

    FF Fp = FFctor(p);
    FFselect(Fp);
    const long d = StdDeg(f);
    DUPFF ans = DUPFFnew(d);
    { FFelem* CoeffVec = ans->coeffs; for (long i=0; i<=d; ++i) CoeffVec[i] = 0; }
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      BigInt c;
      if (!IsInteger(c, coeff(it))) { FFdtor(Fp); CoCoA_THROW_ERROR("coeff not integer","ConvertToDUPFF"); }
      long coeff = ConvertTo<long>(c);
      if (coeff < 0) coeff += p;
      const long exp = StdDeg(PP(it));
      ans->coeffs[exp] = coeff;
    }
    ans->deg = d;
///??? FFdtor(Fp); // AMB 2013-01-29: leaks, but o/w lots of "Invalid read of size 4"
    return ans;
  }


  RingElem ConvertDUPFFToRingElem(DUPFF F, ConstRefRingElem x)
  {
    const SparsePolyRing P = owner(x);
    PPMonoidElem t = LPP(x);
    RingElem ans(P);
    const int D = DUPFFdeg(F);
    for (int i=D; i >= 0; --i)
    {
      if (F->coeffs[i] == 0) continue;
      ans += monomial(P, F->coeffs[i], power(t,i));
    }
    return ans;
  }

/////?????  DUPZ ConvertToDUPZ(ConstRefRingElem f) { ... }

// // 2024-03-20 blindly copied from ConvertToDMPZ below (without really understanding)
// // Exception safe???  You're kidding, right?
//   DMPFF ConvertToDMPFF(ConstRefRingElem f)
//   {
//     CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
//     const SparsePolyRing P = owner(f);
//     const ring R = CoeffRing(P);
//     CoCoA_ASSERT(characteristic(R) != 0);
//     CoCoA_ASSERT(IsFiniteField(R));

//     const long nvars = NumIndets(P);
//     DMPFF g = NULL; // nullptr ???
//     std::vector<long> expv(nvars);
//     for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
//     {
// //      BigInt c;
// //      if (!IsInteger(c, coeff(it)))
// //        CoCoA_THROW_ERROR("coeff not in Z","factorize");
//       exponents(expv, PP(it));
//       int *exps = (int*)malloc(nvars*sizeof(int)); // block absorbed by g
//       std::copy(expv.begin(), expv.end(), exps);
//       g = DMPFFprepend(ConvertTo<long>(coeff(it)), exps, g);
//     }

//     return DMPFFreverse(g); /* restore original order of the terms */
//   }

// Exception safe???  You're kidding, right?
  DMPZ ConvertToDMPZ(ConstRefRingElem f)
  {
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing P = owner(f);
    const ring R = CoeffRing(P);
    CoCoA_ASSERT(characteristic(R) == 0);
    CoCoA_ASSERT(IsZZ(R) || IsFractionFieldOfGCDDomain(R));

    const long nvars = NumIndets(P);
    DMPZ g = NULL; // nullptr ???
    std::vector<long> expv(nvars);
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      BigInt c;
      if (!IsInteger(c, coeff(it)))
        CoCoA_THROW_ERROR("coeff not in Z","factorize");
      exponents(expv, PP(it));
      int *exps = (int*)malloc(nvars*sizeof(int)); // block absorbed by g
      std::copy(expv.begin(), expv.end(), exps);
      g = DMPZprepend(mpzref(c), exps, g);
    }

    return DMPZreverse(g); /* restore original order of the terms */
  }


// // 2024-03-20 blindly copied from ConvertDMPZToRingElem below
//   void ConvertDMPFFToRingElem(RingElem& f, DMPFF poly)
//   {
//     CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
//     const SparsePolyRing P = owner(f);
//     const PPMonoid M = PPM(P);
//     PPMonoidElem t(M);
//     f = 0;
//     const long nvars = NumIndets(P);
//     std::vector<long> exps(nvars);
//     for (DMPFF iter = poly; iter; iter = iter->next)
//     {
// //      if (mpz_sgn(iter->coeff) == 0) continue;
// //      const BigInt c = BigIntFromMPZ(iter->coeff);
//       std::copy(iter->exps, iter->exps+nvars, exps.begin());
//       f += monomial(P, iter->coeff, PPMonoidElem(M, exps));
//     }
//   }

  void ConvertDMPZToRingElem(RingElem& f, DMPZ poly)
  {
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing P = owner(f);
    const PPMonoid M = PPM(P);
    PPMonoidElem t(M);
    f = 0;
    const long nvars = NumIndets(P);
    std::vector<long> exps(nvars);
    for (DMPZ iter = poly; iter; iter = iter->next)
    {
      if (mpz_sgn(iter->coeff) == 0) continue;
      const BigInt c = BigIntFromMPZ(iter->coeff);
      std::copy(iter->exps, iter->exps+nvars, exps.begin());
      f += monomial(P, c, PPMonoidElem(M, exps));
    }
  }

  factorization<RingElem> ConvertToFactorList(DMPZfactors facpows, ring Rx)
  {
    RingElem content = RingElem(Rx, facpows->content); // converts from mpz_t
    RingElem factor(Rx);
    vector<RingElem> facs;
    vector<long> exps;

    DMPZlist iter;
    /*int i;*/
    for (/*i=1,*/ iter=facpows->list; iter; /*++i,*/ iter = iter->next)
    {
      ConvertDMPZToRingElem(factor, iter->poly);
      if (LC(factor) < 0) factor *= -1; // make LC positive
      facs.push_back(factor);
      exps.push_back(iter->deg);
    }
    DMPZfactors_dtor(facpows);  // NASTY - we destroy one of our args!!
    return factorization<RingElem>(facs, exps, content);
  }


  RingElem GCD_DUPFF(ConstRefRingElem f, ConstRefRingElem g)
  {
    CoCoA_ASSERT(owner(f) == owner(g));
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing P = owner(f);
    const vector<long> index = IndetsIn(f);
    const RingElem x = indet(P, index[0]);
    DUPFF F = ConvertToDUPFF(f);
    DUPFF G = ConvertToDUPFF(g);
    DUPFFgcd2(F,G);
    RingElem ans = ConvertDUPFFToRingElem(F, x);
    DUPFFfree(G);
    DUPFFfree(F);
    return ans;
  }


  // RingElem GCD_DMPFF(ConstRefRingElem f, ConstRefRingElem g)
  // {
  //   CoCoA_ASSERT(owner(f) == owner(g));
  //   CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
  //   const SparsePolyRing P = owner(f);
  //   const vector<long> index = IndetsIn(f);
  //   const RingElem x = indet(P, index[0]);
  //   DMPFF F = ConvertToDMPFF(f);
  //   DMPFF G = ConvertToDMPFF(g);
  //   DMPFF H = DMPFFgcd(F,G);
  //   RingElem ans(P);
  //   ConvertDMPFFToRingElem(ans, H);
  //   DMPFFdtor(H);
  //   DMPFFdtor(G);
  //   DMPFFdtor(F);
  //   return ans;
  // }


  RingElem GCD_DMPZ(ConstRefRingElem f, ConstRefRingElem g)
  {
    CoCoA_ASSERT(owner(f) == owner(g));
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing P = owner(f);
    CoCoA_ASSERT(IsZZ(CoeffRing(P)) || IsQQ(CoeffRing(P)));

    if (IsQQ(CoeffRing(P)))
    {
      // case P = QQ[x,y,z,...]
      RingHom coeff = CoeffEmbeddingHom(P);
      RingElem ContF = content(f);
      RingElem ContG = content(g);
      NVARS = NumIndets(P); // BUG BUG BUG:  NVARS is a global!!!
      DMPZ F = ConvertToDMPZ(f/coeff(ContF));
      DMPZ G = ConvertToDMPZ(g/coeff(ContG));
      DMPZ H = DMPZgcd(F,G);
      DMPZdtor(G);
      DMPZdtor(F);

      RingElem ans(P);
      ConvertDMPZToRingElem(ans, H);
      DMPZdtor(H);
      return ans;
    }
    if (IsZZ(CoeffRing(P)))
    {
      // case P = ZZ[x,y,z,...]
      NVARS = NumIndets(P); // BUG BUG BUG:  NVARS is a global!!!
      DMPZ F = ConvertToDMPZ(f);
      DMPZ G = ConvertToDMPZ(g);
      DMPZ H = DMPZgcd(F,G);
      DMPZdtor(G);
      DMPZdtor(F);

      RingElem ans(P);
      ConvertDMPZToRingElem(ans, H);
      DMPZdtor(H);
      return ans;
    }

    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "GCD_DMPZ");
    return f; // just to keep compiler quiet
  }


  factorization<RingElem> factor(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (IsPolyRing(R) && IsQuotientRing(CoeffRing(R))
        && IsPolyRing(BaseRing(CoeffRing(R))) && IsField(CoeffRing(R)))
      return factor_AlgExtn(f);
    if (IsField(R) || !IsTrueGCDDomain(R))
      CoCoA_THROW_ERROR("factorization not supported", "factor(RingElem)");
    if (IsZero(f))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "factor(x)");
    if (IsZZ(R))
    {
      // Factorization is actually over ZZ, so use factor(BigInt) to do the work!
      const factorization<BigInt> IntFacs = factor(ConvertTo<BigInt>(f));
      const vector<BigInt>& prime = IntFacs.myFactors();
      const long nprimes = len(prime);
      vector<RingElem> facs;
      for (long i=0; i < nprimes; ++i)
        facs.push_back(RingElem(R, prime[i]));
      return factorization<RingElem>(facs, IntFacs.myMultiplicities(), RingElem(R, IntFacs.myRemainingFactor()));
    }
    if (!IsPolyRing(R))
      CoCoA_THROW_ERROR(ERR::ReqPolyRing, "factor(RingElem)");
/////    if (!IsSparsePolyRing(R))
/////      CoCoA_THROW_ERROR(ERR::NYI, "factor(RingElem)  NYI over DUP");
    const SparsePolyRing P = owner(f);
    if (IsQQ(CoeffRing(P)))
    {
      // case P = QQ[x,y,z,...]
      const FractionField QQ = CoeffRing(P);
      const RingHom coeff = CoeffEmbeddingHom(P)(EmbeddingHom(QQ));
      const RingElem D = CommonDenom(f);
      NVARS = NumIndets(P); // NVARS is a global!!!
      const DMPZ g = ConvertToDMPZ(f*coeff(D));
      factorization<RingElem> ans = ConvertToFactorList(DMPZfactor(g), P); // BUG BUG BUG leaks the factorlist
      DMPZdtor(g);
      // Fix sign of RemainingFactor -- DMPZfactor sometimes goes wrong :-(
      if (sign(LC(f)) != sign(LC(ans.myRemainingFactor()))) ans.myNewRemainingFactor(-ans.myRemainingFactor());
      if (!IsOne(D)) ans.myNewRemainingFactor(ans.myRemainingFactor()/CanonicalHom(RingZZ(), P)(D));
      return ans;
    }
    if (IsZZ(CoeffRing(P)))
    {
      // We do not attempt to factorize the content in ZZ[x,y,z] as it could easily be very costly.
      PolyRing QQx = NewPolyRing(RingQQ(), NewSymbols(NumIndets(P)));
      RingHom phi = PolyRingHom(P, QQx, ZZEmbeddingHom(RingQQ()), indets(QQx)); // maps ZZ[x,y,z] --> QQ[x,y,z]
      RingElem F_in_QQx = phi(f);
      factorization<RingElem> Factorization_in_QQx = factor(F_in_QQx);
      const vector<RingElem>& Facs_in_QQx = Factorization_in_QQx.myFactors();
      vector<RingElem> FacsInP;
      const RingElem ContentInP = CoeffEmbeddingHom(P)(num(LC(Factorization_in_QQx.myRemainingFactor())));
      RingHom psi = PolyRingHom(QQx, P, QQEmbeddingHom(RingZZ()), indets(P)); // maps QQ[x,y,z] --> ZZ[x,y,z]
      const long NumFacs = len(Facs_in_QQx);
      for (long i=0; i < NumFacs; ++i)
      {
        FacsInP.push_back(psi(Facs_in_QQx[i]));
      }
      return factorization<RingElem>(FacsInP, Factorization_in_QQx.myMultiplicities(), ContentInP);
    }
    if (IsFiniteField(CoeffRing(P)))
    {
      // Case FF(q)[x,y,z,...]
      ring FFq = CoeffRing(P);
      if (LogCardinality(FFq) > 1)
        CoCoA_THROW_ERROR(ERR::NYI, "Factorization over alg extn of finite field");
      long var;
      if ((var = UnivariateIndetIndex(f)) >= 0)
      {
        const RingElem x = indet(P, var);
        SmallFpImpl ModP(ConvertTo<long>(characteristic(P)));
        DUPFp F = ConvertToDUPFp(ModP, f);
        const factorization<DUPFp> facs = factor(F);
        const vector<DUPFp>& DUPfac = facs.myFactors(); // convenient alias
        const long n = len(facs.myFactors());
        vector<RingElem> irreds; irreds.reserve(n);
        vector<long> multiplicity; multiplicity.reserve(n);
        for (long i=0; i < n; ++i)
        {
          irreds.push_back(ConvertFromDUPFp(x, DUPfac[i]));
        }
        return factorization<RingElem>(irreds, facs.myMultiplicities(), CoeffEmbeddingHom(P)(LC(f)));

//         DUPFF F = ConvertToDUPFF(f);
//         DUPFFlist facs = DUPFFfactor(F);
//         DUPFFfree(F);
//         const long n = DUPFFlist_length(facs);
//         vector<RingElem> irreds; irreds.reserve(n);
//         vector<long> multiplicity; multiplicity.reserve(n);
//         for (DUPFFlist curr=facs; curr != nullptr; curr = curr->next)
//         {
//           irreds.push_back(monic(ConvertDUPFFToRingElem(curr->poly, x)));
//           multiplicity.push_back(curr->deg);
//         }
//         DUPFFlist_dtor(facs);
//         return factorization<RingElem>(irreds, multiplicity, CoeffEmbeddingHom(P)(LC(f)));
      }
      CoCoA_THROW_ERROR(ERR::NYI, "Multivariate factorization over finite field");
    }
    CoCoA_THROW_ERROR(ERR::NYI, "Factorization not in QQ[x,y,...]  or in ZZ/(p)[x,y,...]");
    return factorization<RingElem>(vector<RingElem>(),vector<long>(),RingElem(P)); // NEVER EXECUTED, just to keep compiler quiet
  }


  // temporary: by AMB 2008-01-25
  bool IsIrredPoly(ConstRefRingElem f)
  { 
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqPolyRing, "IsIrredPoly");
 // ??? warning: factorize syntax might change: now only for multivariate
    if (IsConstant(f)) return IsIrred(LC(f));
    // ANNA: should do it for ZZ as well (taking care of constant factors)
    const PolyRing P = owner(f);
    if (!IsZZ(CoeffRing(P)) && !IsField(CoeffRing(P)))
      CoCoA_THROW_ERROR(ERR::NYI, "IsIrredPoly: CoeffRing must be ZZ or a field");
    if ( IsSparsePolyRing(P) )
    {
      const factorization<RingElem> facs = factor(f);
      const long NumFacs = len(facs.myFactors());
      const long mult = facs.myMultiplicities()[0];
      return (NumFacs == 1 && mult == 1);  // see  [redmine 1710]
//???      return (IsInvertible(facs.myRemainingFactor()) && NumFacs == 1 && mult == 1); // see  [redmine 1710]
    }
    else // 2022-06-15 ??? presumably DenseUnivariate if we get here
    {
      PolyRing SPR = NewPolyRing(CoeffRing(P), NewSymbols(NumIndets(P)));
      RingHom phi(PolyAlgebraHom(P, SPR, indets(SPR)));
      const factorization<RingElem> facs = factor(phi(f));
      const long NumFacs = len(facs.myFactors());
      const long mult = facs.myMultiplicities()[0];
      return (IsInvertible(facs.myRemainingFactor()) && NumFacs == 1 && mult == 1);
    }

  }


  RingElem RadicalOfPoly(ConstRefRingElem f)
  {
    const ring& P = owner(f);
    const ring& k = CoeffRing(P);
    if (!IsZZ(k) && !IsField(k)) CoCoA_THROW_ERROR(ERR::ReqField, "RadicalOfPoly");
    if (IsZero(f)) return f;
    const bool OverZZ = IsZZ(k);
    const BigInt PseudoContent = OverZZ ? ConvertTo<BigInt>(content(f)) : BigInt(1);
    const BigInt RadicalOfContent = radical(PseudoContent);
    // Special case: constant
    if (IsConstant(f))
      return RadicalOfContent*one(P);
    // Special case: monomial
    if (IsMonomial(f)) return monomial(P, RadicalOfContent, radical(LPP(f)));
    // General case: poly has at least 2 terms.
    const bool CharZero = (characteristic(P) == 0);
    if (!CharZero)
    {
      // Non-zero characteristic:
      // this way is simple, and apparently quite quick in my tests...
      RingElem ans = one(P);
      const factorization<RingElem> FacInfo = SqFreeFactor(f);
      const vector<RingElem>& facs = FacInfo.myFactors();
      const int nfacs = len(facs);
      for (int j=0; j < nfacs; ++j)
        ans *= facs[j];
      return ans;
    }
    // Zero characteristic
#if 1
    // First impl: do ContentFreeFactor, then on each factor do SqFreeFactor... [!FASTER!]
    const factorization<RingElem> FacInfo = ContentFreeFactor(f/PseudoContent);
    const vector<RingElem>& facs = FacInfo.myFactors();
    const int nfacs = len(facs);
    RingElem ans = one(P);
    for (int i=0; i < nfacs; ++i)
    {
      // Pick an indet appearing in facs[i]; currently, just take first we find.
      int x=0; // ok if univariate
      if (NumIndets(P) > 1)
      {
        // must choose if multivariate...
        const vector<long> expv = exponents(LPP(facs[i]));
        // loop below is safe because we know that f is not constant (see above)
        while (expv[x] == 0)
        {
          ++x;
          CoCoA_ASSERT(x < len(expv));
        }
      }
      ans *= facs[i]/monic(gcd(facs[i], deriv(facs[i], x)));
    }
#else
    // Second impl: just do SqFreeFactor then take product of each factors. [!SLOWER!]
    RingElem g = f/PseudoContent;
    RingElem ans = one(P);
    while (!IsConstant(g))
    {
      // Pick an indet appearing in facs[i]; currently, just take first we find.
      int x=0; // ok if univariate
      if (NumIndets(P) > 1)
      {
        // must choose if multivariate...
        const vector<long> expv = exponents(LPP(g));
        while (expv[x] == 0)
        {
          ++x;
          CoCoA_ASSERT(x < len(expv));
        }
      }
      const RingElem contentg = ContentWRT(g, indet(P,x));
      g /= contentg;
      ans *= g/monic(gcd(g, deriv(g, x)));
      g = contentg; // better to swap???
    }
    #endif
    return RadicalOfContent*ans;
  }


  //-------------------------------------------------------

  factorization<RingElem> SqFreeFactor_generic(ConstRefRingElem x)
  {
    return factor(x);
  }


  // Simple rather than super-efficient.
  RingElem IteratedPthRoot(RingElem f)
  {
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    while (IsPthPower(f))
      f = PthRoot(f);
    return f;
  }

  factorization<RingElem> SqFreeFactorPosDerChar0(ConstRefRingElem f, long x)
  {
    const bool OverField = IsField(CoeffRing(owner(f)));
    const RingElem derivf = deriv(f,x);
    RingElem RemainingFacs = gcd(f,derivf);
    if (OverField)
      RemainingFacs = monic(RemainingFacs);
    RingElem rad = f/RemainingFacs;
    RingElem DerivCofactor = derivf/RemainingFacs;

    RingElem U = DerivCofactor - deriv(rad, x);
    long i = 1;

    vector<RingElem> facs;
    vector<long> exps;
    while (!IsZero(U))
    {
      /*const*/ RingElem G = gcd(rad, U); if (OverField) G = monic(G);
      if (!IsConstant(G))
      {
        facs.push_back(G);
        exps.push_back(i);
      }
      ++i;
      rad /= G;
      RemainingFacs /= rad;
      DerivCofactor = U/G;
      U = DerivCofactor - deriv(rad,x);
    }
    if (!IsConstant(rad))
    { if (OverField) rad = monic(rad); facs.push_back(rad); exps.push_back(i); }

    return factorization<RingElem>(facs, exps, RemainingFacs);
  }


  factorization<RingElem> SqFreeFactorChar0(RingElem f)
  {
    // Assume CoeffRing is Field or ZZ
    const ring& P = owner(f);
    const long nvars = NumIndets(P);
    const RingElem PseudoContent = IsField(CoeffRing(P)) ? LC(f) : sign(LC(f))*content(f); // sign(LC(f)) -- see [[redmine 1710]]
    f = f/PseudoContent;
    vector<RingElem> facs;
    vector<long> exps;
    for (long x=0; x < nvars; ++x)
    {
///???      if (deg(f,x) == 0) continue;
      if (deg(f,x) == 0) continue;
      const factorization<RingElem> PartialFacs = SqFreeFactorPosDerChar0(f, x);
      f = PartialFacs.myRemainingFactor();
      facs.insert(facs.end(), PartialFacs.myFactors().begin(), PartialFacs.myFactors().end());
      exps.insert(exps.end(), PartialFacs.myMultiplicities().begin(), PartialFacs.myMultiplicities().end());
    }
    return factorization<RingElem>(facs, exps, CoeffEmbeddingHom(P)(PseudoContent));
  }


  factorization<RingElem> SqFreeFactorPosDerCharP(ConstRefRingElem f, long x)
  {
    const BigInt P = characteristic(owner(f));
    const RingElem derivf = deriv(f,x);
    RingElem RemainingFacs = gcd(f,derivf);
    RingElem rad = f/RemainingFacs;
    RingElem DerivCofactor = derivf/RemainingFacs;

    RingElem U = DerivCofactor - deriv(rad, x);
    long i = 1;

    vector<RingElem> facs;
    vector<long> exps;
    while (i < P && !IsZero(U))
    {
      const RingElem G = monic(gcd(rad, U));
      if (!IsConstant(G))
      {
        facs.push_back(G);
        exps.push_back(i);
      }
      ++i;
      rad /= G;
      RemainingFacs /= rad;
      DerivCofactor = U/G;
      U = DerivCofactor - deriv(rad,x);
    }
    if (!IsConstant(rad))
    { facs.push_back(monic(rad)); exps.push_back(i); }

    return factorization<RingElem>(facs, exps, RemainingFacs);
  }


  factorization<RingElem> SqFreeFactorCharP(RingElem f)
  {
    const BigInt P = characteristic(owner(f));
    const long NumVars = NumIndets(owner(f));
    const RingElem LCF = CoeffEmbeddingHom(owner(f))(LC(f));
    f /= LCF; // good idea or bad???

    vector<RingElem> facs;
    vector<long> mults;
    for (long x=0; x < NumVars; ++x)
    {
////?????      if (deg(f,x) == 0) continue;
      if (deg(f,x) == 0) continue;
      const factorization<RingElem> tmp = SqFreeFactorPosDerCharP(f, x);
      facs.insert(facs.end(), tmp.myFactors().begin(), tmp.myFactors().end());
      mults.insert(mults.end(), tmp.myMultiplicities().begin(), tmp.myMultiplicities().end());
      f = tmp.myRemainingFactor();
    }

    if (deg(f) == 0)
    {
      return factorization<RingElem>(facs, mults, LCF);
    }
    // else deg(f) > 0, so f is a p-th power.
    long q = deg(f);
    f = IteratedPthRoot(f);
    q /= deg(f); // a power of p
    if (q == 1) { facs.push_back(f); mults.push_back(1); return factorization<RingElem>(facs, mults, LCF); } 
    const factorization<RingElem> RootFacs = SqFreeFactorCharP(f);
    vector<RingElem> PthRootFacs = RootFacs.myFactors(); // QUICK HACK now that myFactors is read-only!
    // Now create a sort of CoprimeFactorBasis; we know that the entries of facs
    // are coprime, and also RootsFacs.myFactors are coprime.
    const long n = len(facs);
    for (long i=0; i < n; ++i)
    {
      for (long j=0; j < len(PthRootFacs); ++j)
      {
        if (IsConstant(PthRootFacs[j])) continue;
        const RingElem Q = monic(gcd(facs[i], PthRootFacs[j]));
        if (IsConstant(Q)) continue;
        facs[i] /= Q;
        PthRootFacs[j] /= Q;
        facs.push_back(Q);
        mults.push_back(mults[i]+q*RootFacs.myMultiplicities()[j]);
        if (IsConstant(facs[i])) break;
      }
    }

    // Remove all constant factors
    long last = len(facs)-1;
    while (last >= 0 && IsConstant(facs[last]))
      --last;
    for (long i=0; i < last; ++i)
    {
      if (!IsConstant(facs[i])) continue;
      facs[i] = facs[last]; //swap(facs[i], facs[last]);
      mults[i] = mults[last]; //swap(mults[i], mults[last]);
      --last;
    }
    facs.resize(last+1);
    mults.resize(last+1);


    // Append any (non-constant) new factors
    for (long j=0; j < len(PthRootFacs); ++j)
    {
      if (IsConstant(PthRootFacs[j])) continue;
      facs.push_back(PthRootFacs[j]);
      mults.push_back(q*RootFacs.myMultiplicities()[j]);
    }
    return factorization<RingElem>(facs, mults, LCF);
  }


  factorization<RingElem> SqFreeFactor(ConstRefRingElem f)
  {
    const char* const FnName = "SqFreeFactor";
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR(ERR::NYI, FnName);
    if (!IsTrueGCDDomain(P))
      CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, FnName);
    if (IsZero(f))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, FnName);

    const ring& k = CoeffRing(P);
    //    if (IsQQ(k)) return SqFreeFactorChar0(f);  // ANNA: true for Ch0???
    if (IsZero(characteristic(k))) return SqFreeFactorChar0(f);
    if (IsFiniteField(k)) return SqFreeFactorCharP(f);
    if (IsPrime(characteristic(k))) return SqFreeFactorCharP(f); // BUG BUG  DODGY!!!
    return SqFreeFactor_generic(f);
  }


  namespace // anonymous
  {
    void ContentFreeFactorLoop(vector<RingElem>& ans, RingElem f, vector<bool> skip)
    {
      CoCoA_ASSERT(!IsZero(f));
      CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
      const SparsePolyRing P = owner(f);
      CoCoA_ASSERT(IsTrueGCDDomain(CoeffRing(P)) || IsField(CoeffRing(P))); // coeffs are in field or GCD domain
      const vector<long> var = IndetsIn(f);
      if (var.empty()) return;
      for (int i=0; i < len(var); ++i)
      {
        if (skip[var[i]]) continue;
        RingElem c = ContentWRT(f, indet(P,var[i]));  // not const so we can use std:move
        if (IsOne(c)) continue;
        f /= c;
        skip[i] = true;
        ContentFreeFactorLoop(ans, std::move(c), skip);
      }
      if (StdDeg(f) == 0) return;
      // f is non-trivial; normalize it before adjoining it to the result.
      if (IsField(CoeffRing(P)))
      {
        ans.push_back(monic(f));
        return;
      }
      // otherwise IsTrueGCDDomain(CoeffRing(P))
      {
        const RingHom embed = CoeffEmbeddingHom(P);
        ans.push_back(f/embed(content(f)));
        return;
      }
    }
  } // end of anonymous namespace

  factorization<RingElem> ContentFreeFactor(ConstRefRingElem f)
  {
    const char* const FnName = "ContentFreeFactor";
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName);
    const SparsePolyRing P = owner(f);
    const ring& k = CoeffRing(P);
    if (!IsZZ(k) && !IsField(k) && !IsTrueGCDDomain(k))
      CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, FnName); //???? what error to give here???
    if (IsZero(f))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, FnName);

    const RingElem PseudoContent = (IsZZ(k)) ? CoeffEmbeddingHom(P)(content(f)) : CoeffEmbeddingHom(P)(LC(f));
    if (NumIndets(P) == 1) { factorization<RingElem> ans(PseudoContent); ans.myAppend(f/PseudoContent, 1); return ans; }
    // Idea: could remove numerical content here... (is it worth it?)
    RingElem RemainingFactor = PseudoContent;//LC(f); // BEFORE calling ContentFreeFactorLoop
    vector<RingElem> facs;
    ContentFreeFactorLoop(facs, f/PseudoContent, vector<bool>(NumIndets(P)));
    const vector<long> exps(len(facs), 1);
    // if (!IsZZ(k))  // I think the factors are monic if coeffring is field
    //   for (int i=0; i < len(facs); ++i) // Is This loop useful???
    //     RemainingFactor /= LC(facs[i]);
    return factorization<RingElem>(facs, exps, RemainingFactor);
  }


  RingElem PalindromicFactor(ConstRefRingElem f)
  {
    // assume input is (non-const) univariate
    if (IsPalindromic(f))  return f;  // short-cut because gcd(f,f) is disappointingly slow (2023-07-19)
    return gcd(f, reverse(f));
  }


} // end of namespace CoCoA

// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/factor.C,v 1.41 2024/07/11 15:44:50 bigatti Exp $
// $Log: factor.C,v $
// Revision 1.41  2024/07/11 15:44:50  bigatti
// Summary: error codes: NotField to ReqField
//
// Revision 1.40  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.39  2024/07/02 12:48:13  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.38  2024/03/21 20:07:50  abbott
// Summary: Added some cruft (?!? commented out code -- might be useful one day?)
//
// Revision 1.37  2024/02/29 14:38:29  abbott
// Summary: Commented out apparently useless/unused variable
//
// Revision 1.36  2023/07/26 11:57:21  abbott
// Summary: Minor improvement to PalindromicFactor
//
// Revision 1.35  2023/03/09 22:41:07  abbott
// Summary: Tidied use of exponents (redmine 1694)
//
// Revision 1.34  2023/02/23 20:51:14  abbott
// Summary: Added new fn PalindromicFactor
//
// Revision 1.33  2022/11/16 19:47:42  abbott
// Summary: Modified semantics of IsIrredPoly (redmine 1710); SqFreeFactor allow negative RemainingFactor, so all real factors have positive LC
//
// Revision 1.32  2022/08/08 17:56:26  abbott
// Summary: Revised several fns to work also over ZZ (SqFreeFactor, radical, ...)
//
// Revision 1.31  2022/06/08 15:24:08  abbott
// Summary: IsIrredPoly now works over ZZ too (redmine 1683)
//
// Revision 1.30  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.29  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.28  2020/10/05 19:26:08  abbott
// Summary: Used std:move inside ContentFreeFactorLoop
//
// Revision 1.27  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.26  2020/02/14 12:23:15  abbott
// Summary: Changed fn call name: old "indets", new "IndetsIn"
//
// Revision 1.25  2019/11/14 17:53:39  abbott
// Summary: Removed commented out logging to clog
//
// Revision 1.24  2019/03/19 11:07:08  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.23  2019/03/18 11:22:32  abbott
// Summary: Changed include (after splitting NumTheory)
//
// Revision 1.22  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.21  2018/06/12 13:57:11  abbott
// Summary: Fixed sign of RemainingFactor in multivariate factorization (see redmine #1185)
//
// Revision 1.20  2018/05/22 14:16:40  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.19  2018/05/18 16:39:39  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.18  2018/05/18 12:25:54  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.17  2018/05/17 16:03:04  bigatti
// -- added #include SparsePolyIter
//
// Revision 1.16  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.15  2018/04/19 13:38:10  bigatti
// -- added factor_AlgExtn
//
// Revision 1.14  2018/02/27 17:30:23  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.13  2018/02/27 10:53:51  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.12  2017/09/06 11:56:30  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.11  2017/06/22 11:03:33  abbott
// Summary: Added misin call to LC when converting content over QQ to content over ZZ
//
// Revision 1.10  2016/11/08 13:08:11  bigatti
// -- RadicalOfPoly: changed error mesg for NotField,
//                   changed NumTerms==1 to IsMonomial
//
// Revision 1.9  2016/10/20 18:03:13  abbott
// Summary: Added RadicalOfPoly
//
// Revision 1.8  2016/01/26 14:09:35  bigatti
// -- changed SqFreeFactor QQ --> char 0
//
// Revision 1.7  2014/07/08 15:24:23  abbott
// Summary: Removed AsQuotientRing from an assertion
// Author: JAA
//
// Revision 1.6  2014/07/08 08:40:18  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.5  2014/07/07 14:12:35  abbott
// Summary: Corrected silly typo
// Author: JAA
//
// Revision 1.4  2014/07/07 13:24:48  abbott
// Summary: Removed AsPolyRing; Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.3  2014/05/06 13:20:41  abbott
// Summary: Changed names (my)MaxExponents into (my)Deg
// Author: JAA
//
// Revision 1.2  2014/04/30 16:26:49  abbott
// Summary: Replaced size_t by long; Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.1  2014/04/11 15:06:38  abbott
// Summary: Renamed from TmpFactor.H/C
// Author: JAA
//
// Revision 1.31  2014/04/08 16:35:30  abbott
// Summary: Changed SqfreeFactor to SqFreeFactor
// Author: JAA
//
// Revision 1.30  2014/03/24 12:09:21  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.29  2014/01/16 16:15:13  abbott
// Renamed  SqfrDecomp_generic  into  SqfrFactor_generic  & cleaned the code slightly.
//
// Revision 1.28  2013/10/22 14:02:09  abbott
// Added code for SqFreeFactor.
//
// Revision 1.27  2013/09/12 12:59:52  abbott
// Minor cleaning of factorize fn when applied to a BigInt.
//
// Revision 1.26  2013/04/16 15:51:00  bigatti
// -- added "leak" warning comment
//
// Revision 1.25  2013/02/21 12:55:23  abbott
// Improved factor: now works in Fp[x,y,z] provided given poly is univariate.
//
// Revision 1.24  2013/01/28 13:05:22  abbott
// Fixed two memory leaks in GCD_DMPZ (found by valgrind).
//
// Revision 1.23  2012/10/15 12:38:12  abbott
// Corrected signature for ContentFreeFactorLoop.
//
// Revision 1.22  2012/10/05 09:32:12  abbott
// Changed myExponents into myMultiplicities.
//
// Revision 1.21  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.20  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.19  2012/05/04 20:02:20  abbott
// Added brackets as suggested by newish version of g++.
//
// Revision 1.18  2012/04/27 15:10:00  abbott
// Cleaned ConvertDUPFFToRingElem; inproved factor fn.
//
// Revision 1.17  2012/04/18 14:26:29  abbott
// Corrected an incorrect assertion (found via anna.cocoa5 test).
//
// Revision 1.16  2012/04/15 20:14:33  abbott
// Added several fns:
//   temporary hacks GCD_DUPFF and GCD_DMPZ to use the old CoCoA-4 code
//   ContentFreeFactor.
//
// Revision 1.15  2012/02/24 13:47:29  abbott
// Temporary check-in of a version which compiles (but new fns DO NOT YET WORK).
//
// Revision 1.14  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.13  2012/02/08 15:14:30  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.12  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.11  2011/06/27 12:56:37  bigatti
// -- just a comment Anna --> AMB
//
// Revision 1.10  2009/12/03 17:41:20  abbott
// Added include of cstdlib for malloc.
//
// Revision 1.9  2009/10/26 15:43:25  bigatti
// -- added CopyFromMPZ in ZZ ctor
//
// Revision 1.8  2009/07/24 14:22:40  abbott
// Changed interface to factorizer: new function name (now "factor") with new signature.
// Did some cleaning too.
//
// Revision 1.7  2009/07/02 16:34:32  abbott
// Minor change to keep the compiler quiet.
//
// Revision 1.6  2009/05/14 09:20:03  abbott
// Changed a comment.
//
// Revision 1.5  2008/10/07 15:46:10  abbott
// Added missing log/history lines at the ends of a few files.
//
//
