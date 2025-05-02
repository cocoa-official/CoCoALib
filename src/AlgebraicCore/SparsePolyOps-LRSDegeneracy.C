//   Copyright (c)  2022,2025  John Abbott and Anna M. Bigatti
//   Original author:  Nico Mexis  (2022)

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


#include "CoCoA/SparsePolyOps-LRSDegeneracy.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/factor.H"
#include "CoCoA/factorization.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-misc.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/NumTheory-PrimeModSeq.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/random.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-cyclotomic.H"
#include "CoCoA/SparsePolyOps-Graeffe.H"
#include "CoCoA/SparsePolyOps-resultant.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/verbose.H"
#include "CoCoA/VerificationLevel.H"

#include <vector>
using std::vector;
#include <algorithm>
using std::min;

namespace CoCoA
{

  namespace // anonymous
  {

    // Check input poly: throw if not univariate over char 0
    // Returns the radical (divided by x if possible) -- if this has deg < 2 then obv not LRSdegen
    RingElem LRSDegenerate_check_and_preprocess(ConstRefRingElem f, const ErrorContext& OrigContext)
    {
      const ring &Px = owner(f);
      if (!IsSparsePolyRing(Px))
        CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqSparsePolyRing, OrigContext);
      if (!IsZero(characteristic(Px)))
        CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqChar0, OrigContext);
      if (IsConstant(f))
        CoCoA_THROW_ERROR_WITH_CONTEXT2("Polynomial must be non-constant", OrigContext);
      const long IndetIndex = UnivariateIndetIndex(f);
      if (IndetIndex < 0)
        CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqUnivariate, OrigContext);

      // Divide out factors of x:
      RingElem f_reduced = prim(f);
      if (IsZero(ConstantCoeff(f)))
        f_reduced /= gcd(f, IndetPower(Px, IndetIndex, deg(f)));

      // Extract square-free part
      f_reduced /= gcd(f_reduced, deriv(f_reduced, IndetIndex));

      return f_reduced;
    }


    // simple impl -- not (yet) efficient
    RingElem ScaleX(ConstRefRingElem f, ConstRefRingElem zeta)
    {
      CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
      RingElem g = zero(owner(f));
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      {
        PushBack(g, power(zeta,deg(PP(it)))*coeff(it), PP(it));
      }
      return g;
    }

    /**
     * Used in modular approach to decide whether
     * the n-th power Graeffe transform of f needs
     * to be checked for proving its LRS-degeneracy.
     * For private use.  Hence, no arg checking.
     */
    bool IsLRSDegenerateOrder_ModularCheck(ConstRefRingElem f, const unsigned long n, const long NumPrimesToTry)
    {
      // ASSUME: f univariate, sqfree, f(0) non-zero, deg(f) >= 2, coeffs integer
      CoCoA_ASSERT(n > 2);
      const char* const FnName = "IsLRSDegenerateOrder_ModularCheck";
      VerboseLog VERBOSE(FnName);

      const unsigned long degf = deg(f); // **ASSUME**  8*degf^2 fits into long
      const BigInt LCf = ConvertTo<BigInt>(LC(f));

      const unsigned long gcdBound = n/2;

      static const auto X = NewSymbols(1); // do this once to avoid "burning through" lots of anon symbs

      const unsigned long PrimeLWB = std::max(8*degf*degf, 64*n);
      PrimeSeq1ModN primeSeq(n);  while (*primeSeq < PrimeLWB)  { ++primeSeq; } // skip past "very small primes"
      for (long i = 0; i < NumPrimesToTry; ++i)
      {
        // Jump forward a "random number" of suitable primes
        SmallPrime p = *primeSeq;
        for (int skip = RandomLong(1,10); skip > 0; --skip)  // GUESS suitable range (wider???)
        {
          while (p != 0)  { ++primeSeq; p = *primeSeq; if (p == 0 || LCf%p != 0)  break;}
        }

        if (IsZero(p)) // upper prime bound has been reached
          break;

        // Next commented-out block was orig code, but call to the scaling
        // PolyAlgebraHom was too costly
//         // zeta will be a primitive root of unity in F_p
//         const ring& Px = owner(f);
//         const ring& P = CoeffRing(Px);
//         const QuotientRing Fp = NewZZmod(p);
//         const SparsePolyRing Fpx = NewPolyRing(Fp, X);
//         const RingElem& x = indet(Fpx,0);
//         const std::vector<RingElem> images(NumIndets(Px), x);
//         const RingElem f_p = PolyRingHom(Px, Fpx, CanonicalHom(P, Fp), images)(f);
//         const unsigned long zeta = PowerMod(PrimitiveRoot(p), (p-1)/n, p);

//         for (unsigned long k = 1; k <= gcdBound; ++k)
//         {
//           CheckForInterrupt(FnName);
//           if (!IsCoprime(n, k)) continue;
// //          if (IsConstant(gcd(f_p, PolyAlgebraHom(Fpx, Fpx, {PowerMod(zeta, k, p) * x})(f_p)))) // MOST COSTLY STEP (according to valgrind/callgrind)
//           if (IsConstant(gcd(f_p, ScaleX(f_p, RingElem(Fp,PowerMod(zeta, k, p))))))
//           {
//             return false;
//           }
//         }

          // Equiv to commented-out block, but faster (using DUPFp)
          const SmallFpImpl ModP(p);
          const DUPFp f_p = ConvertToDUPFp(ModP, f);
          const unsigned long zeta = PowerMod(PrimitiveRoot(p), (p-1)/n, p);
          for (unsigned long k = 1; k <= gcdBound; ++k)
          {
            CheckForInterrupt(FnName);
            if (!IsCoprime(n, k)) continue;
            const DUPFp g_p = ScaleX(f_p, ModP.myReduce(PowerMod(zeta, k, p)));
            if (IsConstant(gcd(f_p, g_p)))  return false;
          }
          ++primeSeq;
        }

      return true;
    }



    /**
     * Used in modular approach to decide whether
     * the n-th power Graeffe transform of f needs
     * to be checked for proving its LRS-degeneracy.
     * For private use.  Hence, no arg checking.
     */

    // Input: D is the degree (assumed >= 2)
    // Return list of all indices k s.t. cyclo(k) could divide LRS test resultant
    std::vector<unsigned long> GetKToTryList(const unsigned long D/*, const unsigned long Kbound*/)
    {
      CoCoA_ASSERT(D > 1);
      CoCoA_ASSERT(std::numeric_limits<unsigned long>::max()/D > D); // D*D does not overflow

      const unsigned long Kbound = InvTotientBoundUpto_ulong(D*(D-1));  // conjecture: can use InvTotientBoundUpto_ulong(D) if we want just the first K
      // Create list of values which phi(k) must divide
      // List is repr as a "boolean mask": mask[d] == 1 iff d is in list
      std::vector<unsigned char> mask(D*D);  // vector<bool>  ???
      for (unsigned long d=2; d <= D; ++d)
        mask[d*(d-1)] = 1;

      for (unsigned long a = 1; a < D; ++a)
      {
        unsigned long UPB = min(a, D-a);
        for (unsigned long b = 1; b <= UPB; ++b)
        {
          if (IsEven(a * b))
            mask[a * b] = 1;
        }
      }

      // AllEvenFactorsMask is really a boolean array:
      // entry d is true iff d is even and divides one of the values in "mask"
      std::vector<unsigned char> AllEvenFactorsMask(D*D);  // vector<bool>  ???
      for (unsigned long d=D*(D-1); d >= 2; d -= 2) // start from biggest, only even values.
      {
        if (mask[d] == 0) continue;
        if (AllEvenFactorsMask[d] != 0) continue;
        // Now mark all even factors of d:
        const unsigned long half = d/2;
        for (unsigned long fac=1; fac*fac <= half; ++fac)
          if (half%fac == 0)
          { AllEvenFactorsMask[2*fac] = 1; AllEvenFactorsMask[d/fac] = 1; }
      }

      // Find all k values such that phi(k) divides some elem of NewL
      std::vector<unsigned long> KToTry;
      for (unsigned long k = 3; k <= Kbound; ++k)
      {
        const unsigned long phi = EulerTotient(k);
        if (phi < D*D && AllEvenFactorsMask[phi] != 0)
          KToTry.push_back(k);
      }
      return KToTry;
    }



//     /**
//      * Checks for LRS-degeneracy of the given order, without arg checks
//      */
//     bool IsLRSDegenerateOrderMod_NoArgChecks(ConstRefRingElem f, const unsigned long n, const long IndetIndex, VerificationLevel VerifLev)
//     {
//       // ASSUME: f is elem of ZZ[x], sqfr, f(0) != 0, deg(f) >= 2
//       const char* const FnName = "IsLRSDegenerateOrderMod_NoArgChecks";
//       VerboseLog VERBOSE(FnName);
//       CoCoA_ASSERT(n >= 2);
//       if (n == 2)
//       {
//         const ring &P = owner(f);
//         const RingElem &x = indet(P, IndetIndex);
//         // Equiv to !IsCoprime
//         return (!IsConstant(gcd(f, PolyAlgebraHom(P, P, std::vector<RingElem>(NumIndets(P), -x))(f))));
//       }

//       const long NumPrimesToTry = IsGuaranteed(VerifLev) ? 3 : level(VerifLev); // magic number
//       if (!IsLRSDegenerateOrder_ModularCheck(f, n, NumPrimesToTry))  return false;

// ////      VERBOSE(20) << "Exponent " << n << " passed modular check" << std::endl;
//       return true;
      // if (!IsGuaranteed(VerifLev))  return true;

      // if (!IsPrime(n))
      // {
      //   // SLOW but surely correct -- also works if n is prime, but is slower than GraeffeN
      //   const ring &P = owner(f);
      //   const ring& k = CoeffRing(P);
      //   const ring Kxy = NewPolyRing(k, NewSymbols(2));
      //   const RingElem& x = indet(Kxy,0);  const RingElem& y = indet(Kxy,1);
      //   RingHom x_to_xy = PolyAlgebraHom(P, Kxy, std::vector<RingElem>(NumIndets(P), x*y));
      //   RingHom x_to_y = PolyAlgebraHom(P, Kxy, std::vector<RingElem>(NumIndets(P), y));
      //   const RingElem R = resultant(x_to_xy(f), cyclotomic(n,x), 0); // 0 is IndetIndex of x
      //   return (!IsConstant(gcd(R, x_to_y(f))));
      // }
      // // Might be a quicker way
      // const RingElem Pn = GraeffeN(f, n); // equal to resultant(eval(f, [y]), y^n-x, y);
      // return (!IsSqFree(Pn));
//    }


    enum class FirstOrAll { JustFirst, FindAll };

    /**
     * Modular approach using the given verification level
     * If MULT == 0 then check all possible orders;
     * if MULT != 0 then want to know if f is MULT-LRS-degenerate, so
     * check only orders which divide MULT -- if MULT is trivially excluded, we can return early.
     */
    std::vector<unsigned long> LRSDegeneracyOrders_modular(RingElem f, unsigned long MULT, VerificationLevel VerifLev, FirstOrAll StopCriterion)
    {
      const char* const FnName = "LRSDegeneracyOrders_modular";
      VerboseLog VERBOSE(FnName);

      const ErrorContext context = CoCoA_ERROR_CONTEXT;
      f = LRSDegenerate_check_and_preprocess(f, context);
      vector<unsigned long> ListOfOrders;
      if (deg(f) < 2)  return ListOfOrders;
      vector<RingElem> CommonFac; //  useful only if IsGuaranteed(VerifLev)
      const ring &P = owner(f);
      const ring& K = CoeffRing(P);
      const ring Kxy = NewPolyRing(K, NewSymbols(2));
      const RingElem& x = indet(Kxy,0);  const RingElem& y = indet(Kxy,1);
      RingHom x_to_xy = PolyAlgebraHom(Kxy, Kxy, std::vector<RingElem>{x*y, zero(Kxy)});
      RingHom x_to_y = PolyAlgebraHom(Kxy, Kxy, std::vector<RingElem>{y, zero(Kxy)});
      vector<RingElem> images(NumIndets(P), x);
      const RingHom to_Kxy = PolyAlgebraHom(P, Kxy, images);
      const RingElem fx = to_Kxy(prim(f));
      const RingElem fy = x_to_y(fx);
      const RingElem fxy = x_to_xy(fx);

      vector<RingElem> NegImages{-x,zero(Kxy)};
      RingHom NegX = PolyAlgebraHom(Kxy, Kxy, NegImages);
      if (IsEven(MULT) && !IsConstant(gcd(fx, NegX(fx))))
      {
        ListOfOrders.push_back(2);
        if (MULT == 2 || StopCriterion == FirstOrAll::JustFirst)
          return ListOfOrders;
        CommonFac.push_back(x_to_y(gcd(fx, NegX(fx))));
      }

      const long degf = deg(f);
      if (degf > std::numeric_limits<long>::max()/degf)  // I'd be amazed if this ever happens!
        CoCoA_THROW_ERROR("Deg too big", "LRSDegeneracyOrder_modular");
      /*const*/ std::vector<unsigned long> KToTry_init = GetKToTryList(degf/*, Kbound*/);
      // lines below would be neater with std::erase_if (in C++ 20)
      std::vector<unsigned long> KToTry;
      if (MULT == 0)
      {
        swap(KToTry, KToTry_init);
      }
      else
      {
        for (unsigned long k: KToTry_init)
          if (MULT%k == 0)
            KToTry.push_back(k);
        if (KToTry.empty() || KToTry.back() != MULT)
          return ListOfOrders; // empty or just [2]; either way indicates that f is not MULT-LRS-degenerate
      }
      // Do modular check with order MULT; might be "obviously not" MULT-LRS-degenerate
      if (MULT != 0 && !IsLRSDegenerateOrder_ModularCheck(f, MULT, 3/*VerifLev*/))
        return ListOfOrders; // surely does not contain MULT
      VERBOSE(80) << "How many k to try? " << len(KToTry) << std::endl;
      if (!KToTry.empty())  VERBOSE(80) << "   Largest k is " << KToTry.back() << std::endl;
      for (unsigned long k : KToTry)
      {
//CONJECTURE//if(EulerTotient(k)>degf)continue;
        CheckForInterrupt("LRSDegeneracyOrder_modular");
        if (!IsLRSDegenerateOrder_ModularCheck(f, k, 3/*VerifLev*/))  continue;
        VERBOSE(90) << "Candidate " << k << " passed modular check" << std::endl;
        if (!IsGuaranteed(VerifLev))  { ListOfOrders.push_back(k); continue; }
        RingElem R = resultant(fxy, power(x,k)-1, 0); // 0 is IndetIndex of x
        R /= fy;
        for (int i=0; i < len(CommonFac); ++i)
          if (k % ListOfOrders[i] == 0)
            R /= CommonFac[i];
        R = gcd(R, fy);
        if (!IsConstant(R))
        {
          ListOfOrders.push_back(k);
          CommonFac.push_back(R);
          if (StopCriterion == FirstOrAll::JustFirst)
            return ListOfOrders;
        }
        else { VERBOSE(90) << "False positive for order " << k << std::endl; }
      }
      return ListOfOrders;
    }


  }  // end of namespace anonymous


  /**
   * Checks for LRS-degeneracy of the given order using the given verification level
   */
  bool IsLRSDegenerateOrder(RingElem f, const unsigned long k, VerificationLevel VerifLev)
  {
    const vector<unsigned long> L = LRSDegeneracyOrders_modular(f, k, VerifLev, FirstOrAll::FindAll);
    return (!L.empty() && k == L.back());
  }


  /**
   * Checks for LRS-degeneracy of the given order
   */
  bool IsLRSDegenerateOrder(ConstRefRingElem f, const unsigned long n)
  {
    return IsLRSDegenerateOrder(f, n, guaranteed());
  }




  // Smallest LRS-degeneracy order (or 0 if not LRS-degenerate)
  unsigned long LRSDegeneracyOrder(ConstRefRingElem f)
  {
    const vector<unsigned long> orders = LRSDegeneracyOrders_modular(f, 0ul, guaranteed(), FirstOrAll::JustFirst);
    if (!orders.empty())  return orders[0];
    return 0;
  }

  std::vector<unsigned long> LRSDegeneracyOrders(ConstRefRingElem f)
  {
    return LRSDegeneracyOrders_modular(f, 0ul, guaranteed(), FirstOrAll::FindAll);
  }

  std::vector<unsigned long> LRSDegeneracyOrders(ConstRefRingElem f, VerificationLevel VerLev)
  {
    return LRSDegeneracyOrders_modular(f, 0ul, VerLev, FirstOrAll::FindAll);
  }


  bool IsLRSDegenerate(RingElem f)
  {
    f = LRSDegenerate_check_and_preprocess(f, CoCoA_ERROR_CONTEXT);
    // We do a quick check for cyclo factors before trying general method.
    const auto facs = CyclotomicFactors(f);
    if (!facs.myFactors().empty())
    {
      // Definitely LRS-degenerate unless the only cyclo factor is x-1
      // That's what we're checking for here... terribly long-winded!!
      if (len(facs.myFactors()) > 1 ||
          deg(facs.myFactors()[0]) > 1 ||
          IsOne(ConstantCoeff(facs.myFactors()[0])))
        return true;
    }
    return (LRSDegeneracyOrder(f) != 0);
  }

}
