//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RootBound.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/error.H"
#include "CoCoA/verbose.H"

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// This test checks that RootBound returns a plausble result (CheckUPB)
// and that it is invariant as expected (wrt mult-by-x, mult-by-scalar,
// and subst(x |--> scalar*x) [except result changes by factor of scalar]
//----------------------------------------------------------------------

namespace CoCoA
{

  // Check that A and B are equal to within about 0.5%
  // Since result of RootBound is approx (to within about 0.1%)
  // we cannot check for equality, so we check IsCloseTo instead.
  bool IsCloseTo(const BigRat& A, const BigRat& B)
  {
    if (IsZero(A)) return IsZero(B);
    if (IsZero(B)) return false;
    if (sign(A) != sign(B)) return false;
    
    return ( abs(A-B)/(A+B) < BigRat(1,100) );
  }


  // Given a univariate poly f with rational coeffs, return "f-star".
  RingElem Fstar(const RingElem& f)
  {
    const ring& P = owner(f);
    RingElem ans = monomial(P, abs(LC(f)), LPP(f));
    for (SparsePolyIter it = ++BeginIter(f); !IsEnded(it); ++it)
    {
      ans -= monomial(P, abs(coeff(it)), PP(it));
    }
    return ans;
  }


  void CheckUPB(ConstRefRingElem f)
  {
    const BigRat B = RootBound2(f,0); // *IMPORTANT* no graeffe iterations for this test.
                                      // o/w do not nec have that fstar(B) >= 0
    CoCoA_ASSERT_ALWAYS(B >= 0);
    CoCoA_ASSERT_ALWAYS(IsMonomial(f) || B > 0);

    const PolyRing P(owner(f), CoCoA_ERROR_CONTEXT);
    const ring& k = CoeffRing(P);
    const long index = UnivariateIndetIndex(f);

    const RingElem g = Fstar(f);
    vector<RingElem> images(NumIndets(P), zero(k));
    images[index] = RingElem(k, B);
    RingHom eval = PolyRingHom(P, k, IdentityHom(k), images);
    const RingElem evalg = eval(g);
    CoCoA_ASSERT_ALWAYS(sign(evalg) >= 0);
    const BigRat B1 = RootBound2(f,1);
    CoCoA_ASSERT_ALWAYS(B1 <= B);
    const BigRat B2 = RootBound2(f,2);
    CoCoA_ASSERT_ALWAYS(B2 <= B1);
  }
  

  void CheckInvariance(ConstRefRingElem f)
  {
    const BigRat B = RootBound2(f);

    // Check that bound is invariant if poly is multiplied by a scalar
    for (int k=-10; k <= 10; ++k)
    {
      if (k == 0) continue;
      CoCoA_ASSERT_ALWAYS(RootBound2(k*f) == B);
      CoCoA_ASSERT_ALWAYS(RootBound2(f/k) == B);
    }
    
    // Check that bound is invariant if poly is multiplied by a power of x
    const PolyRing P(owner(f), CoCoA_ERROR_CONTEXT);
    const long index = UnivariateIndetIndex(f);
    const RingElem& x = indet(P,index);
    
    for (int k=1; k < 20; ++k)
      CoCoA_ASSERT_ALWAYS(RootBound2(f*power(x,k)) == B);

    // Check invariance under x |--> k*x, and under x |--> (1/k)*x
    for (int k=-20; k < 20; ++k)
    {
      if (k == 0) continue;
      vector<RingElem> images = indets(P);
      images[index] *= k;
      const RingHom phi = PolyAlgebraHom(P,P, images);
      CoCoA_ASSERT_ALWAYS(IsCloseTo(B, abs(BigInt(k))*RootBound2(phi(f))));

      images[index] /= k*k;
      const RingHom psi = PolyAlgebraHom(P,P, images);
      CoCoA_ASSERT_ALWAYS(IsCloseTo(B, RootBound2(psi(f))/abs(BigInt(k))));
    }
  }


  RingElem RndPoly(long d, const RingElem& x)
  {
    RingElem ans(owner(x));
    for (int i=0; i <= d; ++i)
    {
      ans += RandomLong(-9,9)*power(x,i);
    }
    return ans;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    const PolyRing P = NewPolyRing(RingQQ(), symbols("x,y"));
    const RingElem& x = indet(P,0);

    // Linear polynomials
    CheckUPB(x);
    CheckInvariance(x);

    CheckUPB(x-1);
    CheckInvariance(x-1);

    CheckUPB(x+1);
    CheckInvariance(x+1);

    CheckUPB(21*x-34);
    CheckInvariance(21*x-34);

    CheckUPB(21*x+34);
    CheckInvariance(21*x+34);

    CheckUPB(144*x-89);
    CheckInvariance(144*x-89);

    CheckUPB(144*x+89);
    CheckInvariance(144*x+89);

    // Higher deg polys
    CheckUPB(power(x,2));
    CheckInvariance(power(x,2));

    CheckUPB(power(x,3));
    CheckInvariance(power(x,3));

    CheckUPB(power(x,120));
    CheckInvariance(power(x,120));

    const RingElem f= RingElem(P, "x^2-2");
    CheckUPB(f);
    CheckInvariance(f);

    CheckUPB(power(f,2));
    CheckInvariance(power(f,2));

    CheckUPB(power(f,3));
    CheckInvariance(power(f,3));

    // Check a few random polys
    for (int d=1; d < 10; ++d)
      for (int i=0; i < 8; ++i)
      {
        const RingElem frand = RndPoly(d,x);
        if (IsZero(frand) || deg(frand) == 0) continue;
        CheckUPB(frand);
        CheckInvariance(frand);
      }
  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
