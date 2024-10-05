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
//#include "CoCoA/SparsePolyRing.H"
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
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

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


  void CheckUPB(ConstRefRingElem f)
  {
    const BigRat B = RootBound2(f);
    CoCoA_ASSERT_ALWAYS(B >= 0);
    CoCoA_ASSERT_ALWAYS(IsMonomial(f) || B > 0);

    const PolyRing P(owner(f), CoCoA_ERROR_CONTEXT);
    const ring& k = CoeffRing(P);
    const long index = UnivariateIndetIndex(f);

    const RingElem radf = radical(f);
    RingElem g = monomial(P, LC(radf), LPP(radf));
    for (SparsePolyIter it = ++BeginIter(radf); !IsEnded(it); ++it)
    {
      g -= monomial(P,abs(coeff(it)), PP(it));
    }
    vector<RingElem> images(NumIndets(P), zero(k));
    images[index] = RingElem(k, B);
    RingHom eval = PolyRingHom(P, k, IdentityHom(k), images);
    const RingElem evalg = eval(g);
    CoCoA_ASSERT_ALWAYS(sign(evalg) >= 0);
  }
  

  void CheckInvariance(ConstRefRingElem f)
  {
    const BigRat B = RootBound2(f);

    // Check that bound is invariant if poly is multiplied by a scalar
    for (int k=-20; k < 20; ++k)
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


  void program()
  {
    GlobalManager CoCoAFoundations;

    ring P = NewPolyRing(RingQQ(), symbols("x,y"));

    // Check "correctness" of each bound at certain given polys
    const RingElem f1 = RingElem(P, "x^2+2*x+2"); // Cauchy is best
    const RingElem f2 = RingElem(P, "2*x^2+x+1"); // Lagrange is best
    const RingElem f3 = RingElem(P, "x^2+2*x+3"); // Birkhoff is best
    const RingElem f4 = RingElem(P, "x^2+3*x+2"); // LMS is best

    CoCoA_ASSERT_ALWAYS(RootBound_simple(f1) == RootBound_Cauchy(f1));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Cauchy(f1), BigRat(3,1)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Lagrange(f1), BigRat(4,1)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Birkhoff(f1), BigRat(219,64)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_LMS(f1), BigRat(219,64)));
    
    CoCoA_ASSERT_ALWAYS(RootBound_simple(f2) == RootBound_Lagrange(f2));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Cauchy(f2), BigRat(3,2)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Lagrange(f2), BigRat(1)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Birkhoff(f2), BigRat(219,128)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_LMS(f2), BigRat(155,128)));

    CoCoA_ASSERT_ALWAYS(RootBound_simple(f3) == RootBound_LMS(f3));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Cauchy(f3), BigRat(4,1)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Lagrange(f3), BigRat(5,1)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Birkhoff(f3), BigRat(67,16)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_LMS(f3), BigRat(479,128)));

    CoCoA_ASSERT_ALWAYS(RootBound_simple(f4) == RootBound_Birkhoff(f4));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Cauchy(f4), BigRat(4,1)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Lagrange(f4), BigRat(5,1)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_Birkhoff(f4), BigRat(29,8)));
    CoCoA_ASSERT_ALWAYS(IsCloseTo(RootBound_LMS(f4), BigRat(283,64)));
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
