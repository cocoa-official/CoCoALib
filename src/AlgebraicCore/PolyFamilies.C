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

#include "CoCoA/PolyFamilies.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-resultant.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"

#include <limits>
using std::numeric_limits;
#include <utility>
using std::pair;
using std::make_pair;
#include <vector>
using std::vector;


namespace CoCoA
{

  //-------------------------------------------------------
  // Chebyshev polys of 1st type
  //-------------------------------------------------------

  // Explicit formula from Wikipedia
  RingElem ChebyshevPoly_explicit(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    BigInt c = power(2, n-1);
    RingElem ans = c*power(x,n);
    for (int k=1; k <= n/2; ++k) // integer division!
    {
      c = -(c*(n-2*k+2)*(n-2*k+1))/(4*k*(n-k));
      ans += c*power(x,n-2*k);
    }
    return ans;
  }


  RingElem ChebyshevPoly_iter(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    if (n == 0) return one(owner(x));
    if (n == 1) return x;
    RingElem prev = one(owner(x));
    RingElem curr = x;
    const RingElem twox = 2*x;
    for (int k = 2; k <= n; ++k)
    {
      RingElem tmp = twox*curr - prev;
      swap(prev, curr);
      swap(curr, tmp);
    }
    return curr;
  }


  std::pair<RingElem, RingElem> ChebyshevPoly_recursive(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    if (n == 0) return std::make_pair(one(owner(x)), x);
    if (n == 1) return std::make_pair(x, 2*x*x-1);
    const long halfn = n/2; // integer division
    const std::pair<RingElem,RingElem> RecAns = ChebyshevPoly_recursive(halfn, x);
    const RingElem& Ta = RecAns.first;
    const RingElem& Tb = RecAns.second;
    if (IsEven(n))
      return make_pair(2*power(Ta,2)-1, 2*Ta*Tb - x);
    return make_pair(2*Ta*Tb-x, 2*power(Tb,2)-1);
  }


  // Dispatch function
  RingElem ChebyshevPoly(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "ChebyshevPoly: 1st arg");
    const ring& P = owner(x);
    if (IsPolyRing(P) && IsIndet(x))
      return ChebyshevPoly_explicit(n, x);
    return ChebyshevPoly_recursive(n, x).first;
  }


  //-------------------------------------------------------
  // Chebyshev polys of 2nd type
  //-------------------------------------------------------

  // Explicit formula from Wikipedia
  RingElem ChebyshevPoly2_explicit(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "ChebyshevPoly: 1st arg");
    BigInt c = power(2, n);
    RingElem ans = c*power(x,n);
    for (int k=1; k <= n/2; ++k) // integer division!
    {
      c = -(c*(n-2*k+2)*(n-2*k+1))/(4*k*(n-k+1));
      ans += c*power(x,n-2*k);
    }
    return ans;
  }

  RingElem ChebyshevPoly2_iter(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "ChebyshevPoly2: 1st arg");
    if (n == 0) return one(owner(x));
    if (n == 1) return x;
    RingElem prev = one(owner(x));
    RingElem curr = 2*x;
    const RingElem twox = 2*x;
    for (int k = 2; k <= n; ++k)
    {
      RingElem tmp = twox*curr - prev;
      swap(prev, curr);
      swap(curr, tmp);
    }
    return curr;
  }

  // Dispatch function
  RingElem ChebyshevPoly2(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "ChebyshevPoly2: 1st arg");
    const ring& P = owner(x);
    if (IsPolyRing(P) && IsIndet(x))
      return ChebyshevPoly2_explicit(n, x);
    return ChebyshevPoly2_iter(n, x);
  }


  //=======================================================

  //-------------------------------------------------------
  // Hermite polys  (physics)
  //-------------------------------------------------------

  // Explicit formula from Wikipedia
  RingElem HermitePoly_explicit(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    BigInt c = power(2, n);
    RingElem ans = c*power(x,n);
    for (int k=1; k <= n/2; ++k) // integer division!
    {
      c = -(c*(n-2*k+2)*(n-2*k+1))/(4*k);
      ans += c*power(x,n-2*k);
    }
    return ans;
  }


  // Physicists' Hermite polys:
  RingElem HermitePoly_iter(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    if (n == 0) return one(owner(x));
    const RingElem twox = 2*x;
    if (n == 1) return twox;
    RingElem prev = twox;
    for (int k = 2; k <= n; ++k)
    {
      RingElem curr = twox*prev - deriv(prev,x);
      swap(curr, prev); // really just prev = curr;
    }
    return prev;
  }


  // Dispatch function
  RingElem HermitePoly(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "HermitePoly: 1st arg");
    const ring& P = owner(x);
    if (IsPolyRing(P) && IsIndet(x))
      return HermitePoly_explicit(n, x);
    return HermitePoly_iter(n, x);
  }


  //-------------------------------------------------------
  // Hermite polys  (probability)
  //-------------------------------------------------------

  // Explicit formula from Wikipedia
  RingElem HermitePoly2_explicit(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    BigInt c(1);
    RingElem ans = c*power(x,n);
    for (int k=1; k <= n/2; ++k) // integer division!
    {
      c = -(c*(n-2*k+2)*(n-2*k+1))/(2*k);
      ans += c*power(x,n-2*k);
    }
    return ans;
  }


  RingElem HermitePoly2_iter(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    if (n == 0) return one(owner(x));
    if (n == 1) return x;
    RingElem prev = x;
    for (int k = 2; k <= n; ++k)
    {
      RingElem curr = x*prev - deriv(prev,x);
      swap(curr, prev); // really just assignment prev = curr
    }
    return prev;
  }
  
  // Dispatch function
  RingElem HermitePoly2(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "HermitePoly2: : 1st arg");
    const ring& P = owner(x);
    if (IsPolyRing(P) && IsIndet(x))
      return HermitePoly2_explicit(n, x);
    return HermitePoly2_iter(n, x);
  }

  //=======================================================

  //-------------------------------------------------------
  // Laguerre polys  (mult by factorial so coeffs are integer)
  //-------------------------------------------------------

  // Laguerre poly ***MULT BY FACTORIAL(N)***
  // Inefficient (roughly cubic in n)
  RingElem LaguerrePoly(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "LaguerrePoly: 1st arg");
    if (n == 0) return one(owner(x));
    if (n == 1) return 1-x;
    RingElem prev2 = one(owner(x));
    RingElem prev1 = 1-x;
    for (int k = 2; k <= n; ++k)
    {
      RingElem curr = (2*k-1-x)*prev1 - power(k-1,2)*prev2;
      swap(prev1, prev2); // really just prev2 = prev1;
      swap(curr, prev1);  // really just prev1 = curr;
    }
    return prev1;
  }


  //=======================================================

  //-------------------------------------------------------
  // Bessel polys
  //-------------------------------------------------------

  // Definition copied from Wikipedia
  RingElem BesselPoly(long n, ConstRefRingElem x)
  {
    const ring& P = owner(x);
    if (!IsPolyRing(P)) CoCoA_THROW_ERROR(ERR::ReqPolyRing, "BesselPoly");
    if (n < 0) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "BesselPoly");
    if (n > 32767 && n >= numeric_limits<long>::max()/n) CoCoA_THROW_ERROR(ERR::ArgTooBig, "BesselPoly");
    BigInt coeff(1);
    RingElem ans = one(P);
    RingElem xpow = one(P);
    for (long k = 1; k <= n; ++k)
    {
      xpow *= x;
      coeff = ((((n+k)*(n-k+1))/2)*coeff)/k;  // cannot overflow
      ans += coeff*xpow;
    }
    return ans;
  }

  //=======================================================
  // Dickson polynomials -- defn taken from Wikipedia
  // Impls are simple rather than efficient.

  RingElem DicksonPoly(ConstRefRingElem x, long n, ConstRefRingElem alpha)
  {
    if (n < 0) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "DicksonPoly: 2nd arg");
    if (!IsPolyRing(owner(x))) CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "DicksonPoly");
    const PolyRing& P = owner(x);
    if (owner(alpha) != P) CoCoA_THROW_ERROR(ERR::MixedRings, "DicksonPoly");
    if (n == 0) return one(owner(x));
    const RingElem xsquared = x*x;
    RingElem pwrx = power(x,n);
    RingElem ans = pwrx;
    if (n == 1) return ans;
    RingElem coeff(P,n);
    RingElem pwralpha = -alpha;
    pwrx /= xsquared;
    ans += (coeff*pwrx)*pwralpha;
    for (int i=2; i <= n/2; ++i)
    {
      pwrx /= xsquared;
      pwralpha = (-alpha)*pwralpha;
      coeff = (coeff*(n-2*i+1)*(n-2*i+2))/(i*(n-i)); // exact division!
      ans += (coeff*pwrx)*pwralpha;
    }
    return ans;
  }


  RingElem DicksonPoly2(ConstRefRingElem x, long n, ConstRefRingElem alpha)
  {
    if (n < 0) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "DicksonPoly2: 2nd arg");
    if (!IsPolyRing(owner(x))) CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "DicksonPoly2: 1st arg");
    const PolyRing& P = owner(x);
    if (owner(alpha) != P) CoCoA_THROW_ERROR(ERR::MixedRings, "DicksonPoly2");
    if (n == 0) return one(owner(x));
    const RingElem xsquared = x*x;
    RingElem pwrx = power(x,n);
    RingElem ans = pwrx;
    if (n == 1) return ans;
    RingElem coeff(P,n-1);
    RingElem pwralpha = -alpha;
    pwrx /= xsquared;
    ans += (coeff*pwrx)*pwralpha;
    for (int i=2; i <= n/2; ++i)
    {
      pwrx /= xsquared;
      pwralpha = (-alpha)*pwralpha;
      coeff = (coeff*(n-2*i+1)*(n-2*i+2))/(i*(n-i+1)); // exact division!
      ans += (coeff*pwrx)*pwralpha;
    }
    return ans;
  }


  //=======================================================
  // Swinnerton-Dyer polynomials

  // Just let the BigInt version do the work
  RingElem SwinnertonDyerPoly(ConstRefRingElem x, const std::vector<long>& L)
  {
    const long n = len(L);
    std::vector<BigRat> BigL(n);
    for (int i=0; i < n; ++i)
      BigL[i] =L[i];
    return SwinnertonDyerPoly(x, BigL);
  }

  // This is a relatively simple/naive impl
  RingElem SwinnertonDyerPoly(ConstRefRingElem x, const std::vector<BigRat>& L)
  {
    if (!IsPolyRing(owner(x))) CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "SwinnertonDyerPoly");
    if (!IsIndet(x))  CoCoA_THROW_ERROR(ERR::ReqIndet, "SwinnertonDyerPoly");
    if (L.empty()) return one(owner(x));
    const ring P = NewPolyRing(CoeffRing(owner(x)), symbols("x,y"));
    const RingElem y = indet(P,1);
    RingElem ans = indet(P,0);
    const int n = len(L);
    std::vector<RingElem> XplusYimg; XplusYimg.push_back(ans+y); XplusYimg.push_back(zero(P));
    const RingHom XplusY = PolyAlgebraHom(P,P, XplusYimg);
    for (int i=0; i < n; ++i)
    {
      ans = XplusY(ans);
      ans = resultant(ans, power(y,2)-L[i], 1);
    }
    std::vector<RingElem> images; images.push_back(x); images.push_back(zero(owner(x)));
    const RingHom phi = PolyAlgebraHom(P, owner(x), images);
    return phi(ans);
  }


} // end of namespace CoCoA
