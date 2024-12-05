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

#include "CoCoA/SparsePolyOps-SturmSeq.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H

//#include<vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {
    
    RingElem remainderZZ(RingElem f, const RingElem& g, const RingElem& x)
    {
      const int degg = deg(g);
      while (!IsZero(f) && deg(f) >= degg)
      {
        const BigInt LCf = ConvertTo<BigInt>(LC(f));
        const BigInt LCg = ConvertTo<BigInt>(LC(g));
        const BigInt h = gcd(LCf, LCg);
        if (sign(LCg) > 0)
          f = (LCg/h)*f - (LCf/h) * power(x,deg(f)-degg) * g;
        else
          f = (-LCg/h)*f + (LCf/h) * power(x,deg(f)-degg) * g;
      }
      if (IsZero(f)) return f;
      const RingHom coeff = CoeffEmbeddingHom(owner(f));
      return f/coeff(content(f));
    }

    // Naive slow version
    // RingElem remainder(RingElem f, const RingElem& g, const RingElem& x)
    // {
    //   const int degg = deg(g);
    //   RingHom phi = CoeffEmbeddingHom(owner(f));
    //   while (!IsZero(f) && deg(f) >= degg)
    //   {
    //     f -= phi(LC(f)/LC(g))*power(x,deg(f)-degg) * g;
    //   }
    //   return f;
    // }

  } // end of namespace anonymous

  
  // Naive Impl: simple rather than fast
  // Taken from Wikipedia page for "Sturm's Theorem"
  // We multiply entries so that coeffs are integer.
  // ASSUMES coeffs are in QQ
  std::vector<RingElem> SturmSeq(const RingElem& f)
  {
    const PolyRing P(owner(f), CoCoA_ERROR_CONTEXT);
    const char* const FnName = "SturmSeq";
////    if (!IsPolyRing(P)) CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName);
    if (!IsQQ(CoeffRing(P))) CoCoA_THROW_ERROR(ERR::BadRing, FnName);
    if (!IsField(CoeffRing(P))) CoCoA_THROW_ERROR(ERR::ReqField, FnName);
    if (!IsOrderedDomain(CoeffRing(P))) CoCoA_THROW_ERROR(ERR::ReqOrdDom, FnName);
    if (IsZero(f)) CoCoA_THROW_ERROR(ERR::ReqNonZero, FnName);
    const int var = UnivariateIndetIndex(f);
    if (var < 0) CoCoA_THROW_ERROR(ERR::ReqUnivariate, FnName);
    const RingElem& x = indet(P,var);
    vector<RingElem> seq;  seq.reserve(1+deg(f));
    RingElem prev = f*ConvertTo<BigInt>(CommonDenom(f));
    // ??? SLUG:  remove content from prev???
    seq.push_back(prev);
    RingElem curr = deriv(prev,x);
    { const BigInt c = ConvertTo<BigInt>(content(curr)); if (!IsOne(c)) curr /= c; }    
    while (!IsZero(curr))
    {
      seq.push_back(curr);
      RingElem next = -remainderZZ(prev, curr, x);
      prev = curr;
      curr = next;
    }
    return seq;
  }


  long NumRealRoots(const RingElem& f)
  {
    const ring& P = owner(f);
    const char* const FnName = "NumRealRoots";
    if (!IsPolyRing(P)) CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName);
    if (!IsQQ(CoeffRing(P))) CoCoA_THROW_ERROR(ERR::BadRing, FnName);
    if (!IsField(CoeffRing(P))) CoCoA_THROW_ERROR(ERR::ReqField, FnName);
    if (!IsOrderedDomain(CoeffRing(P))) CoCoA_THROW_ERROR(ERR::ReqOrdDom, FnName);
    if (IsZero(f)) CoCoA_THROW_ERROR(ERR::ReqNonZero, FnName);
    if (deg(f) == 0) CoCoA_THROW_ERROR("Poly must have deg >= 1", FnName);
    const int var = UnivariateIndetIndex(f);
    if (var < 0) CoCoA_THROW_ERROR(ERR::ReqUnivariate, FnName);
    const vector<RingElem> Sseq = SturmSeq(f);
    long CountPlusInf = 0;
    long CountMinusInf = 0;
    int SignPlusInf = 0;   // 3-way value!!  -1, 0 or +1
    int SignMinusInf = 0;  // 3-way value!!  -1, 0 or +1
    const int n = len(Sseq);
    for (int i=0; i < n; ++i)
    {
      const int SignLC = sign(LC(Sseq[i]));
//      if (SignLC == 0) continue;
      if (SignLC != SignPlusInf) { SignPlusInf = SignLC; ++CountPlusInf; }
      if (IsEven(deg(Sseq[i])))
      { if (SignLC != SignMinusInf) { SignMinusInf = SignLC; ++CountMinusInf; } }
      else { if (-SignLC != SignMinusInf) { SignMinusInf = -SignLC; ++CountMinusInf; } }
    }
    return CountMinusInf - CountPlusInf;
  }



} // end of namespace CoCoA
