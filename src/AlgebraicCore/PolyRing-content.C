//   Copyright (c)  2021  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/PolyRing.H"

#include "CoCoA/BigInt.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/VectorOps.H"    // for HasUniqueOwner
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"  // for len


#include <vector>
// using std::vector;


namespace CoCoA
{

  // Functions dealing with scaling polynomials by a constant:
  // content, monic, CommonDenom, ClearDenom, prim

  RingElem content(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "content(f)");
    const PolyRing Rx = owner(f);
    const ring R = CoeffRing(Rx);

    if (IsFractionFieldOfGCDDomain(R))
    {
      if (IsZero(f)) return zero(R);
      RingElem ans(R);
      Rx->myContentFrF(raw(ans), raw(f));
      return ans;
    }
    if (IsTrueGCDDomain(R))
    {
      if (IsZero(f)) return zero(R);
      RingElem ans(R);
      Rx->myContent(raw(ans), raw(f));
      return ans;
    }
    CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "content(f)");
    return zero(R); // never reached, just to keep compiler quiet
  }


  RingElem CommonDenom(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "CommonDenom(f)");
    const PolyRing Qx = owner(f);
    const ring Q = CoeffRing(Qx);
    if (!IsFractionField(Q))
      CoCoA_THROW_ERROR(ERR::NotElemFrF, "CoeffRing inside CommonDenom(f)");
    const ring R = BaseRing(Q);
    if (!IsTrueGCDDomain(R))
      CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "BaseRing of CoeffRing inside CommonDenom(f)");
//     if (IsField(R))  // see documentation (Bugs section)
//       CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "content(f)");
    if (IsZero(f)) return one(R);

    RingElem ans(R);
    Qx->myCommonDenom(raw(ans), raw(f));
    return ans;
  }


  RingElem CommonDenom(const std::vector<RingElem>& v)
  {
    if (v.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(v))
      CoCoA_THROW_ERROR(ERR::MixedRings, "CommonDenom(vector<RingElem> v)");
    RingElem q = CommonDenom(v[0]);
    for (long i=1; i<len(v); ++i)  q = lcm(q, CommonDenom(v[i]));
    return q;
  }


  RingElem ClearDenom(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "ClearDenom(f)");
    const PolyRing Qx = owner(f);
    const ring Q = CoeffRing(Qx);
    if (!IsFractionField(Q))
      CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "CoeffRing inside ClearDenom(f)");
    const ring R = BaseRing(Q);
    if (!IsTrueGCDDomain(R))
      CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "BaseRing of CoeffRing inside ClearDenom(f)");
//     if (IsField(R))  // see documentation (Bugs section)
//       CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "content(f)");
    if (IsZero(f)) return zero(Qx);

    RingElem ans(Qx);
    Qx->myClearDenom(raw(ans), raw(f));
    return ans;
  }


  RingElem monic(ConstRefRingElem f)
  {
    const PolyRing Rx = owner(f);
    RingElem ans(Rx);
    Rx->myMonic(raw(ans), raw(f));
    return ans;
  }


  namespace // anonymous
  {

    void prim_NumDen(BigInt& N, BigInt& D, ConstRefRingElem f)
    {
      const ring& P = owner(f);
      if (IsZZ(P)) { N = ConvertTo<BigInt>(f); D = 1; return; }
      if (IsQQ(P)) { const BigRat q = ConvertTo<BigRat>(f); N = num(q); D = den(q); return; }
//      if (IsFractionField(P)) ???
      if (!IsPolyRing(P) || !IsZero(characteristic(P)))
        CoCoA_THROW_ERROR("Expected elem of poly ring of char 0", "prim_NumDen");
      N = 0; D = 1;
      BigInt n,d; // used in loops below
      if (IsDenseUPolyRing(P))
      {
        const long degf = deg(f);
        for (long k=0; k <= degf; ++k)
        {
          RingElemAlias c = coeff(f,k);
          if (IsZero(c)) continue;
          prim_NumDen(n, d, c);
          N = gcd(N, n);
          D = lcm(D, d);
        }
      }
      else // Sparse poly repr
      {
        for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
        {
          prim_NumDen(n, d, coeff(it));
          N = gcd(N, n);
          D = lcm(D, d);
        }
      }
    }
    
  } // end of anonymous namespace


  // Essentially same as g = monic(f); return g*CommonDenom(g);
  // Chooses sign so that LC(ans) > 0
  RingElem prim(ConstRefRingElem f)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P) || characteristic(P) != 0)
      CoCoA_THROW_ERROR("Expected elem of poly ring of char 0", "prim");
    if (IsZero(f)) return f;
    BigInt N, D;
    prim_NumDen(N,D, f);
    if (LC(f) < 0) { D = -D; }
    // return (D*f)/N;  // lines below avoid wasteful mult/div by 1
    if (IsOne(N)) { if (IsOne(D)) return f; else return D*f; }
    if (IsOne(D)) return f/N; else return (D*f)/N;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PolyRing-content.C,v 1.5 2024/07/15 17:10:00 bigatti Exp $
// $Log: PolyRing-content.C,v $
// Revision 1.5  2024/07/15 17:10:00  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.4  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.3  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2021/07/19 11:14:00  abbott
// Summary: Minor: improved err mesg, removed cruft
//
// Revision 1.1  2021/03/01 19:39:43  abbott
// Summary: Split off from PolyRing.C; added prim (redmine 1580)
//
//
