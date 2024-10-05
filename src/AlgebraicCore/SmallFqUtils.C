//   Copyright (c)  2015  John Abbott and Anna M. Bigatti

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

#include "CoCoA/SmallFqUtils.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SmallFqVecImpl.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/convert.H"
#include "CoCoA/ring.H"

//#include <vector>
using std::vector;

namespace CoCoA
{

  std::vector<SmallFpImpl::value> PolyToVec(const RingElem& m, const SmallFpImpl& FFp)
  {
    const int d = deg(m);
    vector<SmallFpImpl::value> coeffs(1+d, zero(SmallFp)/*FFp.myReduce(0)*/);
    const vector<RingElem> C = CoeffVecWRT(m, indet(owner(m),0));
    for (int i=0; i <= d; ++i)
      coeffs[i] = FFp.myReduce(ConvertTo<long>(C[i]));
    return coeffs;
  }

  RingElem VectorToPoly(const ring& Fpx, const std::vector<SmallFpImpl::value>& v)
  {
    const long p = ConvertTo<long>(characteristic(Fpx)); // BUG BUG BUG   cheating!!!
    SmallFpImpl ModP(p);
    const long d = len(v)-1;
    const RingElem x = indet(Fpx,0);
    RingElem f(Fpx);
    for (int i=d; i >= 0; --i)
    {
      f *= x;
      f += ModP.myExport(v[i]); // inefficient because it forces reduction mod p of result of myExport!
    }
    return f;
  }

  void next(const SmallFpImpl& ModP, vector<SmallFpImpl::value>& v)
  {
    int i=0;
    while (i < len(v))
    {
      v[i] = ModP.myAdd(v[i],one(SmallFp));
      if (!IsZero(v[i])) return;
      ++i;
    }
  }

  bool IsOne(const vector<SmallFpImpl::value>& x) noexcept
  {
    if (!IsOne(x[0])) return false;
    const int n = len(x);
    for (int i=1; i < n; ++i)
      if (!IsZero(x[i])) return false;
    return true;
  }

  long HasOrder(long n, const FFqImpl_vec& extn, const std::vector<SmallFpImpl::value>& x)
  {
    // const long p = ModP.myModulus();

    // long q = SmallPower(p,d);
//    factorization<long> qfacs = factor_TrialDiv(n,n);
    const int d = extn.myExtnDeg();
    const vector<long> primes = factor_TrialDiv(n,n).myFactors();
///    clog<<"HasOrder: n="<<n<<"   primes="<<primes<<endl;
//???      vector<SmallFpImpl::value> X(d); X[1]=1;
vector<SmallFpImpl::value> pwr(d); 
    for (int i=0; i < len(primes); ++i)
    {
      const long e = n/primes[i];
//      clog<<"Checking pwr e="<<e<<"   x="<<X<<endl;
      extn.myPower(&pwr[0],&x[0],e);
      if (IsOne(pwr)) return false;
    }
    return true;
}


  bool IsGroupGen(const SmallFpImpl& ModP, const std::vector<SmallFpImpl::value>& MinPoly)
  {
    const long d = len(MinPoly)-1;
    CoCoA_ASSERT(IsOne(MinPoly[d]));
    if (d > 1 && IsZero(MinPoly[0])) return false;
    const long p = ModP.myModulus();
    ring Fp = NewRingFp(p);
    ring Fpx = NewPolyRing(Fp, NewSymbols(1));  // wasteful to create new ring each time
    const RingElem x = indet(Fpx,0);
    RingElem f = VectorToPoly(Fpx, MinPoly);
    if (!IsIrred(f)) return false;
///    clog<<"IRRED OK"<<endl;

//??    FFqImpl_vec extn(ModP,MinPoly);
    FFqImpl_vec extn(p,MinPoly);
    vector<SmallFpImpl::value> X(d, zero(SmallFp)); X[1]=one(SmallFp);
    return HasOrder(SmallPower(p,d)-1,extn,X);
}

    std::vector<SmallFpImpl::value> FindGroupGenerator(long p, int d)
{
  SmallFpImpl ModP(p);
  typedef SmallFpImpl::value FpElem;
  vector<FpElem> f(1+d);
  f[d] = one(SmallFp);
  f[1] = one(SmallFp); // since gen cannot be of form x^d + const

//  std::clog << "Looking for group gen..." << endl;
  while (!IsGroupGen(ModP, f))
    next(ModP,f);

//  std::clog << "Group gen is " << f << endl;
  return f;
}


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SmallFqUtils.C,v 1.12 2022/02/18 14:11:58 abbott Exp $
// $Log: SmallFqUtils.C,v $
// Revision 1.12  2022/02/18 14:11:58  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.11  2022/02/02 09:25:15  abbott
// Summary: Renamed SmoothFactor to factor_TrialDiv (redmine 950)
//
// Revision 1.10  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.9  2021/01/07 15:16:53  abbott
// Summary: Corrected copyright
//
// Revision 1.8  2019/03/18 11:20:00  abbott
// Summary: Changed include (after changing NumTheory)
//
// Revision 1.7  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.6  2018/05/18 16:42:11  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.5  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.4  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.3  2018/02/27 10:55:52  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.2  2016/11/03 12:28:16  abbott
// Summary: Added comment
//
// Revision 1.1  2015/12/18 15:25:07  abbott
// Summary: Added impls of non-prime finite fields
//
//
