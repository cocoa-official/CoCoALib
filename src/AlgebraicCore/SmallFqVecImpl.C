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

#include "CoCoA/SmallFqVecImpl.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/SmallFqUtils.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/convert.H"
#include "CoCoA/utils.H"

//#include<vector>
using std::vector;

namespace CoCoA
{


  FFqImpl_vec::FFqImpl_vec(long p /*const SmallFpImpl& FFp*/, const vector<FpElem>& m):
      myFFp(p),
      myMinPoly(m),
      myDeg(len(m)-1),
      myTailDeg(myDeg-1) // updated inside body
  {
    while (myTailDeg >= 0 && IsZero(m[myTailDeg]))
      --myTailDeg;
    CoCoA_ASSERT(myTailDeg > 0);
    MakeLogExpTbls();
}


  long FFqImpl_vec::compress(const std::vector<FpElem>& v) const noexcept
  {
//    std::clog<<"compress: arg="<<v<<std::endl;
    long ans = 0;
    for (int i=myDeg-1; i >= 0; --i)
      ans = ans*myModulus() + myFFp.myExportNonNeg(v[i]);
//    std::clog<<"compress: ans="<<ans<<std::endl;
    return ans;
  }
  
  void FFqImpl_vec::MakeLogExpTbls() const
  {
    const long card = SmallPower(myModulus(), myDeg);
//    std::clog<<"card="<<card<<std::endl;
//    std::clog<<"minpoly="<<myMinPoly<<std::endl;
    vector<long> ExpTbl(card);
    vector<long> LogTbl(card);
    vector<long> Add1Tbl(card);
    vector<FpElem> pwr(myDeg);
    pwr[0] = one(SmallFp);
//    std::clog<<"MakeLogExpTbls: pwr="<<pwr<<std::endl;
    for (long exp =0; exp < card; ++exp)
    {
      const long CompressedVal = compress(pwr);
      LogTbl[CompressedVal] = exp;
      ExpTbl[exp] = CompressedVal;
//    std::clog<<"MakeLogExpTbls: pwr="<<pwr<<std::endl;
      myMulByGen(pwr);
//    std::clog<<"MakeLogExpTbls: pwr="<<pwr<<std::endl;
    }
//    std::clog<<"LogTbl="<<LogTbl<<std::endl;
//    std::clog<<"ExpTbl="<<ExpTbl<<std::endl;

    for (int i=0; i < card-1; ++i)
    {
      const long val = ExpTbl[i];
      if (val == myModulus()-1) continue; // skip -1
      if (val%myModulus() != myModulus()-1)
      {
//        std::clog<<"val+1="<<val+1<<std::endl;
        Add1Tbl[i] = LogTbl[val+1];
      }
      else
      {
//        std::clog<<"val+1-p="<<val+1-myModulus<<std::endl;
        Add1Tbl[i] = LogTbl[val+1-myModulus()];
      }
    }
//    std::clog<<"Add1Tbl="<<Add1Tbl<<std::endl;
  }

  void FFqImpl_vec::myMulByGen(std::vector<FpElem>& v) const noexcept
  {
    const FpElem c = myFFp.myNegate(v[myDeg-1]);
//    std::clog<<"myMulByGen: v="<<v<<std::endl;
//    std::clog<<"myMulByGen: c="<<c<<std::endl;
    for (int i=myDeg-1; i > 0; --i)
    {
      v[i] = myFFp.myAddMul(v[i-1], c, myMinPoly[i]);
    }
    v[0] = myFFp.myMul(c, myMinPoly[0]);
//    std::clog<<"myMulByGen: v="<<v<<std::endl;
  }


  FFqImpl_vec::FFqImpl_vec(const ideal& I):
      myFFp(ConvertTo<long>(characteristic(RingOf(I)))),
    myMinPoly(PolyToVec(gens(I)[0], myFFp)),
      myDeg(deg(gens(I)[0]))
  {
    myTailDeg = myDeg-1;
    while (myTailDeg >= 0 && IsZero(myMinPoly[myTailDeg])) --myTailDeg;
    CoCoA_ASSERT(myTailDeg > 0);
  }


  std::ostream& operator<<(std::ostream& out, const FFqImpl_vec& arith);
//??  bool operator==(const FFqImpl_vec& arith1, const FFqImpl_vec& arith2);
//??  bool operator!=(const FFqImpl_vec& arith1, const FFqImpl_vec& arith2);


  void FFqImpl_vec::myGen(OUTvalue_t ans) const noexcept
  {
    for (int i=0; i < myDeg; ++i)
      ans[i] = zero(SmallFp);
    ans[1] = one(SmallFp);
  }

  void FFqImpl_vec::myAdd(OUTvalue_t ans, value_t a, value_t b) const noexcept
  {
    // CoCoA_ASSERT(len(ans) == myDeg-1);
    // CoCoA_ASSERT(len(a) == myDeg-1);
    // CoCoA_ASSERT(len(b) == myDeg-1);

    for (int i=0; i < myDeg; ++i)
      ans[i] = myFFp.myAdd(a[i],b[i]);
  }

  void FFqImpl_vec::mySub(OUTvalue_t ans, value_t a, value_t b) const noexcept
  {
    // CoCoA_ASSERT(len(ans) == myDeg-1);
    // CoCoA_ASSERT(len(a) == myDeg-1);
    // CoCoA_ASSERT(len(b) == myDeg-1);

    for (int i=0; i < myDeg; ++i)
      ans[i] = myFFp.mySub(a[i],b[i]);
  }


  void FFqImpl_vec::myMul(OUTvalue_t ans, value_t a, value_t b) const
  {
    // CoCoA_ASSERT(len(ans) == myDeg-1);
    // CoCoA_ASSERT(len(a) == myDeg-1);
    // CoCoA_ASSERT(len(b) == myDeg-1);
    vector<SmallFpImpl::NonRedValue> tmp(2*myDeg-1); // temp workspace
  for (int i=0; i < myDeg; ++i)
    if (!IsZero(a[i]))
    for (int j=0; j < myDeg; ++j)
      tmp[i+j] += a[i]*b[j];// BUG!  myFFp.myAdd(..., myFFp.myMul(...))

  for (int i=2*myDeg-2; i >= myDeg; --i)
  {
    FpElem c = myFFp.myNormalize(tmp[i]);
    if (IsZero(c)) continue;
    c = myFFp.myNegate(c);

    // Add just the tail
    for (int j=0; j <= myTailDeg; ++j)
      tmp[i-myDeg+j] += c*myMinPoly[j];
  }
  for (int i=0; i < myDeg; ++i)
    ans[i] = myFFp.myNormalize(tmp[i]);
}


  void FFqImpl_vec::myDiv(OUTvalue_t ans, value_t a, value_t b) const
{
//  const int D = len(MinPoly);
//  const ring& Fp = owner(a[0]);
  vector<FpElem> A = myMinPoly;
  vector<FpElem> B(myDeg); for (int i=0; i<myDeg; ++i) { B[i] = b[i]; }
  vector<FpElem> Acofac; Acofac.reserve(myDeg);
  vector<FpElem> Bcofac; Bcofac.reserve(myDeg);
  Bcofac.push_back(one(SmallFp));

  int DegA = myDeg;
  int DegB = myDeg-1; while (DegB >= 0 && IsZero(B[DegB])) --DegB;
  while (DegB >= 0)
  {
    // clog << "OUTER LOOP" << endl;
    // clog << "A = " << A << endl;
    // clog << "B = " << B << endl;
    // clog << "Acofac = " << Acofac << endl;
    // clog << "Bcofac = " << Bcofac << endl;
    while (DegA >= DegB)
    {
//      cout << "INNER LOOP:" << endl;
      const int shift = DegA-DegB;
//      cout << "shift = " << shift << endl;
      const FpElem c = myFFp.myNegate(myFFp.myDiv(A[DegA],B[DegB]));
//      cout << "c = " << c << endl;
      for (int i=0; i < DegB; ++i)
        A[i+shift] = myFFp.myAdd(A[i+shift], myFFp.myMul(c,B[i]));
//        A[i+shift] += c*B[i]; // NORMALIZE!
      if (len(Acofac) < len(Bcofac)+shift) Acofac.resize(len(Bcofac)+shift);
//      clog << "extnded Acofac = " << Acofac << endl;
      for (int j=0; j < len(Bcofac); ++j)
        Acofac[j+shift] = myFFp.myAdd(Acofac[j+shift], myFFp.myMul(c,Bcofac[j]));
//        Acofac[j+shift] += c*Bcofac[j];  // NORMALIZE!
//      clog << "updated Acofac = " << Acofac << endl;
      int NewDegA = DegA-1;
      while (NewDegA >= 0 && IsZero(A[NewDegA])) --NewDegA;
      A.resize(1+NewDegA);
      DegA = NewDegA;
    }
    std::swap(A,B);
    std::swap(DegA,DegB);
    std::swap(Acofac,Bcofac);
  }
  if (DegA != 0) CoCoA_THROW_ERROR("Not a field","Fq::div");

  const FpElem scale = myFFp.myRecip(A[0]);
  for (int i=0; i < len(Acofac); ++i)
    Acofac[i] = myFFp.myMul(scale, myFFp.myNormalize(Acofac[i]));
  
  myMul(ans,&Acofac[0],a);
}


  // n must be strictly positive, and no aliasing between lhs and x.
  void FFqImpl_vec::myBinaryPowerLoop(OUTvalue_t ans, value_t x, long n) const  // assumes n >= 1
  {
    CoCoA_ASSERT(n >= 1);
    if (n == 1) { for (int i=0; i<myDeg; ++i) { ans[i] = x[i]; } return; }
    if (n == 2) { myMul(ans,x,x); return; }
    myBinaryPowerLoop(ans, x, n/2);  // floor division!
    myMul(ans, ans, ans);//???    mySquare(rawlhs, rawlhs);
    if (n&1) myMul(ans, ans, x);
  }



  void FFqImpl_vec::myPower(OUTvalue_t ans, value_t x, long n) const
  {
    CoCoA_ASSERT(n >= 0);
    // Dispose of the trivial cases; note that 0^0 gives 1
//??    if (n == 0 || myIsOne(x)) { ans =myAssign(rawlhs, raw(myOne())); return; }
//??    if (n == 1 || myIsZero(x)) { ans = x; return; }

    myBinaryPowerLoop(ans, x, n);
  }


  
  std::vector<long> FFqImpl_vec::myExport(value_t x) const
  {
    vector<long> ans(myDeg);
    for (int i=0; i < myDeg; ++i)
      ans[i] = myFFp.myExport(x[i]);
    return ans;
  }


//   bool IsOne(const vector<SmallFpImpl::value>& x)
//   {
//     if (!IsOne(x[0])) return false;
//     const int n = len(x);
//     for (int i=1; i < n; ++i)
//       if (!IsZero(x[i])) return false;
//     return true;
//   }

//   long HasOrder(long n, const FFqImpl_vec& extn, const std::vector<SmallFpImpl::value>& x)
//   {
//     // const long p = ModP.myModulus();

//     // long q = SmallPower(p,d);
// //    factorization<long> qfacs = SmoothFactor(n,n);
//     const int d = extn.myExtnDeg();
//     const vector<long> primes = SmoothFactor(n,n).myFactors();
// ///    clog<<"HasOrder: n="<<n<<"   primes="<<primes<<endl;
// //???      vector<SmallFpImpl::value> X(d); X[1]=1;
// vector<SmallFpImpl::value> pwr(d); 
//     for (int i=0; i < len(primes); ++i)
//     {
//       const long e = n/primes[i];
// //      clog<<"Checking pwr e="<<e<<"   x="<<X<<endl;
//       extn.myPower(&pwr[0],&x[0],e);
//       if (IsOne(pwr)) return false;
//     }
//     return true;
// }

//   RingElem VectorToPoly(const ring& Fpx, const std::vector<SmallFpImpl::value>& v)
//   {
//     const long p = ConvertTo<long>(characteristic(Fpx)); // BUG BUG BUG   cheating!!!
//     SmallFpImpl ModP(p);
//     const long d = len(v)-1;
//     const RingElem x = indet(Fpx,0);
//     RingElem f(Fpx);
//     for (int i=d; i >= 0; --i)
//     {
//       f *= x;
//       f += ModP.myExport(v[i]);
//     }
//     return f;
//   }

//   bool IsGroupGen(const SmallFpImpl& ModP, const std::vector<SmallFpImpl::value>& MinPoly)
//   {
//     const long d = len(MinPoly)-1;
//     CoCoA_ASSERT(IsOne(MinPoly[d]));
//     if (d > 1 && IsZero(MinPoly[0])) return false;
//     const long p = ModP.myModulus();
//     ring Fp = NewZZmod(p);
//     ring Fpx = NewPolyRing(Fp, 1);
//     const RingElem x = indet(Fpx,0);
//     RingElem f = VectorToPoly(Fpx, MinPoly);
//     if (!IsIrred(f)) return false;
// ///    clog<<"IRRED OK"<<endl;

// //??    FFqImpl_vec extn(ModP,MinPoly);
//     FFqImpl_vec extn(p,MinPoly);
//     vector<SmallFpImpl::value> X(d, zero(SmallFp)); X[1]=one(SmallFp);
//     return HasOrder(SmallPower(p,d)-1,extn,X);
// }


//   void next(const SmallFpImpl& ModP, vector<SmallFpImpl::value>& v)
//   {
//     int i=0;
//     while (i < len(v))
//     {
//       v[i] = ModP.myAdd(v[i],one(SmallFp));
//       if (!IsZero(v[i])) return;
//       ++i;
//     }
//   }


//   std::vector<SmallFpImpl::value> FindGroupGenerator(long p, int d)
// {
//   SmallFpImpl ModP(p);
//   typedef SmallFpImpl::value FpElem;
//   vector<FpElem> f(1+d);
//   f[d] = one(SmallFp);
//   f[1] = one(SmallFp); // since gen cannot be of form x^d + const

// //  cout << "Looking for group gen..." << endl;
//   while (!IsGroupGen(ModP, f))
//     next(ModP,f);

// //  cout << "Group gen is " << f << endl;
//   return f;
// }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SmallFqVecImpl.C,v 1.8 2022/02/18 14:11:58 abbott Exp $
// $Log: SmallFqVecImpl.C,v $
// Revision 1.8  2022/02/18 14:11:58  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.7  2021/07/19 11:14:40  abbott
// Summary: Layout
//
// Revision 1.6  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.5  2021/01/07 15:16:53  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2020/06/17 15:49:27  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.3  2020/01/26 14:41:59  abbott
// Summary: Revised includes after splitting NumTheory (redmine 1161)
//
// Revision 1.2  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.1  2015/12/18 15:25:07  abbott
// Summary: Added impls of non-prime finite fields
//
//
