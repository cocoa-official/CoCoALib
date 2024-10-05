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

#include "CoCoA/SmallFqLogImpl.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/SmallFqUtils.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"

//#include <vector>
using std::vector;

namespace CoCoA
{


  FFqLogImpl::FFqLogImpl(const ideal& I):
      myModulus(ConvertTo<long>(characteristic(RingOf(I)))),
      myDeg(deg(gens(I)[0])),
      myCard(SmallPower(myModulus, myDeg)),
      myAdd1Tbl(myCard),
      myExpTbl(myCard),
      mySmallLogTbl(myModulus)
  {
    MakeArithTbls(I);
  }

  void FFqLogImpl::myMulByGen(const SmallFpImpl& FFp, std::vector<SmallFpImpl::value>& v, std::vector<SmallFpImpl::value>& MonicMinPoly) const
  {
    const SmallFpImpl::value c = FFp.myNegate(v[myDeg-1]);
//    std::clog<<"myMulByGen: v="<<v<<std::endl;
//    std::clog<<"myMulByGen: c="<<c<<std::endl;
    for (int i=myDeg-1; i > 0; --i)
    {
      v[i] = FFp.myAddMul(v[i-1], c, MonicMinPoly[i]);
    }
    v[0] = FFp.myMul(c, MonicMinPoly[0]);
//    std::clog<<"myMulByGen: v="<<v<<std::endl;
  }


  long FFqLogImpl::compress(const SmallFpImpl& FFp, const std::vector<SmallFpImpl::value>& v) const noexcept
  {
//    std::clog<<"compress: arg="<<v<<std::endl;
    long ans = 0;
    for (int i=myDeg-1; i >= 0; --i)
      ans = ans*myModulus + FFp.myExportNonNeg(v[i]);
//    std::clog<<"compress: ans="<<ans<<std::endl;
    return ans;
  }
  // std::vector<SmallFpImpl::value> PolyToVec(const RingElem& m, const SmallFpImpl& FFp)
  // {
  //   const int d = deg(m);
  //   vector<SmallFpImpl::value> coeffs(1+d, zero(SmallFp)FFp.myReduce(0));
  //   const vector<RingElem> C = CoeffVecWRT(m, indet(owner(m),0));
  //   for (int i=0; i <= d; ++i)
  //     coeffs[i] = FFp.myReduce(ConvertTo<long>(C[i]));
  //   return coeffs;
  // }

  void FFqLogImpl::MakeArithTbls(const ideal& I)
  {
    SmallFpImpl FFp(myModulus);
    vector<SmallFpImpl::value> MonicMinPoly = PolyToVec(gens(I)[0], FFp); // **ASSUME** principal
    
    vector<long> LogTbl(myCard);
    vector<SmallFpImpl::value> pwr(myDeg);
    pwr[0] = one(SmallFp);
//    std::clog<<"MakeLogExpTbls: pwr="<<pwr<<std::endl;
    // First build log and exp tables.
    for (long exp = 0; exp < myCard-1; ++exp)
    {
      const long CompressedVal = compress(FFp, pwr);
      CoCoA_ASSERT(CompressedVal != 0);
      CoCoA_ASSERT(LogTbl[CompressedVal] == 0); // slot has not yet been assigned
      LogTbl[CompressedVal] = exp;
      if (CompressedVal < myModulus) mySmallLogTbl[CompressedVal] = exp;
      myExpTbl[exp] = CompressedVal;
//    std::clog<<"MakeLogExpTbls: pwr="<<pwr<<std::endl;
      myMulByGen(FFp, pwr, MonicMinPoly);
//    std::clog<<"MakeLogExpTbls: pwr="<<pwr<<std::endl;
    }
//    std::clog<<"LogTbl="<<LogTbl<<std::endl;
//    std::clog<<"ExpTbl="<<ExpTbl<<std::endl;

    // Now build the "add 1" table; skip entry for -1 which will never be used.
    for (repr_t i=0; i < myCard-1; ++i)
    {
      const long val = myExpTbl[i];
      if (val == myModulus-1) continue; // skip -1
//      if (val == myModulus-1) std::clog<<"SURPRISE i="<<i<<std::endl;
      if (val%myModulus != myModulus-1)
      {
//        std::clog<<"val+1="<<val+1<<std::endl;
        myAdd1Tbl[i] = LogTbl[val+1];
      }
      else
      {
//        std::clog<<"val+1-p="<<val+1-myModulus<<std::endl;
        myAdd1Tbl[i] = LogTbl[val+1-myModulus];
      }
    }
//    std::clog<<"Add1Tbl="<<myAdd1Tbl<<std::endl;
  }

  // ***READ THIS***
  // Repr is 1+log(val); and 0 is repr as 0!
  
  FFqLogImpl::repr_t FFqLogImpl::myGen() const noexcept
  {
    return 1+1;
  }
  
  FFqLogImpl::repr_t FFqLogImpl::myReduce(long n) const noexcept
  {
    const long r = LeastNNegRemainder(n, myModulus);
    if (r == 0) return 0;
    return 1+mySmallLogTbl[r];
  }

  FFqLogImpl::repr_t FFqLogImpl::myNegate(repr_t x) const noexcept
  {
    if (myIsZero(x) || myModulus == 2) return x;
    const repr_t half = myCard/2; // integer division!
    if (x > half)
      return x - half;
    return x + half;
  }


  FFqLogImpl::repr_t FFqLogImpl::myAdd(repr_t x, repr_t y) const noexcept
  {
    if (myIsZero(x)) return y;
    if (myIsZero(y)) return x;
    const repr_t q = myDiv(y,x);
    if (myIsMinusOne(q)) return 0;
    return x + myAdd1Tbl[q-1];
  }

  FFqLogImpl::repr_t FFqLogImpl::mySub(repr_t x, repr_t y) const noexcept
  {
    if (myIsZero(x)) return y;
    if (myIsZero(y)) return x;
    if (x == y) return 0;
    const repr_t q = myNegate(myDiv(y,x));
    return x + myAdd1Tbl[q-1];
  }
    
  FFqLogImpl::repr_t FFqLogImpl::myMul(repr_t x, repr_t y) const noexcept
  {
    if (myIsZero(x) || myIsZero(y)) return 0;
    const repr_t ans = x+y-1;
    if (ans < myCard) return ans;
    return ans - (myCard-1);
  }
  
  FFqLogImpl::repr_t FFqLogImpl::myDiv(repr_t x, repr_t y) const
  {
    if (myIsZero(y)) CoCoA_THROW_ERROR(ERR::DivByZero, "FFqLogImpl::myDiv");
    if (myIsZero(x)) return 0;
    if (x >= y) return x-y+1;
    return x + (myCard-y);
  }
  
  FFqLogImpl::repr_t FFqLogImpl::myPower(repr_t x, long n) const noexcept // assumes n >= 0
  {
    CoCoA_ASSERT(n >= 0);
    if (n == 0) return 1;
    if (myIsZero(x)) return x;
    const long q1 = myCard-1;
    return 1 + (static_cast<long>(x-1)*(n%q1))%q1; // no overflow provided q1^2 < MaxLong
  }

//  long FFqLogImpl::myExtnDeg() const { return myDeg; }
  long FFqLogImpl::myExport(repr_t x) const noexcept
  {
    if (myIsZero(x)) return 0;
    return myExpTbl[x-1];
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SmallFqLogImpl.C,v 1.8 2022/10/25 17:18:41 abbott Exp $
// $Log: SmallFqLogImpl.C,v $
// Revision 1.8  2022/10/25 17:18:41  abbott
// Summary: Changed int to repr_t when running over field elems
//
// Revision 1.7  2022/02/18 14:11:58  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
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
// Revision 1.3  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.2  2016/01/27 14:03:26  abbott
// Summary: Added some comments and some assertions
//
// Revision 1.1  2015/12/18 15:25:07  abbott
// Summary: Added impls of non-prime finite fields
//
//
