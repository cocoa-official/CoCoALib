//   Copyright (c)  2005-2008,2014  John Abbott and Anna M. Bigatti

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

#include "CoCoA/OrdvArith.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/PPOrdering.H"
#include "CoCoA/PREPROCESSOR_DEFNS.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/assert.H"
#include "CoCoA/config.H"
#include "CoCoA/convert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/verbose.H"


#include <iostream>
using std::ostream;
using std::endl;
//#include <vector>
using std::vector;
#include <limits>
using std::numeric_limits;
#include <algorithm>
using std::find;
using std::copy;
//using std::swap;


namespace CoCoA
{

  OrdvArith::base::base(long NumIndets, long GradingDim, long NumOrdvEntries):
    IntrusiveReferenceCount(),
    myNumIndets(NumIndets),
    myGradingDim(GradingDim),
    myOrdvBuffer(NumOrdvEntries),
    myExpvBuffer(NumIndets)
  {
    CoCoA_ASSERT(NumIndets > 0);
    CoCoA_ASSERT(NumIndets < 1000000); // complain about ridiculously large number of indets
    CoCoA_ASSERT(GradingDim <= NumIndets);
    CoCoA_ASSERT(NumOrdvEntries >= NumIndets);
    myBitsPerOrdvEntry = numeric_limits<SmallExponent_t>::digits;
    myPackingDensity = numeric_limits<OrdvElem>::digits / myBitsPerOrdvEntry;
    CoCoA_ASSERT(myPackingDensity >= 1);
    // Recompute myBitsPerOrdvEntry; this may increase the value (safely, of course)...
    myBitsPerOrdvEntry = numeric_limits<OrdvElem>::digits / myPackingDensity;
    CoCoA_ASSERT(myPackingDensity == 1 || myBitsPerOrdvEntry < numeric_limits<OrdvElem>::digits);
    if (myPackingDensity == 1) // special case because shift op not defined if shift >= wordsize
      myOrdvMask = numeric_limits<OrdvElem>::max();
    else
      myOrdvMask = (static_cast<OrdvElem>(1) << myBitsPerOrdvEntry) - 1;
    CoCoA_ASSERT(myOrdvMask != 0);
    // Reset myOrdvWords to the correct value...
    myOrdvWords = 1 + (NumOrdvEntries-1)/myPackingDensity;
    myOrdvWordsForCmp = 1 + (myNumIndets-1)/myPackingDensity;

    myRefCountZero();
  }


  OrdvArith::base::~base()
  {}



  void OrdvArith::base::myAssignZero(OrdvElem* ordv) const
  {
    for (long i=0; i < myOrdvWords; ++i)
      ordv[i] = 0;
  }


  void OrdvArith::base::myAssign(OrdvElem* dest, const OrdvElem* src) const
  {
//    std::copy(&src[0], &src[myOrdvWords], &dest[0]); // seems slightly slower
    for (long i=0; i < myOrdvWords; ++i)
      dest[i] = src[i];
  }


  void OrdvArith::base::mySwap(OrdvElem* ordv1, OrdvElem* ordv2) const
  {
//    if (ordv1 == ordv2) return; // worth checking this special case???
    for (long i=0; i < myOrdvWords; ++i)
      std::swap(ordv1[i], ordv2[i]);
  }


  void OrdvArith::base::myMulIndetPower(OrdvElem* ordv, long var, long exp) const
  {
    CoCoA_ASSERT(exp >= 0);
    CoCoA_ASSERT(0 <= var && var < myNumIndets);
    vector<long> PowerExpv(myNumIndets);  // should be member of OrdvArith::base?  Use myExpvBuffer???
    PowerExpv[var] = exp;
    vector<OrdvElem> PowerOrdv(myOrdvWords);  // should be member of OrdvArith::base?  Use myOrdvBuffer????
    myAssignFromExpv(&PowerOrdv[0], PowerExpv);
    myMul(ordv, ordv, &PowerOrdv[0]);
  }


  void OrdvArith::base::myPower(OrdvElem* ordv, const OrdvElem* ordv1, long LongExp) const
  {
    CoCoA_ASSERT(LongExp >= 0);
#ifdef CoCoA_DEBUG
    myPowerOverflowCheck(ordv1, LongExp);
#endif
    const OrdvElem exp = static_cast<OrdvElem>(LongExp);
    if (exp > myOrdvMask)  CoCoA_THROW_ERROR1(ERR::ExpTooBig);
    for (long i=0; i < myOrdvWords; ++i)
      ordv[i] = exp*ordv1[i];
  }


  void OrdvArith::base::myPowerOverflowCheck(const OrdvElem* ordv, long LongExp) const
  {
    if (LongExp == 0 || LongExp == 1) return;
    CoCoA_ASSERT(LongExp >= 0);
    if (static_cast<unsigned long>(LongExp) > myOrdvMask)
      CoCoA_THROW_ERROR1(ERR::ExpTooBig);
    const OrdvElem exp = static_cast<OrdvElem>(LongExp);
    // ??? Is it worth uncommenting these two shortcuts???
    // if (pow == 0) { myAssignZero(ordv); return; }
    // if (pow == 1) { myAssign(ordv, ordv1); return; }
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
    {
      // Check for ordv element overflow.
      if (myOrdvBuffer[i] > 0 && myOrdvMask/myOrdvBuffer[i] < exp)
        CoCoA_THROW_ERROR1(ERR::ExpTooBig);
    }
  }


  long OrdvArith::base::myStdDeg(const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpvBuffer(myNumIndets); // NB hides data member!
#endif
    myComputeExpv(myExpvBuffer, ordv);
    long d=0;
    for (long i=0; i < myNumIndets; ++i)
      d += myExpvBuffer[i];  // ignore possible overflow
    return d;
  }


  void OrdvArith::base::myWDeg(degree& d, const OrdvElem* ordv) const
  {
    CoCoA_ASSERT(GradingDim(d) == myGradingDim);
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myGradingDim);
    for (long i=0; i < myGradingDim; ++i)
      SetComponent(d, i, myOrdvBuffer[i]);
  }


  int OrdvArith::base::myCmpWDegPartial(const OrdvElem* ordv1, const OrdvElem* ordv2, long PartialGrDim) const  // assumes GrDim >= 0
  {
    CoCoA_ASSERT(0 <= PartialGrDim && PartialGrDim <= myGradingDim);
    const long last = PartialGrDim/myPackingDensity;
    for (long i=0; i < last; ++i)
      if (ordv1[i] != ordv2[i]) return (ordv1[i] > ordv2[i])?1:-1;
    const long ShiftCount = (last+1)*myPackingDensity - PartialGrDim;
    if (ShiftCount == myPackingDensity) return 0;
    CoCoA_ASSERT(myPackingDensity > 1 && myBitsPerOrdvEntry < numeric_limits<OrdvElem>::digits);
    // Reach here only if myPackingDensity > 1, so myBitsPerOrdvEntry < wordsize
    const OrdvElem last1 = ordv1[last] >> (ShiftCount*myBitsPerOrdvEntry);
    const OrdvElem last2 = ordv2[last] >> (ShiftCount*myBitsPerOrdvEntry);
    if (last1 == last2) return 0;
    if (last1 > last2) return 1;
    return -1;
  }


  bool OrdvArith::base::myIsZero(const OrdvElem* ordv) const
  {
    for (long i=0; i < myOrdvWordsForCmp; ++i)
      if (ordv[i] != 0) return false;
    return true;
  }


  // Simple rather than efficient.
  bool OrdvArith::base::myIsIndet(long& index, const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpvBuffer(myNumIndets); // NB hides data member!
#endif
    myComputeExpv(myExpvBuffer, ordv);

    long j = myNumIndets;
    for (long i = 0; i < myNumIndets; ++i)
    {
      if (myExpvBuffer[i] == 0) continue;
      if (j != myNumIndets || myExpvBuffer[i] != 1) return false;
      j = i;
    }
    if (j == myNumIndets) return false;
    index = j;
    return true;
  }


  OrdvArith::OrdvElem OrdvArith::base::myOrdvGetNth(const OrdvElem* ordv, long n) const
  {
    CoCoA_ASSERT(n < myOrdvWords*myPackingDensity);
    const long posn = n/myPackingDensity;
    const long n_shifts = myPackingDensity-1 - (n%myPackingDensity);
    const OrdvElem tmp = ordv[posn] >> (n_shifts * myBitsPerOrdvEntry);  // NB shift amount is less than word width!
    return (tmp & myOrdvMask);
  }


  void OrdvArith::base::myCompress(OrdvElem* ordv, const vector<OrdvElem>& buffer) const
  {
    if (myPackingDensity == 1)
    {
      std::copy(buffer.begin(), buffer.end(), ordv);
      return;
    }
    long posn = 0;
    for (long i=0; i < myOrdvWords; ++i)
    {
      OrdvElem word = 0; // this value is totally irrelevant, it gets shifted into "hyperspace"
      for (long j=0; j < myPackingDensity; ++j)
      {
        word <<= myBitsPerOrdvEntry; // ok because myBitsPerOrdvEntry < wordsize!!!
	if (posn < myNumIndets) word += buffer[posn];
	++posn;
      }
      ordv[i] = word;
    }
  }


  void OrdvArith::base::myDecompress(vector<OrdvElem>& buffer, const OrdvElem* ordv, long NumCompts) const
  {
    if (myPackingDensity == 1)
    {
      std::copy(&ordv[0], &ordv[NumCompts], buffer.begin());
      return;
    }
    long BasePosn = 0;
    for (long i=0; i < myOrdvWords; ++i)
    {
      OrdvElem word = ordv[i];
      for (long j=myPackingDensity; j-- > 0;)
      {
	if (BasePosn + j < NumCompts)
          buffer[BasePosn + j] = (word & myOrdvMask);
        word >>= myBitsPerOrdvEntry;  // ok because myBitsPerOrdvEntry < wordsize!!!
      }
      BasePosn += myPackingDensity;
    }
  }



  //---------------------------------------------------------------------------
  // LexImpl

  OrdvArith::LexImpl::LexImpl(long NumIndets):
    base(NumIndets, 0, NumIndets)
  {}


  void OrdvArith::LexImpl::myAssignFromExpv(OrdvElem* ordv, const vector<long>& expv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    for (long i=0; i < myNumIndets; ++i)
      myOrdvBuffer[i] = expv[i];
    myCompress(ordv, myOrdvBuffer);
  }


  void OrdvArith::LexImpl::myComputeExpv(vector<long>& expv, const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
      expv[i] = IntegerCast<long>(myOrdvBuffer[i]);
  }


  long OrdvArith::LexImpl::myExponent(const OrdvElem* ordv, long var) const
  {
    return IntegerCast<long>(myOrdvGetNth(ordv, var));
  }


  void OrdvArith::LexImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "OrdvArith::LexImpl(" << myNumIndets << ")";
  }



  //---------------------------------------------------------------------------
  // XelImpl

  OrdvArith::XelImpl::XelImpl(long NumIndets):
    base(NumIndets, 0, NumIndets)
  {}


  void OrdvArith::XelImpl::myAssignFromExpv(OrdvElem* ordv, const vector<long>& expv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    for (long i=0; i < myNumIndets; ++i)
      myOrdvBuffer[i] = expv[myNumIndets-i-1];
    myCompress(ordv, myOrdvBuffer);
  }


  void OrdvArith::XelImpl::myComputeExpv(vector<long>& expv, const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
      expv[i] = IntegerCast<long>(myOrdvBuffer[myNumIndets-i-1]);
  }


  long OrdvArith::XelImpl::myExponent(const OrdvElem* ordv, long var) const
  {
    return IntegerCast<long>(myOrdvGetNth(ordv, myNumIndets-var-1));
  }


  void OrdvArith::XelImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "OrdvArith::XelImpl(" << myNumIndets << ")";
  }



  //---------------------------------------------------------------------------
  // StdDegLexImpl

  OrdvArith::StdDegLexImpl::StdDegLexImpl(long NumIndets):
    base(NumIndets, 1, NumIndets)
  {}


  void OrdvArith::StdDegLexImpl::myAssignFromExpv(OrdvElem* ordv, const vector<long>& expv) const
  {
    OrdvElem deg = expv[0];
    for (long i=1; i < myNumIndets; ++i)
    {
      CoCoA_ASSERT("Exponent overflow" && deg <= myOrdvMask-expv[i]);
      deg += expv[i];
    }
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myOrdvBuffer[0] = deg;
    for (long i=1; i < myNumIndets; ++i)
      myOrdvBuffer[i] = expv[i-1];

    myCompress(ordv, myOrdvBuffer);
  }


  void OrdvArith::StdDegLexImpl::myComputeExpv(vector<long>& expv, const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    OrdvElem expN = myOrdvBuffer[0];
    for (long i=1; i < myNumIndets; ++i)
    {
      const OrdvElem& ordvi = myOrdvBuffer[i];
      expN -= ordvi;
      expv[i-1] = ordvi;
    }
    expv[myNumIndets-1] = expN;
  }


  long OrdvArith::StdDegLexImpl::myStdDeg(const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, 1);
    return myOrdvBuffer[0];
  }


  long OrdvArith::StdDegLexImpl::myExponent(const OrdvElem* ordv, long var) const
  {
    if (var < myNumIndets-1) return IntegerCast<long>(myOrdvGetNth(ordv, var+1));
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    OrdvElem ans = myOrdvBuffer[0];
    for (long i=1; i < myNumIndets; ++i)
      ans -= myOrdvBuffer[i];  // NB this cannot underflow if degree has not overflowed
    return IntegerCast<long>(ans);
  }


  void OrdvArith::StdDegLexImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "OrdvArith::StdDegLexImpl(" << myNumIndets << ")";
  }



  //---------------------------------------------------------------------------
  // StdDegRevLexImpl

  OrdvArith::StdDegRevLexImpl::StdDegRevLexImpl(long NumIndets):
    base(NumIndets, 1, NumIndets)
  {}


  void OrdvArith::StdDegRevLexImpl::myAssignFromExpv(OrdvElem* ordv, const vector<long>& expv) const
  {
    OrdvElem deg = expv[0];
    for (long i=1; i < myNumIndets; ++i)
    {
      CoCoA_ASSERT("Negative exponent" && expv[i] >= 0);
      CoCoA_ASSERT("Exponent overflow" && static_cast<unsigned long>(expv[i]) <= myOrdvMask-deg);
      deg += expv[i];
    }
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myOrdvBuffer[0] = deg;
    for (long i=1; i < myNumIndets; ++i)
      myOrdvBuffer[i] = deg - expv[myNumIndets - i];

    myCompress(ordv, myOrdvBuffer);
  }


  void OrdvArith::StdDegRevLexImpl::myComputeExpv(vector<long>& expv, const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    const OrdvElem deg = myOrdvBuffer[0];
    OrdvElem exp0 = (2-myNumIndets)*deg; // may HARMLESSLY become "negative" or "overflow"
    for (long i=1; i < myNumIndets; ++i)
    {
      const OrdvElem ordvi = myOrdvBuffer[i];
      exp0 += ordvi;
      expv[myNumIndets-i] = deg - ordvi;
    }
    expv[0] = exp0;
  }


  long OrdvArith::StdDegRevLexImpl::myStdDeg(const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, 1);
    return myOrdvBuffer[0];
  }


  long OrdvArith::StdDegRevLexImpl::myExponent(const OrdvElem* ordv, long var) const
  {
    if (var != 0) return IntegerCast<long>(myOrdvGetNth(ordv, 0) - myOrdvGetNth(ordv, myNumIndets-var));

#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    OrdvElem ans = (2-myNumIndets)*myOrdvBuffer[0]; // may HARMLESSLY become "negative" or overflow
    for (long i=1; i < myNumIndets; ++i)
      ans += myOrdvBuffer[i];
    return IntegerCast<long>(ans);
  }


  void OrdvArith::StdDegRevLexImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "OrdvArith::StdDegRevLexImpl(" << myNumIndets << ")";
  }



//   //---------------------------------------------------------------------------
//   // WDegRevLexImpl

//   OrdvArith::WDegRevLexImpl::WDegRevLexImpl(std::vector<long> W):
//     base(len(W), 1, len(W))
//   { myW = W; }


//   void OrdvArith::WDegRevLexImpl::myAssignFromExpv(OrdvElem* ordv, const vector<long>& expv) const
//   {
//     OrdvElem wdeg = expv[0];
//     for (long i=1; i < myNumIndets; ++i)
//     {
//       CoCoA_ASSERT("Negative exponent" && expv[i] >= 0);
//       CoCoA_ASSERT("Exponent overflow" && static_cast<unsigned long>(myW[i]*expv[i]) <= myOrdvMask-wdeg);
//       wdeg += W*expv[i];
//     }
// #ifdef CoCoA_THREADSAFE_HACK
//     vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
// #endif
//     myOrdvBuffer[0] = wdeg;
//     for (long i=1; i < myNumIndets; ++i)
//       myOrdvBuffer[i] = wdeg - myW[i]*expv[myNumIndets - i];

//     myCompress(ordv, myOrdvBuffer);
//   }


//   void OrdvArith::WDegRevLexImpl::myComputeExpv(vector<long>& expv, const OrdvElem* ordv) const
//   {
// #ifdef CoCoA_THREADSAFE_HACK
//     vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
// #endif
//     myDecompress(myOrdvBuffer, ordv, myNumIndets);
//     const OrdvElem wdeg = myOrdvBuffer[0];
//     OrdvElem WDegWithout0 = 0;
//     for (long i=1; i < myNumIndets; ++i)
//     {
//       const OrdvElem ordvi = myOrdvBuffer[i];
//       expv[myNumIndets-i] = (wdeg - ordvi)/myW[i];
//       WDegWithout0 += ordvi;
//     }
//     expv[0] = (wdeg-WDegWithout0)/myW[0];
//   }


//   long OrdvArith::WDegRevLexImpl::myWDeg(const OrdvElem* ordv) const
//   {
// #ifdef CoCoA_THREADSAFE_HACK
//     vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
// #endif
//     myDecompress(myOrdvBuffer, ordv, 1);
//     return myOrdvBuffer[0];
//   }


//   long OrdvArith::WDegRevLexImpl::myExponent(const OrdvElem* ordv, long var) const
//   {
//     if (var != 0) return IntegerCast<long>((myOrdvGetNth(ordv, 0) - myOrdvGetNth(ordv, myNumIndets-var))/myW[0]);

// #ifdef CoCoA_THREADSAFE_HACK
//     vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
// #endif
//     myDecompress(myOrdvBuffer, ordv, myNumIndets);
//     OrdvElem WDegWithout0 = 0;
//     for (long i=1; i < myNumIndets; ++i)
//       WDegWithout0 += myOrdvBuffer[i];
//     return IntegerCast<long>((wdeg-WDegWithout0)/myW[0]);
//   }


//   void OrdvArith::WDegRevLexImpl::myOutputSelf(ostream& out) const
//   {
//     out << "OrdvArith::WDegRevLexImpl(" << myW << ")";
//   }


 
  //---------------------------------------------------------------------------
  // StdDegRevLexImpl2

  OrdvArith::StdDegRevLexImpl2::StdDegRevLexImpl2(long NumIndets):
    base(NumIndets, 1, NumIndets)
  {}


  void OrdvArith::StdDegRevLexImpl2::myAssignFromExpv(OrdvElem* ordv, const vector<long>& expv) const
  {
    OrdvElem PartialSum = 0;
    long j = myNumIndets-1;
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    for (long i=0; i < myNumIndets; ++i, --j)
    {
      CoCoA_ASSERT("Exponent overflow" && IntegerCast<OrdvElem>(expv[i]) <= myOrdvMask-PartialSum);
      PartialSum += expv[i];
      myOrdvBuffer[j] = PartialSum;
    }

    myCompress(ordv, myOrdvBuffer);
  }


  void OrdvArith::StdDegRevLexImpl2::myComputeExpv(vector<long>& expv, const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    expv[0] = myOrdvBuffer[myNumIndets-1];
    long i = myNumIndets-1;
    for (long j=1; j < myNumIndets; ++j, --i)
      expv[i] = myOrdvBuffer[j-1] - myOrdvBuffer[j];
  }


  long OrdvArith::StdDegRevLexImpl2::myStdDeg(const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, 1);
    return myOrdvBuffer[0];
  }


  long OrdvArith::StdDegRevLexImpl2::myExponent(const OrdvElem* ordv, long var) const
  {
    if (var == 0) return IntegerCast<long>(myOrdvGetNth(ordv, myNumIndets-1));
    return IntegerCast<long>(myOrdvGetNth(ordv, myNumIndets-var-1) - myOrdvGetNth(ordv, myNumIndets-var));
  }


  void OrdvArith::StdDegRevLexImpl2::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "OrdvArith::StdDegRevLexImpl2(" << myNumIndets << ")";
  }



  //---------------------------------------------------------------------------
  // MatrixOrderingImpl

  OrdvArith::MatrixOrderingImpl::MatrixOrderingImpl(long NumIndets, long GradingDim, const ConstMatrixView& /*OrderMatrix*/):
    base(NumIndets, GradingDim, NumIndets)
  {
    CoCoA_THROW_ERROR1(ERR::NYI);
//     CoCoA_ASSERT(myGradingDim < NumRows(OrderMatrix));
//     CoCoA_ASSERT(NumRows(OrderMatrix) == NumIndets);
//     CoCoA_ASSERT(NumCols(OrderMatrix) == NumIndets);
//       // Check that the matrix entries are non-negative
//     for (long i=0; i < myNumIndets; ++i)
//       for (long j=0; j < myNumIndets; ++j)
//         CoCoA_ASSERT(myOrderMatrix[i][j] >= 0);
   }


  void OrdvArith::MatrixOrderingImpl::myAssignFromExpv(OrdvElem* ordv, const vector<long>& expv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    for (long i=0; i < myNumIndets; ++i)
    {
      myOrdvBuffer[i] = 0;
      for (long j=0; j < myNumIndets; ++j)
        myOrdvBuffer[i] += myOrderMatrix[i][j]*expv[j];
    }

    myCompress(ordv, myOrdvBuffer);
  }


  void OrdvArith::MatrixOrderingImpl::myComputeExpv(vector<long>& expv, const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
    {
      long deg = 0;
      for (long j=0; j < myNumIndets; ++j)
        deg += myAdjointOrderMatrix[i][j] * myOrdvBuffer[j];
      expv[i] = deg/myOrderMatrixDet;
    }
  }


  long OrdvArith::MatrixOrderingImpl::myExponent(const OrdvElem* ordv, long var) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    OrdvElem ans = 0;
    for (long j=0; j < myNumIndets; ++j)
      ans += myAdjointOrderMatrix[var][j] * myOrdvBuffer[j];
    ans /= myOrderMatrixDet;
    return IntegerCast<long>(ans);
  }



  void OrdvArith::MatrixOrderingImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    out << "OrdArith::MatrixOrdering(GradingDim=" << myGradingDim << ", ";
//     out << "PPOrdering(GRADING=matrix([";
//     for (long i=0; i < myGradingDim; ++i)
//     {
//       if (i > 0) out << ", ";
//       out << "[";
//       for (long j=0; j < myNumIndets; ++j)
//       {
//         if (j > 0) out << ", ";
//         out << myOrderMatrix[i][j];
//       }
//       out << "]";
//     }
//     out << "]), ";
    out << "matrix([";
    for (long i=0; i < myNumIndets; ++i) //??? start from myGradingDim???
    {
      if (i > 0) out << ", ";
      out << "[";
      for (long j=0; j < myNumIndets; ++j)
      {
        if (j > 0) out << ", ";
        out << myOrderMatrix[i][j];
      }
      out << "]";
    }
    out << "]))";
  }




  //--------------------  MatrixOrderingMod32749Impl --------------------

//   // copy ctor
//   OrdvArith::MatrixOrderingMod32749Impl::MatrixOrderingMod32749Impl(const MatrixOrderingMod32749Impl* copy):
//     base(copy->myNumIndets, copy->myGradingDim, len(copy->myOrdvBuffer))
//   {
//     //std::cout << "copy-ctor -- called --" << std::endl;
//     CoCoA_ASSERT(copy != NULL);
//     myOrderMatrix.resize(myNumIndets, vector<int>(myNumIndets));
//     myInverseOrderMatrix.resize(myNumIndets, vector<int>(myNumIndets));
//     for (long i=0; i < myNumIndets; ++i)
//       for (long j=0; j < myNumIndets; ++j)
//       {
//         myOrderMatrix[i][j] = copy->myOrderMatrix[i][j];
//         myInverseOrderMatrix[i][j] = copy->myInverseOrderMatrix[i][j];
//       }
//     myRefCountZero();
//   }


  OrdvArith::MatrixOrderingMod32749Impl::MatrixOrderingMod32749Impl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix):
    base(NumIndets, GradingDim, NumIndets)
  {
    VerboseLog VERBOSE("MatrixOrderingMod32749Impl ctor");
    VERBOSE(99) << "-- called --" << std::endl;
    CoCoA_ASSERT(NumRows(OrderMatrix) == NumIndets);
    CoCoA_ASSERT(NumCols(OrderMatrix) == NumIndets);
    CoCoA_ASSERT(IsZZ(RingOf(OrderMatrix))||IsQQ(RingOf(OrderMatrix)));
    if (!IsTermOrdering(OrderMatrix))  CoCoA_THROW_ERROR1(ERR::ReqTermOrdering);

    const QuotientRing Fp = NewRingFp(ourMODULUS, GlobalSettings::ResidueRepr::NonNegative);
    matrix M(NewDenseMat(Fp, NumIndets, NumIndets));
    const RingHom phi = CanonicalHom(RingOf(OrderMatrix), Fp);

    for (long i=0; i < GradingDim; ++i)
      for (long j=0; j < NumCols(OrderMatrix); ++j)
        if (sign(OrderMatrix(i,j)) < 0)
          CoCoA_THROW_ERROR2(ERR::NYI, "temporarily requiring weights to be non-negative");
    matrix PosOrdMat = MakeTermOrdMat(OrderMatrix);
    for (long i=0; i < myNumIndets; ++i)
      for (long j=0; j < myNumIndets; ++j)
      {
        if (PosOrdMat(i,j)>= ourMODULUS)  CoCoA_THROW_ERROR1(ERR::ArgTooBig);
        SetEntry(M, i, j, phi(PosOrdMat(i, j)));
      }
    matrix InvM = inverse(M);

    BigInt tmp;
    myOrderMatrix.resize(myNumIndets, vector<int>(myNumIndets));
    myInverseOrderMatrix.resize(myNumIndets, vector<int>(myNumIndets));
    for (long i=0; i < myNumIndets; ++i)
      for (long j=0; j < myNumIndets; ++j)
      {
        myOrderMatrix[i][j] = ConvertTo<long>(PosOrdMat(i,j));
        myInverseOrderMatrix[i][j] = ConvertTo<long>(InvM(i,j));
      }

#ifdef CoCoA_DEBUG
    // Verify that myOrderMatrix is all non-negative
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    for (long i=0; i < nrows; ++i)
      for (long j=0; j < ncols; ++j)
        CoCoA_ASSERT(myOrderMatrix[i][j] >= 0);
    // Verify that myOrderMatrix*myInverseOrderMatrix is the identity
    const matrix prod = M*InvM;
    for (long i=0; i < myNumIndets; ++i)
      for (long j=0; j < myNumIndets; ++j)
      {
        if (i == j) CoCoA_ASSERT("BAD INVERSE" && IsOne(prod(i,j)));
        else        CoCoA_ASSERT("BAD INVERSE" && IsZero(prod(i,j)));
      }
#endif
    myRefCountZero();
  }


  void OrdvArith::MatrixOrderingMod32749Impl::myAssignFromExpv(OrdvElem* ordv, const vector<long>& expv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    for (long i=0; i < myNumIndets; ++i)
    {
      myOrdvBuffer[i] = 0;
      for (long j=0; j < myNumIndets; ++j)
        myOrdvBuffer[i] += myOrderMatrix[i][j]*expv[j];
    }

    myCompress(ordv, myOrdvBuffer);
  }


  void OrdvArith::MatrixOrderingMod32749Impl::myComputeExpv(vector<long>& expv, const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
    {
      unsigned long deg = 0;
      for (long j=0; j < myNumIndets; ++j)
      {
        deg += (myInverseOrderMatrix[i][j] * myOrdvBuffer[j]);
        if (deg>46336*ourMODULUS) deg -= 46336*ourMODULUS;
      }
      expv[i] = deg%ourMODULUS;
    }
  }


  long OrdvArith::MatrixOrderingMod32749Impl::myExponent(const OrdvElem* ordv, long var) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    OrdvElem ans = 0;
    for (long j=0; j < myNumIndets; ++j)
    {
      ans += (myInverseOrderMatrix[var][j] * myOrdvBuffer[j]);
      if (ans > 46336*ourMODULUS) ans -= 46336*ourMODULUS;
    }
    return ans%ourMODULUS; // no need to use IntegerCast here, overflow can never occur
  }


  void OrdvArith::MatrixOrderingMod32749Impl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    out << "OrdvArith::MatrixOrdering(GradingDim=" << myGradingDim << ", ";
//     out << "PPOrdering(GRADING=matrix([";
//     for (long i=0; i < myGradingDim; ++i)
//     {
//       if (i > 0) out << ", ";
//       out << "[";
//       for (long j=0; j < myNumIndets; ++j)
//       {
//         if (j > 0) out << ", ";
//         out << myOrderMatrix[i][j];
//       }
//       out << "]";
//     }
//     out << "]), ";
    out << "matrix([";
    for (long i=0; i < myNumIndets; ++i)
    {
      if (i > 0) out << ", ";
      out << "[";
      for (long j=0; j < myNumIndets; ++j)
      {
        if (j > 0) out << ", ";
        out << myOrderMatrix[i][j];
      }
      out << "]";
    }
    out << "]))";
  }


//   void OrdvArith::MatrixOrderingMod32749Impl::mySetMatrix(const matrix& M)
//   {
//     CoCoA_ASSERT(NumRows(M) == myNumIndets);
//     CoCoA_ASSERT(NumCols(M) == myNumIndets);
//     BigInt tmp;

//     for (long i=0; i < myNumIndets; ++i)
//       for (long j=0; j < myNumIndets; ++j)
//       {
//         if (!IsInteger(tmp, M(i, j)))
//           CoCoA_THROW_ERROR1("entry of MatrixOrdering is not integer");
//         if (!convert(myInverseOrderMatrix[i][j], tmp))
//           CoCoA_THROW_ERROR1("entry of MatrixOrdering is not integer");
//       }
//   }


//   void OrdvArith::MatrixOrderingMod32749Impl::mySetInverseMatrixTmp(const matrix& /*M*/)
//   { /*???*/  }



  //--------------------  MatrixOrdering32bitImpl --------------------

  OrdvArith::MatrixOrdering32bitImpl::MatrixOrdering32bitImpl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix):
    base(NumIndets, GradingDim, NumIndets)
  {
    VerboseLog VERBOSE("MatrixOrdering32bitImpl ctor");
    VERBOSE(99) << "-- called --" << std::endl;
    CoCoA_ASSERT(NumRows(OrderMatrix) == NumIndets);
    CoCoA_ASSERT(NumCols(OrderMatrix) == NumIndets);
    CoCoA_ASSERT(IsZZ(RingOf(OrderMatrix))||IsQQ(RingOf(OrderMatrix)));
    if (!IsTermOrdering(OrderMatrix))  CoCoA_THROW_ERROR1(ERR::ReqTermOrdering);

    const QuotientRing Fp = NewRingFp(ourMODULUS, GlobalSettings::ResidueRepr::NonNegative);
    matrix M(NewDenseMat(Fp, NumIndets, NumIndets));
    const RingHom phi = CanonicalHom(RingOf(OrderMatrix), Fp);

    for (long i=0; i < GradingDim; ++i)
      for (long j=0; j < NumCols(OrderMatrix); ++j)
        if (sign(OrderMatrix(i,j)) < 0)
          CoCoA_THROW_ERROR2(ERR::NYI, "temporarily requiring weights to be non-negative");
    matrix PosOrdMat = MakeTermOrdMat(OrderMatrix);
    for (long i=0; i < myNumIndets; ++i)
      for (long j=0; j < myNumIndets; ++j)
      {
        if (PosOrdMat(i,j) >= ourMODULUS)  CoCoA_THROW_ERROR1(ERR::ArgTooBig);
        SetEntry(M, i, j, phi(PosOrdMat(i, j)));
      }
    const matrix InvM = inverse(M);

    myOrderMatrix.resize(myNumIndets, vector<INT>(myNumIndets));
    myInverseOrderMatrix.resize(myNumIndets, vector<INT>(myNumIndets));
    for (long i=0; i < myNumIndets; ++i)
      for (long j=0; j < myNumIndets; ++j)
      {
        myOrderMatrix[i][j] = ConvertTo<INT>(PosOrdMat(i,j));
        myInverseOrderMatrix[i][j] = ConvertTo<INT>(InvM(i,j));
      }

#ifdef CoCoA_DEBUG
    // Verify that myOrderMatrix is all non-negative
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    for (long i=0; i < nrows; ++i)
      for (long j=0; j < ncols; ++j)
        CoCoA_ASSERT(myOrderMatrix[i][j] >= 0);
    // Verify that myOrderMatrix*myInverseOrderMatrix is the identity
    const matrix prod = M*InvM;
    for (long i=0; i < myNumIndets; ++i)
      for (long j=0; j < myNumIndets; ++j)
      {
        if (i == j) CoCoA_ASSERT("BAD INVERSE" && IsOne(prod(i,j)));
        else        CoCoA_ASSERT("BAD INVERSE" && IsZero(prod(i,j)));
      }
#endif  // CoCoA_DEBUG
    myRefCountZero();
  }


  void OrdvArith::MatrixOrdering32bitImpl::myAssignFromExpv(OrdvElem* ordv, const vector<long>& expv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    for (long i=0; i < myNumIndets; ++i)
    {
      myOrdvBuffer[i] = 0;
      for (long j=0; j < myNumIndets; ++j)
        myOrdvBuffer[i] += myOrderMatrix[i][j]*expv[j];
    }

    myCompress(ordv, myOrdvBuffer);
  }


  void OrdvArith::MatrixOrdering32bitImpl::myComputeExpv(vector<long>& expv, const OrdvElem* ordv) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
    {
      unsigned long deg = 0;
      for (long j=0; j < myNumIndets; ++j)
      {
        deg += (myInverseOrderMatrix[i][j] * myOrdvBuffer[j]);
        if (deg > 46336*ourMODULUS)   deg -= 46336*ourMODULUS;
      }
      expv[i] = deg%ourMODULUS;
    }
  }


  long OrdvArith::MatrixOrdering32bitImpl::myExponent(const OrdvElem* ordv, long var) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<OrdvElem> myOrdvBuffer(myNumIndets); // NB hides data member!
#endif
    myDecompress(myOrdvBuffer, ordv, myNumIndets);
    OrdvElem ans = 0;
    for (long j=0; j < myNumIndets; ++j)
    {
      ans += (myInverseOrderMatrix[var][j] * myOrdvBuffer[j]);
      if (ans > 46336*ourMODULUS)   ans -= 46336*ourMODULUS;
    }
    return ans%ourMODULUS; // no need to use IntegerCast here, overflow can never occur
  }


  void OrdvArith::MatrixOrdering32bitImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    out << "OrdvArith::MatrixOrdering(GradingDim=" << myGradingDim << ", ";
//     out << "PPOrdering(GRADING=matrix([";
//     for (long i=0; i < myGradingDim; ++i)
//     {
//       if (i > 0) out << ", ";
//       out << "[";
//       for (long j=0; j < myNumIndets; ++j)
//       {
//         if (j > 0) out << ", ";
//         out << myOrderMatrix[i][j];
//       }
//       out << "]";
//     }
//     out << "]), ";
    out << "matrix([";
    for (long i=0; i < myNumIndets; ++i)
    {
      if (i > 0) out << ", ";
      out << "[";
      for (long j=0; j < myNumIndets; ++j)
      {
        if (j > 0) out << ", ";
        out << myOrderMatrix[i][j];
      }
      out << "]";
    }
    out << "]))";
  }


//   void OrdvArith::MatrixOrdering32bitImpl::mySetMatrix(const matrix& M)
//   {
//     CoCoA_ASSERT(NumRows(M) == myNumIndets);
//     CoCoA_ASSERT(NumCols(M) == myNumIndets);
//     BigInt tmp;

//     for (long i=0; i < myNumIndets; ++i)
//       for (long j=0; j < myNumIndets; ++j)
//       {
//         if (!IsInteger(tmp, M(i, j)))
//           CoCoA_THROW_ERROR1("entry of MatrixOrdering is not integer");
//         if (!convert(myInverseOrderMatrix[i][j], tmp))
//           CoCoA_THROW_ERROR1("entry of MatrixOrdering is not integer");
//       }
//   }


//   void OrdvArith::MatrixOrdering32bitImpl::mySetInverseMatrixTmp(const matrix& /*M*/)
//   { /*???*/  }



  //---------------------------------------------------------------------------

  OrdvArith::reference NewOrdvArith(const PPOrdering& PPO)
  {
    VerboseLog VERBOSE("NewOrdvArith(ord)");
    VERBOSE(99) << "-- called --" << std::endl;
    if (IsLex(PPO))
      return OrdvArith::reference(new OrdvArith::LexImpl(NumIndets(PPO)));
    if (IsStdDegLex(PPO))
      return OrdvArith::reference(new OrdvArith::StdDegLexImpl(NumIndets(PPO)));
    if (IsStdDegRevLex(PPO))
      return OrdvArith::reference(new OrdvArith::StdDegRevLexImpl(NumIndets(PPO)));

    // If we get here, we have a matrix ordering.

    const long n = NumIndets(PPO);
    const long g = GradingDim(PPO);

    ConstMatrixView M(OrdMat(PPO));

    CoCoA_ASSERT(NumRows(M) == n);
    CoCoA_ASSERT(NumCols(M) == n);
    CoCoA_ASSERT(g <= n);

    return OrdvArith::reference(new OrdvArith::MatrixOrdering32bitImpl(n, g, M));

    // (1) Get matrix M out of the ordering.
    // (2) Make an equivalent matrix M2 which is strictly positive
    // (3) Build a MatrixOrderingImpl object with M2, but also need
    //     the transformation matrix to be able to calculate degrees!!
  }


  std::ostream& operator<<(std::ostream& out, const OrdvArith::reference& OA)
  {
    if (!out) return out;  // short-cut for bad ostreams
    OA->myOutputSelf(out);
    return out;
  }


} // end of namespace CoCoA
