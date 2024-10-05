//   Copyright (c)  2014  John Abbott and Anna M. Bigatti

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

#include "CoCoA/UPoly.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/utils.H"


#include<algorithm>
using std::max;
#include <vector>
using std::vector;

// FOR DEBUGGING!!!
#include <iostream>
using std::cout;
using std::endl;

namespace CoCoA
{

  // test whether poly is of form a*x^n-b
  bool IsBinomial(const vector<BigInt>& c)
  {
    const long deg = len(c)-1;
    if (IsZero(c[0])) return false;
    for (long i=1; i < deg; ++i)
      if (!IsZero(c[i])) return false;
    // ASSUME !IsZero(c[deg])
    return true;
  }


  // Assumes poly is binomial, i.e. only first and last entries in c are non-zero.
  // Assumes n,d are coprime.
  void EvalBinomial(BigInt& AnsN, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    const long deg = len(c)-1;
    AnsN = c[deg]*power(n,deg) + c[0]*power(d,deg); // would be faster to use GMP fns directly!
  }

  BigInt horner(const vector<BigInt>& c, const BigInt& a)
  {
    const long deg = len(c)-1;
    BigInt ans = c[deg];
    for (long i=deg-1; i >= 0; --i)
      ans = a*ans+c[i];
    return ans;
  }

  BigInt HornerMPZ(const vector<BigInt>& c, const BigInt& a)
  {
    const long deg = len(c)-1;
    BigInt ans = c[deg];
    for (long i=deg-1; i >= 0; --i)
    {
      mpz_mul(mpzref(ans), mpzref(ans), mpzref(a));
      mpz_add(mpzref(ans), mpzref(ans), mpzref(c[i]));
    }
    return ans;
  }

  // BigRat HornerQQ(const vector<BigInt>& c, const BigRat& a)
  // {
  //   const long deg = len(c)-1;
  //   BigRat ans(c[deg],1);
  //   for (long i=deg-1; i >= 0; --i)
  //     ans = a*ans+c[i];
  //   return ans;
  // }

  void HornerQQ(BigInt& AnsN, BigInt& AnsD, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    const long deg = len(c)-1;
    AnsN = c[deg];
    AnsD = BigInt(1);
    for (long i=deg-1; i >= 0; --i)
    {
      AnsD *= d;
      AnsN = AnsN*n+c[i]*AnsD;
    }
  }

//just numerator
  void HornerQQ(BigInt& AnsN, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    const long deg = len(c)-1;
    AnsN = c[deg];
    BigInt AnsD(1);
    for (long i=deg-1; i >= 0; --i)
    {
      AnsD *= d;
      AnsN = AnsN*n+c[i]*AnsD;
    }
  }

//just numerator
  void HornerGMPQQ(BigInt& AnsN, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    const long deg = len(c)-1;
    AnsN = c[deg];
    BigInt AnsD(1);
    BigInt tmp; // workspace, to avoid repeated new-delete
    for (long i=deg-1; i >= 0; --i)
    {
      mpz_mul(mpzref(AnsD), mpzref(AnsD), mpzref(d));
      mpz_mul(mpzref(AnsN), mpzref(AnsN), mpzref(n));
      mpz_mul(mpzref(tmp), mpzref(c[i]), mpzref(AnsD));
      mpz_add(mpzref(AnsN), mpzref(AnsN), mpzref(tmp));
//    AnsD *= d;
//    AnsN = AnsN*n+c[i]*AnsD;
    }
  }

  void HornerRangeQQ(BigInt& AnsN, BigInt& AnsD, long lo, long hi, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    AnsN = c[hi];
    AnsD = BigInt(1);
    for (long i=hi-1; i >= lo; --i)
    {
      AnsD *= d;
      AnsN = AnsN*n+c[i]*AnsD;
    }
  }

// just numerator!!
  void HornerRangeQQ(BigInt& AnsN, long lo, long hi, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    AnsN = c[hi];
    BigInt PwrD(1);
    BigInt tmp; // workspace, to avoid repeated new-delete
    for (long i=hi-1; i >= lo; --i)
    {
//    D *= d;
//    AnsN = AnsN*n+c[i]*D;
      mpz_mul(mpzref(PwrD), mpzref(PwrD), mpzref(d));
      mpz_mul(mpzref(AnsN), mpzref(AnsN), mpzref(n));
      mpz_mul(mpzref(tmp), mpzref(c[i]), mpzref(PwrD));
      mpz_add(mpzref(AnsN), mpzref(AnsN), mpzref(tmp));
    }
  }


  BigInt HornerRange(long lo, long hi, const vector<BigInt>& c, const BigInt& a)
  {
    BigInt ans = c[hi];
    for (long i=hi-1; i >= lo; --i)
    {
//  EQUIV  ans = a*ans+c[i];
      mpz_mul(mpzref(ans), mpzref(ans), mpzref(a));
      mpz_add(mpzref(ans), mpzref(ans), mpzref(c[i]));
    }

    return ans;
  }

  void HornerRange(BigInt& ans, long lo, long hi, const vector<BigInt>& c, const BigInt& a)
  {
    ans = c[hi];
    for (long i=hi-1; i >= lo; --i)
    {
//  EQUIV  ans = a*ans+c[i];
      mpz_mul(mpzref(ans), mpzref(ans), mpzref(a));
      mpz_add(mpzref(ans), mpzref(ans), mpzref(c[i]));
    }
  }

  BigInt horner2(const vector<BigInt>& c, const BigInt& a)
  {
    if (len(c) < 4) return horner(c, a);
    const long half = len(c)/2;
    return HornerRange(0, half-1, c, a) + power(a, half)*HornerRange(half, len(c)-1, c, a);
  }

  BigInt HornerChunk1(long ChunkSize, const vector<BigInt>& c, const BigInt& a)
  {
    if (ChunkSize >= len(c)) return horner(c,a);
// assert(ChunkSize < len(c));
    long hi = len(c);
    long lo = ChunkSize*((hi-1)/ChunkSize);
    const BigInt pwr = power(a, ChunkSize);
    BigInt ans;
    while (lo >= 0)
    {
//    clog<<"Doing range: "<<lo<<" to "<<hi-1<<endl;
      ans = ans*pwr + HornerRange(lo,hi-1, c,a);
//    clog<<"ans="<<ans<<endl;
      hi = lo;
      lo -= ChunkSize;
    }
    return ans;
  }

  BigInt HornerChunk2(long ChunkSize, const vector<BigInt>& c, const BigInt& a)
  {
    if (ChunkSize >= len(c)) return HornerMPZ(c,a);
// assert(ChunkSize < len(c));
    long hi = len(c);
    long lo = ChunkSize*((hi-1)/ChunkSize);
    const BigInt pwr = power(a, ChunkSize);
    BigInt ans;
    BigInt tmp;
    while (lo >= 0)
    {
//    clog<<"Doing range: "<<lo<<" to "<<hi-1<<endl;
      mpz_mul(mpzref(ans), mpzref(ans), mpzref(pwr));
      HornerRange(tmp,lo,hi-1, c,a);
      mpz_add(mpzref(ans), mpzref(ans), mpzref(tmp));

//    clog<<"ans="<<ans<<endl;
      hi = lo;
      lo -= ChunkSize;
    }
    return ans;
  }

  void HornerChunkQQ(BigInt& AnsN, BigInt& AnsD, long ChunkSize, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    if (ChunkSize >= len(c)) return HornerQQ(AnsN, AnsD, c,n,d);
// assert(ChunkSize < len(c));
    long hi = len(c);
    long lo = ChunkSize*((hi-1)/ChunkSize);
    const BigInt PwrN = power(n, ChunkSize);
    const BigInt PwrD = power(d, ChunkSize);
    AnsN = 0;
    AnsD = 1;
    BigInt TmpN;
    BigInt TmpD;
    while (lo >= 0)
    {
//    clog<<"Doing range: "<<lo<<" to "<<hi-1<<endl;
      HornerRangeQQ(TmpN, TmpD, lo,hi-1, c,n,d);
      TmpN *= AnsD; // exact division!
      if (hi-lo == ChunkSize) AnsD *= PwrD;
      else AnsD = power(d, hi-lo);
      AnsN *= PwrN;
      AnsN += TmpN;
      hi = lo;
      lo -= ChunkSize;
    }
  }

// just numerator!!!
  void HornerChunkQQ(BigInt& AnsN, long ChunkSize, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    if (ChunkSize >= len(c)) return HornerGMPQQ(AnsN, c,n,d);
// assert(ChunkSize < len(c));
    long hi = len(c);
    long lo = ChunkSize*((hi-1)/ChunkSize);
    const BigInt PwrN = power(n, ChunkSize);
    const BigInt PwrD = power(d, ChunkSize);
    AnsN = 0;
    BigInt AnsD(1);
    BigInt TmpN;
    while (lo >= 0)
    {
//    clog<<"Doing range: "<<lo<<" to "<<hi-1<<endl;
      HornerRangeQQ(TmpN, lo,hi-1, c,n,d);
///    AnsN = AnsN*PwrN + AnsD*TmpN;
      mpz_mul(mpzref(AnsN), mpzref(AnsN), mpzref(PwrN));
      mpz_mul(mpzref(TmpN), mpzref(TmpN), mpzref(AnsD));
      mpz_add(mpzref(AnsN), mpzref(AnsN), mpzref(TmpN));
      if (lo == 0) return;
      if (hi-lo != ChunkSize) AnsD = power(d, hi-lo);
      else AnsD *= PwrD;
      // AnsD *= PwrD;
      // TmpN *= AnsD/TmpD; // exact division!
      // AnsN *= PwrN;
      // AnsN += TmpN;
      hi = lo;
      lo -= ChunkSize;
    }
  }

  long NearestPwr2(long n)
  {
    long ans=1;
    while (n > 3)
    {
      n /= 2;
      ans *= 2;
    }
    if (n == 3) return 4*ans;
    return n*ans;
  }


  BigInt HornerRecursiveIter(const vector<BigInt>& c, const BigInt& a)
  {
    const long NumBits = FloorLog2(a);
    const long TargetChunkSize = 2+ 1600/NumBits;
    if (len(c) < 2*TargetChunkSize) return HornerMPZ(c, a);
    if (len(c) < 4*TargetChunkSize) return HornerChunk2((len(c)+1)/2,c, a);
//???  const long NumChunks = 1 + (len(c)-1)/TargetChunkSize;
    const long NumChunks = NearestPwr2(1 + (len(c)-1)/TargetChunkSize);
//DEBUGGING    if(FirstTime){cout<<"NumChunks="<<NumChunks<<endl;FirstTime=false;}
    const long ChunkSize = 1 + (len(c)-1)/NumChunks;
    vector<BigInt> tmp(NumChunks);
//  clog<<"ChunkSize="<<ChunkSize<<endl;
//  clog<<"NumChunks="<<NumChunks<<endl;

    long hi = len(c);
    long i = ((hi-1)/ChunkSize);
    long lo = i*ChunkSize;
    while (lo >= 0)
    {
      HornerRange(tmp[i], lo,hi-1, c,a);
      hi=lo; lo -= ChunkSize;
      --i;
    }

//  clog<<"tmp="<<tmp<<endl;
    BigInt pwr = power(a,ChunkSize);
    while (len(tmp) > 1)
    {
      const long last = len(tmp)/2;
      for (long i=0; i < last; ++i)
      {
        mpz_mul(mpzref(tmp[2*i+1]), mpzref(tmp[2*i+1]), mpzref(pwr));
        mpz_add(mpzref(tmp[i]), mpzref(tmp[2*i]), mpzref(tmp[2*i+1]));
      }
      if (IsOdd(len(tmp))) mpz_swap(mpzref(tmp[last]), mpzref(tmp.back()));
      tmp.resize((1+len(tmp))/2);
      if (len(tmp) == 1) break;
      pwr *= pwr;
    }
    return tmp[0];
  }

// same as HornerRecursiveIter but uses double size initial chunks
  BigInt HornerRecursiveIter2(const vector<BigInt>& c, const BigInt& a)
  {
    const long NumBits = FloorLog2(a);
    const long TargetChunkSize = 2*(2+ 1600/NumBits);
    if (len(c) < 2*TargetChunkSize) return HornerMPZ(c, a);
    if (len(c) < 4*TargetChunkSize) return HornerChunk2((len(c)+1)/2,c, a);
//???  const long NumChunks = 1 + (len(c)-1)/TargetChunkSize;
    const long NumChunks = NearestPwr2(1 + (len(c)-1)/TargetChunkSize);
//DEBUGGING    if(FirstTime){cout<<"NumChunks="<<NumChunks<<endl;FirstTime=false;}
    const long ChunkSize = 1 + (len(c)-1)/NumChunks;
    vector<BigInt> tmp(NumChunks);
//  clog<<"ChunkSize="<<ChunkSize<<endl;
//  clog<<"NumChunks="<<NumChunks<<endl;

    long hi = len(c);
    long i = ((hi-1)/ChunkSize);
    long lo = i*ChunkSize;
    while (lo >= 0)
    {
      HornerRange(tmp[i], lo,hi-1, c,a);
      hi=lo; lo -= ChunkSize;
      --i;
    }

//  clog<<"tmp="<<tmp<<endl;
    BigInt pwr = power(a,ChunkSize);
    while (len(tmp) > 1)
    {
      const long last = len(tmp)/2;
      for (long i=0; i < last; ++i)
      {
        mpz_mul(mpzref(tmp[2*i+1]), mpzref(tmp[2*i+1]), mpzref(pwr));
        mpz_add(mpzref(tmp[i]), mpzref(tmp[2*i]), mpzref(tmp[2*i+1]));
      }
      if (IsOdd(len(tmp))) mpz_swap(mpzref(tmp[last]), mpzref(tmp.back()));
      tmp.resize((1+len(tmp))/2);
      if (len(tmp) == 1) break;
      pwr *= pwr;
    }
    return tmp[0];
  }

// Same as HornerRecursiveIter but handles odd length differently in main loop.
  BigInt HornerRecursiveIter3(const vector<BigInt>& c, const BigInt& a)
  {
    const long NumBits = FloorLog2(a);
    const long TargetChunkSize = 2+ 1600/NumBits;
    if (len(c) < 2*TargetChunkSize) return HornerMPZ(c, a);
    if (len(c) < 4*TargetChunkSize) return HornerChunk2((len(c)+1)/2,c, a);
    const long NumChunks = 1 + (len(c)-1)/TargetChunkSize;
//DEBUGGGING    if(FirstTime){cout<<"NumChunks="<<NumChunks<<endl;FirstTime=false;}
    const long ChunkSize = 1 + (len(c)-1)/NumChunks;
    vector<BigInt> tmp(NumChunks);
//  clog<<"ChunkSize="<<ChunkSize<<endl;
//  clog<<"NumChunks="<<NumChunks<<endl;

    long hi = len(c);
    long lo = ChunkSize*((hi-1)/ChunkSize);
    long i = NumChunks-1;
    while (lo >= 0)
    {
      HornerRange(tmp[i], lo,hi-1, c,a);
      hi=lo; lo -= ChunkSize;
      --i;
    }

//  clog<<"tmp="<<tmp<<endl;
    BigInt pwr = power(a,ChunkSize);
    while (len(tmp) > 1)
    {
      const long last = len(tmp)/2;
      if (IsOdd(len(tmp)))
      {
        mpz_mul(mpzref(tmp[2*last]), mpzref(tmp[2*last]), mpzref(pwr));
        mpz_add(mpzref(tmp[2*last-1]), mpzref(tmp[2*last-1]), mpzref(tmp[2*last]));
      }
      for (long i=0; i < last; ++i)
      {
        mpz_mul(mpzref(tmp[2*i+1]), mpzref(tmp[2*i+1]), mpzref(pwr));
        mpz_add(mpzref(tmp[i]), mpzref(tmp[2*i]), mpzref(tmp[2*i+1]));
      }
//    if (IsOdd(len(tmp))) mpz_swap(mpzref(tmp[last]), mpzref(tmp.back()));
//    tmp.resize((1+len(tmp))/2);
      tmp.resize(last);
      if (len(tmp) == 1) break;
      pwr *= pwr;
    }
    return tmp[0];
  }


  BigInt HornerRecursive(const vector<BigInt>& c, const BigInt& a)
  {
    if (len(c) <= 33) return HornerChunk2(2,c,a);
//  if (len(c) < 65) return HornerMPZ(c,a);
//  if (len(c) < 21) return HornerChunk1(4, c,a);
    const long ChunkSize = 4;
//  const long ChunkSize = (1+len(c))/2;  // forces 2 chunks
//  const long ChunkSize = 1+len(c)/(1+len(c)/32); // force chunks of size just larger than 32
    const long VecSize = 1+(len(c)-1)/ChunkSize;
    vector<BigInt> tmp(VecSize);
    long hi = len(c);
    long lo = ChunkSize*((hi-1)/ChunkSize);
    long i=VecSize-1;
    while (lo >= 0)
    {
//    tmp[i] = HornerRange(lo,hi-1, c,a);
      HornerRange(tmp[i], lo,hi-1, c,a);
      hi=lo; lo -= ChunkSize;
      --i;
    }

    return HornerRecursive(tmp, power(a,ChunkSize));
  }


  void HornerRecursiveQQ(BigInt& AnsN, BigInt& AnsD, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    if (len(c) < 21) return HornerChunkQQ(AnsN,AnsD, 4, c,n,d);
    const long ChunkSize = 4;
    const long VecSize = 1+(len(c)-1)/ChunkSize;
    vector<BigInt> tmp(VecSize);
    long hi = len(c);
    long lo = ChunkSize*((hi-1)/ChunkSize);
    long i=VecSize-1;
    BigInt TmpN;
    BigInt TmpD;
    while (lo >= 0)
    {
      HornerRangeQQ(TmpN, TmpD, lo,hi-1, c,n,d);
      tmp[i] = TmpN;
      if (hi-lo != ChunkSize) tmp[i] *= power(d, ChunkSize - (hi-lo));
      hi=lo; lo -= ChunkSize;
      --i;
    }

    HornerRecursiveQQ(AnsN, AnsD, tmp, power(n,ChunkSize), power(d,ChunkSize));
    AnsD *= power(d,ChunkSize);
  }

//just numerator!!
  void HornerRecursiveQQ(BigInt& AnsN, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    if (len(c) < 21) return HornerChunkQQ(AnsN, 1+(len(c)-1)/4, c,n,d);
    const long ChunkSize = 4;
    const long VecSize = 1+(len(c)-1)/ChunkSize;
    vector<BigInt> tmp(VecSize);
    long hi = len(c);
    long lo = ChunkSize*((hi-1)/ChunkSize);
    const BigInt ExtraFactor = power(d, ChunkSize - (hi-lo));
    long i=VecSize-1;
    BigInt TmpN;
//  const BigInt PwrD = power(d, ChunkSize);
    BigInt TmpD;
    while (lo >= 0)
    {
      HornerRangeQQ(TmpN, lo,hi-1, c,n,d);
      tmp[i] = TmpN;
      if (hi-lo != ChunkSize) tmp[i] *= ExtraFactor;
      hi=lo; lo -= ChunkSize;
      --i;
    }

    HornerRecursiveQQ(AnsN, tmp, power(n,ChunkSize), power(d,ChunkSize));
    AnsN /= ExtraFactor;
  }


  void HornerRecursiveIterQQ(BigInt& AnsN, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
//  const long NumBits = FloorLog2(n)+FloorLog2(d);
    const long NumBits = max(exponent(n),exponent(d));
    const long TargetChunkSize = (NumBits>999)?2:(1+2000/NumBits);
    const long LenC = len(c);
    if (LenC <= 2*TargetChunkSize) return HornerGMPQQ(AnsN, c, n,d);
    if (LenC <= 4*TargetChunkSize) return HornerChunkQQ(AnsN, 1+(len(c)-1)/4,c, n,d);
    const long NumChunks = 1 + (LenC-1)/TargetChunkSize;
    const long ChunkSize = 1 + (LenC-1)/NumChunks;
    vector<BigInt> tmp(NumChunks);
//  clog<<"ChunkSize="<<ChunkSize<<endl;
//  clog<<"NumChunks="<<NumChunks<<endl;

    BigInt ExcessFactor(1);
    long hi = LenC;
    long i = (hi-1)/ChunkSize;
    long lo = i*ChunkSize;
    if (hi-lo != ChunkSize) ExcessFactor = power(d, ChunkSize-(hi-lo));
    while (lo >= 0)
    {
      HornerRangeQQ(tmp[i], lo,hi-1, c, n,d);
      if (hi-lo != ChunkSize) tmp[i] *= ExcessFactor;
      hi=lo; lo -= ChunkSize;
      --i;
    }

//  clog<<"tmp="<<tmp<<endl;
    BigInt PwrN = power(n,ChunkSize);
    BigInt PwrD = power(d,ChunkSize);
    while (true)
    {
      const long last = len(tmp)/2;
      for (long i=0; i < last; ++i)
      {
        mpz_mul(mpzref(tmp[2*i+1]), mpzref(tmp[2*i+1]), mpzref(PwrN));
        mpz_mul(mpzref(tmp[2*i]), mpzref(tmp[2*i]), mpzref(PwrD));
        mpz_add(mpzref(tmp[i]), mpzref(tmp[2*i]), mpzref(tmp[2*i+1]));
      }
      if (IsOdd(len(tmp))) { mpz_mul(mpzref(tmp[last]), mpzref(tmp[2*last]), mpzref(PwrD)); ExcessFactor *= PwrD; }
      tmp.resize((1+len(tmp))/2);
      if (len(tmp) == 1) break;
      PwrN *= PwrN;
      PwrD *= PwrD;
    }
//???  mpz_divexact(mpzref(AnsN), mpzref(tmp[0]), mpzref(ExcessFactor));
    AnsN = tmp[0]/ExcessFactor;
  }


  void HornerRecursiveIterQQ2(BigInt& AnsN, const vector<BigInt>& c, const BigInt& n, const BigInt& d)
  {
    if (IsZero(n)) { AnsN = c[0]; return; }
    if (IsBinomial(c)) { EvalBinomial(AnsN, c, n, d); return; }
//    const long NumBits = max(exponent(n),exponent(d));
    const long TargetChunkSize = 2;
    const long LenC = len(c);
    if (LenC <= 7) return HornerGMPQQ(AnsN, c, n,d);
    if (LenC <= 23) return HornerChunkQQ(AnsN, 1+(len(c)-1)/4,c, n,d);
    const long NumChunks = 1 + (LenC-1)/TargetChunkSize;
    const long ChunkSize = 1 + (LenC-1)/NumChunks;
    vector<BigInt> tmp(NumChunks);
//  clog<<"ChunkSize="<<ChunkSize<<endl;
//  clog<<"NumChunks="<<NumChunks<<endl;

    BigInt ExcessFactor(1);
    long hi = LenC;
    long i = (hi-1)/ChunkSize;
    long lo = i*ChunkSize;
    if (hi-lo != ChunkSize) ExcessFactor = power(d, ChunkSize-(hi-lo));
    while (lo >= 0)
    {
      HornerRangeQQ(tmp[i], lo,hi-1, c, n,d);
      if (hi-lo != ChunkSize) tmp[i] *= ExcessFactor;
      hi=lo; lo -= ChunkSize;
      --i;
    }

//  clog<<"tmp="<<tmp<<endl;
    BigInt PwrN = power(n,ChunkSize);
    BigInt PwrD = power(d,ChunkSize);
    while (true)
    {
      const long last = len(tmp)/2;
      for (long i=0; i < last; ++i)
      {
        mpz_mul(mpzref(tmp[2*i+1]), mpzref(tmp[2*i+1]), mpzref(PwrN));
        mpz_mul(mpzref(tmp[2*i]), mpzref(tmp[2*i]), mpzref(PwrD));
        mpz_add(mpzref(tmp[i]), mpzref(tmp[2*i]), mpzref(tmp[2*i+1]));
      }
      if (IsOdd(len(tmp))) { mpz_mul(mpzref(tmp[last]), mpzref(tmp[2*last]), mpzref(PwrD)); ExcessFactor *= PwrD; }
      tmp.resize((1+len(tmp))/2);
      if (len(tmp) == 1) break;
      PwrN *= PwrN;
      PwrD *= PwrD;
    }
//???  mpz_divexact(mpzref(AnsN), mpzref(tmp[0]), mpzref(ExcessFactor));
    AnsN = tmp[0]/ExcessFactor;
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/UPoly.C,v 1.10 2022/02/18 14:12:01 abbott Exp $
// $Log: UPoly.C,v $
// Revision 1.10  2022/02/18 14:12:01  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.9  2021/01/07 15:16:53  abbott
// Summary: Corrected copyright
//
// Revision 1.8  2018/06/12 13:56:02  abbott
// Summary: Added comment
//
// Revision 1.7  2018/05/18 12:24:47  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.6  2017/10/16 19:57:51  abbott
// Summary: Deleted cruft
//
// Revision 1.5  2017/10/16 19:53:50  abbott
// Summary: Added new fns ChebyshevPoly, HermitePoly, LaguerrePoly
//
// Revision 1.4  2017/02/22 12:35:51  abbott
// Summary: Added comments
//
// Revision 1.3  2015/11/23 18:22:45  abbott
// Summary: Renamed ILogBase -> FloorLogBase; added FloorLog2, FloorLog10
//
// Revision 1.2  2014/05/06 13:14:40  abbott
// Summary: Commented out unused variable
// Author: JAA
//
// Revision 1.1  2014/03/07 14:20:21  abbott
// Summary: Added new UPoly (DUPZZ??) code
// Author: JAA
//
//
