//   Copyright (c)  2022  John Abbott and Anna M. Bigatti
//   With contributions from Mohanad Baraghith (Uni Passau)

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

#include "CoCoA/MatrixOps.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/convert.H"
#include "CoCoA/long32or64.H"
#include "CoCoA/MatrixFp.H"

#include <algorithm>
using std::count;
//#include <vector>
using std::vector;

#ifdef CoCoA_32BIT_LONG
#define RND_PRIME_BITS 14
#else
#define RND_PRIME_BITS 30
#endif


namespace CoCoA
{

  bool HasLargeZeroBlock(const ConstMatrixView& M)
  {
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    const long n = NumRows(M);
    // Dispose of 2 easy small cases
    if (n == 0) return false;
    if (n == 1) return IsZero(M(0,0));

    // Make a "boolean copy" Mbool of M with entries just 0 or 1
    // if M(i,j) == 0 then Mbool(i,j) == 0; otherwise Mbool(i,j) == 1
    using BOOL = unsigned char; //  just needs to repr 0 or 1
 // using BOOL = bool;          //  which type is better??? 0 <--> false,  1 <--> true
    vector<vector<BOOL>> Mbool(n,vector<BOOL>(n));
    for (long i=0; i < n; ++i)
      for (long j=0; j < n; ++j)
        Mbool[i][j] = !IsZero(M(i,j));

    // Count how many zeroes in each row: ZeroCount[k] = num zeroes in k-th row
    vector<long> ZeroCount(n);
    for (long i=0; i < n; ++i)
      ZeroCount[i] = std::count(Mbool[i].begin(), Mbool[i].end(), BOOL(0));

    // Fill "row" with row indices such that
    // if i<j then ZeroCount[row[i]] >= ZeroCount[row[j]]
    vector<long> row(n);
    // First fill row = [0,1,2,3,4,...] then sort it
    for (long i=0; i < n; ++i) row[i] = i; // easier to read than std::generate(..., lambda)
    sort(row.begin(), row.end(), // sort by decreasing ZeroCount!
         [/*const*/ &ZeroCount](long i1, long i2) { return ZeroCount[i1] > ZeroCount[i2]; });

    if (ZeroCount[row[0]] == n) return true; // there is a zero row

    // Look for largest index which could corr to zero-block:
    long MaxCheckIndex = 0;
    for (long i=1; i < n; ++i)
      if (i+ZeroCount[row[i]] >= n)
        MaxCheckIndex = i+1;
    if (MaxCheckIndex == 0) return false;


    vector<BOOL> NZ(n); // at start of each iter below, NZ[k] is non-zero iff we have seen a non-zero in the k-th col
    // Iterate through rows in order dictated by array row
    for (long i=0; i < MaxCheckIndex; ++i)
    {
      long Zcount = 0;
      for (long j=0; j < n; ++j)
      {
        if (NZ[j] != 0) continue;
        if (Mbool[row[i]][j] == 0)  ++Zcount;
        else NZ[j] = 1;
      }
      if (Zcount + i > n) return true;
//???      if (Zcount + MaxCheckIndex < n) return false;
    }
    return false;


    // More efficient but harder to comprehend!!!  [not tested??]
    // Vector col contains col indices which we must check
    vector<long> col(n);
    for (long i=0; i < n; ++i)  // easier to read than std::generate(..., lambda)
      col[i] = i;

    long ncols = n;
    for (long i=0; i < n; ++i)
    {
      long OUT=0;
      long J = 0;
      while (J < ncols)
      {
        if (Mbool[row[i]][col[J]] == 0) { col[OUT] = col[J]; ++OUT; }
        ++J;
      }
      if (OUT == 0) return false;
      if (OUT+i > n) return true;
      ncols = OUT;
    }
    /*NEVER GET HERE*/    return false;
  }


  // Checks mod p for several large small primes p
  // Result: false3 -- det(M) is definitely not 0
  //         uncertain3 -- det(M) may be 0 or not
  //         true3 -- det(M) is 0 (can only happen if M is 1x1)
  bool3 IsZeroDet3_modular(const ConstMatrixView& M)
  {
    // Assume square
    // Assume all entries are INT/RAT
    if (!IsZero(characteristic(RingOf(M)))) return uncertain3;
    const size_t n = NumRows(M);
    if (n == 0) return false3;
    if (n == 1) return (IsZero(M(0,0)) ? true3 : false3);

    for (int trial=0; trial < 3; ++trial)
    {
      SmallPrime p = RandomNBitPrime(RND_PRIME_BITS);
      SmallFpImpl ModP(p);
      MatrixFp Mp(ModP,n,n);

      bool IsBadPrime = false;
      BigRat q;
      for (size_t i=0; i < n; ++i)
        for (size_t j=0; j < n; ++j)
        {
          if (!IsConvertible(q, M(i,j))) return uncertain3;
          if (IsZero(ModP.myReduce(den(q))))
            IsBadPrime = true;
          else
            Mp(i,j) = ModP.myReduce(ConvertTo<BigRat>(M(i,j)));
        }
      if (IsBadPrime) { --trial; continue; }
      if (det(Mp) != 0) return false3;
    }
    return uncertain3;
  }


  // Always fast; may simply return uncertain3
  bool3 IsZeroDet3(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);

    const size_t n = NumRows(M);
    if (n == 0) return false3;
    if (n == 1) return (IsZero(M(0,0)) ? true3 : false3);
    if (!IsZero(characteristic(RingOf(M))))
      return uncertain3;

    if (HasLargeZeroBlock(M)) return true3;
    return IsZeroDet3_modular(M);
  }


  // May be slow (if has to compute det, i.e. if IsZeroDet3 gives uncertain3)
  bool IsZeroDet(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    const bool3 QuickAns = IsZeroDet3(M);
    if (!IsUncertain3(QuickAns)) return (IsTrue3(QuickAns));
    return IsZero(det(M));
  }


}  // end of namespace CoCoA
