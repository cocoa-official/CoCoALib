//   Copyright (c)  2015,2025  John Abbott and Anna M. Bigatti

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


#include "CoCoA/combinatorics.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/random.H"
#include "CoCoA/utils.H"

#include <algorithm>
//using std::swap;
using std::sort;
#include <limits>
using std::numeric_limits;
#include <set>
using std::set;
//#include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {

    // Procedure fills slots 0 to N (included!)
    // Source: http://www.mathpages.com/home/kmath383.htm
    // IDEA: could easily write a version which EXTENDS an existing table.
    void NumPartitionsTbl(vector<BigInt>& tbl, long N)
    {
      tbl.resize(N+1);
      tbl[0] = 1;
      for (long j=1; j <= N; ++j)
      {
        CheckForInterrupt("NumPartitionsTbl");
        int sign = 1;
        long k = 1;
        long i = 1;
        BigInt sum;
        while (true)
        {
          if (i > j) break;
          if (sign == 1) sum += tbl[j-i]; else sum -= tbl[j-i];
          i += k; // should not overflow if tbl fits into RAM
          if (i > j) break; 
          if (sign == 1) sum += tbl[j-i]; else sum -= tbl[j-i];
          sign = -sign;
          ++k;
          i += 2*k-1; // should not overflow if tbl fits into RAM
        }
        tbl[j] = sum;
      }
    }

  } // end of anonymous

  BigInt NumPartitions(const MachineInt& n)
  {
    if (IsNegative(n)) return BigInt(0); // or give error???
    if (IsZero(n)) return BigInt(1);
    const long N = AsSignedLong(n);
    vector<BigInt> tbl(N+1);
    NumPartitionsTbl(tbl, N);
    return tbl[N];
  }


  //-------------------------------------------------------

  BigInt CatalanNumber(const MachineInt& n)
  {
    if (IsNegative(n))
      CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    const unsigned long N = uabs(n);
    if (N > std::numeric_limits<unsigned long>::max()/2)
      CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    return binomial(2*N,N)/(N+1);
  }
  
  //-------------------------------------------------------

  std::vector<long> RandomSubsetIndices(const MachineInt& n)
  {
    if (IsNegative(n))  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    const long N = AsSignedLong(n);
    vector<long> ans;
    for (int i=0; i < N; ++i)
      if (RandomBool())
        ans.push_back(i);
    return ans;
  }


  namespace // anonymous
  {
    // This version is slow when r is not small, but requires
    // little memory when r is small and n is large.
    // It could also easily be adapted to make "small" random
    // subsets of BigInt values.
    std::vector<long> RandomSubsetIndices_slow(long n, long r)
    {
      std::set<long> S;
      long size = 0;
      while (size < r)
      {
        const long val = RandomLong(0,n-1);
        if (S.find(val) != S.end())  continue;
        ++size;
        S.insert(val);
      }
      std::vector<long> elems; elems.reserve(r);
      for (const auto& val: S)
        elems.push_back(val);
      return elems; // seems to be sorted
    }

  } // end of namespace anonymous


  std::vector<long> RandomSubsetIndices(const MachineInt& n, const MachineInt& r)
  {
    if (IsNegative(n))  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    const long N = AsSignedLong(n);
    long R = AsSignedLong(r);
    if (R < 0 || R > N)  CoCoA_THROW_ERROR1(ERR::OutOfRange);
    if (R <= 2+ N/512) // reducing 512 makes the program require less RAM (but also rather more time)
      return RandomSubsetIndices_slow(N,R); // more memory efficient
    std::vector<long> elems; elems.reserve(R); // do this early to trigger std::bad_alloc quickly if the answer is too large
    const bool invert = (R > N/2);
    if (invert)  R = N-R;
    std::vector<bool> member(N);
    long size = 0;
    while (size < R)
    {
      const long val = RandomLong(0,N-1);
      if (member[val])  continue;
      ++size;
      member[val] = true;
    }
    for (long i=0; i < N; ++i)
      if (invert^member[i])
        elems.push_back(i);
    return elems;
  }


  std::vector<long> RandomTupleIndices(const MachineInt& n, const MachineInt& r)
  {
    if (IsNegative(n) || IsNegative(r))  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    const long N = AsSignedLong(n);
    const long R = AsSignedLong(r);
    vector<long> ans(R);
    for (long i=0; i < R; ++i)
      ans[i] = RandomLong(0,N-1); // both extremes included
    return ans;
  }


  // Taken from Wikipedia (Random permutation); algorithm "Knuth Shuffle"
  std::vector<long> RandomPermutation(const MachineInt& n)
  {
    if (IsNegative(n))  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    const long N = AsSignedLong(n);
    vector<long> ans(N);
    for (int i=0; i < N; ++i)
      ans[i] = i;
    for (int i=N-1; i > 0; --i)
    {
      const int j = RandomLong(0,i);
      if (j != i)
        std::swap(ans[i], ans[j]);
    }
    
    return ans;
  }


  // BUG? Should this be a template fn?
  int signature(const std::vector<int>& perm) noexcept
  {
    int ans = 1;
    const int n = len(perm);
    for (int i=0; i < n; ++i)
      for (int j=i+1; j < n; ++j)
        if (perm[i] > perm[j])  ans = -ans;
    return ans;
  }

  int signature(const std::vector<long>& perm) noexcept
  {
    int ans = 1;
    const int n = len(perm);
    for (int i=0; i < n; ++i)
      for (int j=i+1; j < n; ++j)
        if (perm[i] > perm[j])  ans = -ans;
    return ans;
  }


  //--------------------------------------------
  // Iterator for subsets of {0,1,2,...,n-1}

  // All subsets
  SubsetIter::SubsetIter(long n):
      myN(n),
      myCard(-1),
      myCurrSubset()
  { if (n < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative); }


  // Only subsets of specified cardinality
  SubsetIter::SubsetIter(long n, long card):
      myN(n),
      myCard(card),
      myCurrSubset()
  {
    if (n < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    if (card < 0 || card > n)  CoCoA_THROW_ERROR1(ERR::OutOfRange);
    // Now place initial subset of given card in myCurrSubset
    myCurrSubset.resize(myCard);
    for (int j=0; j < myCard; ++j)
      myCurrSubset[j] = j;
  }


  inline void SubsetIter::MarkAsEnded()
  {
    myCurrSubset.clear();
    myN = -1;
  }

  SubsetIter& SubsetIter::operator++()
  {
    if (myN == -1) { return *this; } // do nothing if ended
    const long size = myCurrSubset.size();
    if (size == myN)  { MarkAsEnded();  return *this; }
    if (size == 0) { myCurrSubset.push_back(0); return *this; }
    if (myCurrSubset.back() != myN-1)  { ++myCurrSubset.back(); return *this; }
    int i = size-2;
    while (i >= 0)
    {
      if (myCurrSubset[i]+1 != myCurrSubset[i+1]) break;
      --i;
    }
    if (i >= 0)
    {
      ++myCurrSubset[i];
      while (++i < size)
      { myCurrSubset[i] = 1+myCurrSubset[i-1]; }
      return *this;
    }
    if (myCard >= 0)  { MarkAsEnded();  return *this; }
    myCurrSubset.push_back(size);
    for (int j=0; j < size; ++j)
      myCurrSubset[j] = j;
    return *this;
  }

  
  SubsetIter SubsetIter::operator++(int)
  {
    SubsetIter copy(*this);
    operator++();
    return copy;
  }
  

  std::ostream& operator<<(std::ostream& out, const SubsetIter& it)
  {
    if (!out) return out;
    out << "SubsetIter(myN=" << it.myN << ")"; // print also myCurrSubset???
    return out;
  }


  //--------------------------------------------
  // Iterator for k-tuples of {0,1,2,...,n-1}

  // Only subsets of specified cardinality
  TupleIter::TupleIter(long n, long card):
      myN(n),
      myCard(card),
      myCurrTuple()
  {
    if (n < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    if (card < 0)  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    myCurrTuple.resize(myCard);     // Fill with initial zero tuple
  }


  inline void TupleIter::MarkAsEnded()
  {
    myCurrTuple.clear();
    myN = -1;
  }


  TupleIter& TupleIter::operator++()
  {
    if (myN == -1) { return *this; } // do nothing if ended
    int i = 0;
    while (i < myCard)
    {
      if (myCurrTuple[i] != myN-1)  break;
      myCurrTuple[i] = 0;
      ++i;
    }
    if (i == myCard)  { MarkAsEnded(); return *this; }
    ++myCurrTuple[i];
    return *this;
  }

  TupleIter TupleIter::operator++(int)
  {
    TupleIter copy(*this);
    operator++();
    return copy;
  }
  

  std::ostream& operator<<(std::ostream& out, const TupleIter& it)
  {
    if (!out) return out;
    out << "TupleIter(myN=" << it.myN << ")"; // print also myCurrTuple???
    return out;
  }


} // end of namespace CoCoA
