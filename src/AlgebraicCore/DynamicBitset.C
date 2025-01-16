//   Copyright (c)  2006-2010  Anna Bigatti, Massimo Caboara

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


#include "CoCoA/DynamicBitset.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/VectorOps.H" // for printing vectors


#include <algorithm>
using std::transform;
//#include <bitset>
using std::bitset;
#include <functional>
//??
#include <iterator>
using std::back_inserter;  // for transform
#include <iostream>
//using std::ostream;
//#include <vector>
using std::vector;

//static const bool MAX_DEBUG = false;

namespace CoCoA
{
  DynamicBitset::OutputStyle DynamicBitset::ourOutputStyle; // clean;

// class DynamicBitset


  void DynamicBitset::myResize(long n) // only for ctors
  {
    CoCoA_ASSERT(n >= 0);
    myLenValue = n;
    if (n == 0) return;
    myVec.resize((n-1)/ourNumBitsInBlock + 1); // all 0s
  }
  

  DynamicBitset::DynamicBitset(ConstRefPPMonoidElem pp)
  {
    const long n = NumIndets(owner(pp));
    vector<long> expv = exponents(pp);
    // same as DynamicBitset(expv)
    myResize(n);
    for (long i=0; i!=n; ++i)
      if (expv[i] != 0)  mySet(i);
  }


//   DynamicBitset::DynamicBitset(const vector<long>& v)
//   {
//     myResize(len(v));
//     for (long i=0; i!=len(v); ++i)
//       if (v[i] != 0)  mySet(i);
//   }


  DynamicBitset::DynamicBitset(long n)
  {
    if (n < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    myResize(n);
  }


  DynamicBitset::DynamicBitset(const DynamicBitset& rhs)
  {
    myLenValue = rhs.myLenValue;
    myVec = rhs.myVec;
  }


  // Anna 15Apr2010:
  // if ">>=" appears to be too slow make vector of precomputed masks (static)
  bool DynamicBitset::IamAll1s() const noexcept
  {
    constexpr BitBlock mask = ~0; // 111111111
   for (vector<BitBlock>::const_iterator it=myVec.begin(); it!=myVec.end(); ++it)
      if ((std::operator^(*it, mask)).any())  // xor for bitset
      {
        if (it+1 != myVec.end() || myLenValue%ourNumBitsInBlock==0) return false;
        const int shift = (ourNumBitsInBlock - myLenValue%ourNumBitsInBlock);
//JAA        mask >>= (ourNumBitsInBlock - myLenValue%ourNumBitsInBlock); // 00000111111
        if (std::operator^(*it, mask>>shift).any()) return false;
      }
    return true;
  }
  

  void DynamicBitset::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    if (myLenValue==0)  return;
    vector<BitBlock>::const_reverse_iterator rit=myVec.rbegin();
    // first block without spurious 0s
    for (long n=(myLenValue-1)%ourNumBitsInBlock; n>=0; --n)
      out << rit->test(n);
    //    for (++rit; rit!=myVec.rend(); ++rit)  out << '-' << *rit;
    for (++rit; rit!=myVec.rend(); ++rit) out << *rit;
  }


  void DynamicBitset::myOutputSelf8(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    if (myLenValue==0)  return;
    vector<BitBlock>::const_reverse_iterator rit=myVec.rbegin();
    // first block without spurious 0s
    for (long n=(myLenValue-1)%ourNumBitsInBlock; n>=0; --n)
    {
      out << rit->test(n);
      if (n%8==0 && n!=0)  out << '.';
    }
    for (++rit; rit!=myVec.rend(); ++rit)
    {
      out << '-';
      for (long n=ourNumBitsInBlock-1; n>=0; --n)
      {
        out << rit->test(n);
        if (n%8==0 && n!=0)  out << '.';
      }
    }
  }


  namespace // Anna 15Apr2010: maybe not the best way.  Just practising STL.
  {
    inline unsigned long ToULong(const DynamicBitset::BitBlock& b) noexcept
    { return b.to_ulong(); }
  }
  
  void DynamicBitset::myOutputSelfLong(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    vector<unsigned long> v;
    v.reserve(len(myVec));
    transform(myVec.rbegin(), myVec.rend(), back_inserter(v), ToULong);
    out << v;
  }


  long count(const DynamicBitset& b) noexcept
  {
    long NumOnes = 0;
    // for (vector<DynamicBitset::BitBlock>::const_iterator it=b.myVec.begin(); it!=b.myVec.end(); ++it)
    //   NumOnes += it->count();
    for (const auto& block: b.myVec)
      NumOnes += block.count();
    return NumOnes;
  }


  DynamicBitset flip(DynamicBitset DB)
  {
    const DynamicBitset::BitBlock All1s = ~0;
    const long w = DynamicBitset::ourNumBitsInBlock;
    const long lenDB = DB.myLenValue;
    const long n = lenDB/w; // integer division!
    for (long i=0; i < n; ++i)
      DB.myVec[i] ^= All1s;

    const int shift = (w - lenDB%w);
    if (shift != 0)
      DB.myVec.back() ^= (All1s >> shift);
    return DB;
  }


  std::ostream& operator<<(std::ostream& out, const DynamicBitset& DB)
  {
    if (!out) return out;  // short-cut for bad ostreams
    switch (DynamicBitset::ourOutputStyle)
    {
    case DynamicBitset::OutputStyle::clean:          DB.myOutputSelf(out);  break;
    case DynamicBitset::OutputStyle::AsRevVecOfLong: DB.myOutputSelfLong(out);  break;
    case DynamicBitset::OutputStyle::WithSeparators: DB.myOutputSelf8(out);  break;
    default: CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    }
    return out;
  }


  // can be improved making vector of exponents (for dense PPMonoid)
  PPMonoidElem NewPP(const PPMonoid& PPM, const DynamicBitset& b)
  {
    if (NumIndets(PPM) != len(b))  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    PPMonoidElem pp(PPM);
    for (long i=0; i!=len(b); ++i)
      if (b.Iam1At(i))  pp *= indet(PPM, i);
    return pp;
  }


}// end namespace cocoa





/*

Some future optimization:

for IsFLesser: proceed in the test word by word:
compute the first word of g1-f, the first word of g2-f, the check

ConnectionBlock: use ptr and not iterators.

sparse representation for facets? Use vector, (sort?).
Only reasonable if density is much much lower than #VARS!

*/
