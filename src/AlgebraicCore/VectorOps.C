//   Copyright (c)  2016  John Abbott and Anna M. Bigatti

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

#include "CoCoA/VectorOps.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/utils.H"

#include<list>
using std::list;
#include<vector>
using std::vector;

namespace CoCoA
{

  // This should be a template fn?
  BigInt BinaryProduct(vector<BigInt>& v)
  {
    while (true)
    {
      const int n = len(v);
      if (n < 4) break;
      if (n&1) v[n-2] *= v[n-1];
      const int n2 = n/2; // integer division!
      for (int i=0; i < n2; ++i)
        v[i] = v[2*i]*v[2*i+1];
      v.resize(n2);
    }
    const int n = len(v);
    if (n == 1) return v[0];
    if (n == 2) return v[1]*v[1];
    return v[0]*v[1]*v[2];
  }

  BigInt BinaryProduct(const list<BigInt>& L)
  {
    if (L.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (len(L) == 1) return L.front();
    if (len(L) == 2) return L.front()*L.back();
    list<BigInt> tmp;
    const auto END = L.end();
    for (auto it = L.begin(); it != END; ++it)
    {
      const auto prev = it;
      ++it;
      if (it == END) { tmp.back() *= *it; break; }
      tmp.push_back((*prev)*(*it));
    }
    return BinaryProduct(tmp);
  }
  
} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/VectorOps.C,v 1.8 2024/07/15 17:12:42 bigatti Exp $
// $Log: VectorOps.C,v $
// Revision 1.8  2024/07/15 17:12:42  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.7  2022/02/18 14:12:01  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2021/01/07 15:16:53  abbott
// Summary: Corrected copyright
//
// Revision 1.5  2020/06/17 15:49:29  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.4  2020/02/18 11:30:16  abbott
// Summary: Some tidying; change var name l in to L
//
// Revision 1.3  2018/05/22 14:16:40  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.2  2018/05/18 12:24:47  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.1  2018/05/17 15:24:30  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.1  2017/02/01 12:17:57  abbott
// Summary: Added new fn BinaryProduct -- never called currently
//
//
