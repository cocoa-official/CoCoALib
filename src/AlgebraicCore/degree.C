//   Copyright (c)  2002-2009,2011  John Abbott and Anna M. Bigatti

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

#include "CoCoA/degree.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
#include<algorithm>
using std::max;
//#include <vector>
using std::vector;

namespace CoCoA
{

  // inline void degree::CheckCompatible(const degree& d1, const degree& d2, const char* fn)
  // {
  //   if (len(d1.myCoords) != len(d2.myCoords))
  //     CoCoA_THROW_ERROR(ERR::MixedDegrees, fn);
  // }

  inline void degree::CheckCompatible(const degree& d1, const degree& d2, const ErrorContext& ErrCtx)
  {
    if (len(d1.myCoords) != len(d2.myCoords))
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::MixedDegrees, ErrCtx);
  }


  const BigInt& degree::operator[](long index) const
  {
    if (index < 0 || index >= len(myCoords))
      CoCoA_THROW_ERROR(ERR::BadIndex, "degree::operator[]");
    return myCoords[index];
  }


  void degree::mySetComponent(long index, const BigInt& N)
  {
    CoCoA_ASSERT(0 <= index && index < len(myCoords));
    myCoords[index] = N;
  }


  void degree::mySetComponent(long index, const MachineInt& n)
  {
    CoCoA_ASSERT(0 <= index && index < len(myCoords));
    myCoords[index] = BigInt(n);
  }


  degree& degree::operator+=(const degree& d)
  {
    degree::CheckCompatible(*this, d, CoCoA_ERROR_CONTEXT);
    const long dim = len(myCoords);
    for (long i=0; i < dim; ++i)
      myCoords[i] += d[i];
    return *this;
  }


  degree& degree::operator-=(const degree& d)
  {
    degree::CheckCompatible(*this, d, CoCoA_ERROR_CONTEXT);
    const long dim = len(myCoords);
    for (long i=0; i < dim; ++i)
      myCoords[i] -= d[i];
    return *this;
  }


  bool IsZero(const degree& d) noexcept
  {
    const long dim = GradingDim(d);
    //??? return find_if(coords.begin(), coords.end(), NonZero) != coords.end();
    for (long i=0; i < dim; ++i)
      if (d[i] != 0) return false;
    return true;
  }


  degree operator+(const degree& d1, const degree& d2)
  {
    degree::CheckCompatible(d1, d2, CoCoA_ERROR_CONTEXT);
    const long dim = GradingDim(d1);
    degree ans(dim);
    for (long i=0; i < dim; ++i)
      ans.mySetComponent(i, d1[i] + d2[i]);
    return ans;
  }


  degree operator-(const degree& d1, const degree& d2)
  {
    degree::CheckCompatible(d1, d2, CoCoA_ERROR_CONTEXT);
    const long dim = GradingDim(d1);
    degree ans(dim);
    for (long i=0; i < dim; ++i)
      ans.mySetComponent(i, d1[i] - d2[i]);
    return ans;
  }


  degree top(const degree& d1, const degree& d2)
  {
    degree::CheckCompatible(d1, d2, CoCoA_ERROR_CONTEXT);
    const long dim = GradingDim(d1);
    degree ans(dim);
    for (long i=0; i < dim; ++i)
      ans.mySetComponent(i, max(d1[i], d2[i]));
    return ans;
  }


  int cmp(const degree& d1, const degree& d2)
  {
    degree::CheckCompatible(d1, d2, CoCoA_ERROR_CONTEXT);
    // The rest is the same as the body of FastCmp; writing it explicitly
    // here avoids compiler foibles (i.e. choosing not to make FastCmp inline).
//     const long dim = GradingDim(d1);
    return FastCmp(d1, d2);
//     return LexCmp3(&d1.myCoords[0], &d1.myCoords[dim],
//                    &d2.myCoords[0], &d2.myCoords[dim]);
//     for (long i=0; i < dim; ++i)
//       if (d1[i] != d2[i])
//         return (d1[i] > d2[i] ? 1 : -1);
//     return 0;
  }


  ostream& operator<<(ostream& out, const degree& d)
  {
    if (!out) return out;  // short-cut for bad ostreams
    const long dim = GradingDim(d);
    if (dim == 0) return out << "()";     // no grading -- graded over N^0
    //  if (dim == 1) return out << d[0]; // omit parens if grading is over N
    // General case
    out << "(";
    for (long i=0; i < dim-1; ++i)
      out << d[i] << ", ";
    out << d[dim-1] << ")";
    return out;
  }

  void SetComponent(degree& d, long index, const BigInt& N)
  {
    if (index < 0 || index >= GradingDim(d))
      CoCoA_THROW_ERROR(ERR::BadIndex, "SetComponent(degree, index, N)");
    d.mySetComponent(index, N);
  }

  void SetComponent(degree& d, long index, const MachineInt& n)
  {
    if (index < 0 || index >= GradingDim(d))
      CoCoA_THROW_ERROR(ERR::BadIndex, "SetComponent(degree, index, n)");
    d.mySetComponent(index, n);
  }


} // end of namespace CoCoA
