#ifndef CoCoA_degree_H
#define CoCoA_degree_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti

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


// Header file for the class degree


#include "CoCoA/BigInt.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/utils.H"
// WARNING: next line must be last include (see redmine 591)
#include "CoCoA/VectorOps.H"  // for LexCmp3


#include <vector>
//using std::vector;
#include <iosfwd>
// using std::ostream;


namespace CoCoA
{

  // NOTE: most fns here are NOT noexcept because of arg checks

  namespace OrdvArith { class base; } // forward declaration for friend declaration

  class degree
  {
  public:
    explicit degree(long dim);
//???    friend long dim(const degree& deg);  // unnecessary because of GradingDim???
    const BigInt& operator[](long index) const;                      ///< deg[index] (READ ONLY), assumes 0 <= index < dim
    long myGradingDim() const noexcept;
    void mySetComponent(long index, const BigInt& N);                ///< deg[index] = N, assumes 0 <= index < dim
    void mySetComponent(long index, const MachineInt& n);            ///< deg[index] = n, assumes 0 <= index < dim
    degree& operator+=(const degree& d);                             ///< Computes deg1+=deg2
    degree& operator-=(const degree& d);                             ///< Computes deg1-=deg2
    friend int FastCmp(const degree& d1, const degree& d2) noexcept; ///< like cmp but inline and no checks
    friend inline bool IsZero(const degree& d) noexcept;
//    static void CheckCompatible(const degree& d1, const degree& d2, const char* fn); ///< checks dims are equal, throws if not.
    static void CheckCompatible(const degree& d1, const degree& d2, const ErrorContext& ErCtx); ///< checks dims are equal, throws if not.
  private: // data members
    std::vector<BigInt> myCoords;
  };


  long GradingDim(const degree& d) noexcept;                         ///< Computes "dimension" of the degree
  degree operator+(const degree& d1, const degree& d2);              ///< Computes d1+d2
  degree operator-(const degree& d1, const degree& d2);              ///< Computes d1-d2
  degree top(const degree& d1, const degree& d2);                    ///< i-th compt is max(d1[i], d2[i])
  int cmp(const degree& d1, const degree& d2);                       ///< lex comparison, result is <0,=0,>0
  void SetComponent(degree& d, long index, const BigInt& N);
  void SetComponent(degree& d, long index, const MachineInt& n);
  std::ostream& operator<<(std::ostream& out, const degree& d);      ///< prints out the degree, handles dim=1 specially.

  //---------------------------------------------------------------------------
  // inline functions

  inline degree::degree(long dim):
      myCoords(dim)
  {}


  inline long degree::myGradingDim() const noexcept
  {
    return len(myCoords);
  }


  // FastCmp is the same as cmp except it is inline and does not check compatibility
  inline int FastCmp(const degree& d1, const degree& d2) noexcept
  {
    const long dim = d1.myGradingDim();
    return LexCmp3(&d1.myCoords[0], &d1.myCoords[dim],
                   &d2.myCoords[0], &d2.myCoords[dim]);
  }


  inline bool operator<(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) < 0;
  }


  inline bool operator<=(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) <= 0;
  }


  inline bool operator>(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) > 0;
  }


  inline bool operator>=(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) >= 0;
  }


  inline bool operator==(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) == 0;
  }


  inline bool operator!=(const degree& d1, const degree& d2)
  {
    return cmp(d1, d2) != 0;
  }


  inline long GradingDim(const degree& d) noexcept
  {
    return d.myGradingDim();
  }

} // end of namespace CoCoA

#endif
