#ifndef CoCoA_combinatorics_H
#define CoCoA_combinatorics_H

//   Copyright (c)  2015,2020  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/MachineInt.H"

#include <vector>
// using std::vector;

namespace CoCoA
{

  BigInt NumPartitions(const MachineInt& n);
  BigInt CatalanNumber(const MachineInt& n);

  // Random subsets of {0,1,2,...,n-1}
  std::vector<long> RandomSubsetIndices(const MachineInt& n);
  std::vector<long> RandomSubsetIndices(const MachineInt& n, const MachineInt& r);

  // Random r-tuple of elems from {0,1,2,...,n-1}
  std::vector<long> RandomTupleIndices(const MachineInt& n, const MachineInt& r);

  // Random permutation of tuple (0,1,2,...,n-1)
  std::vector<long> RandomPermutation(const MachineInt& n);

  int signature(const std::vector<int>& perm) noexcept;
  int signature(const std::vector<long>& perm) noexcept;


  // Iterator for subsets of {0,1,2,...,n-1}
  // 2023-01-01: Currently by cardinality then "lex"
  class SubsetIter
  {
  public:
    explicit SubsetIter(long n);
    SubsetIter(long n, long card);
    SubsetIter(const SubsetIter&) = default;
    SubsetIter& operator++();
    SubsetIter operator++(int);
    bool IamEnded() const noexcept { return (myN==-1); }
    const std::vector<long>& operator*() const noexcept { return myCurrSubset; }  // returns empty subset when ended
    friend std::ostream& operator<<(std::ostream& out, const SubsetIter& it);
  private: // data members
    long myN;    // -1 means iter has ended
    long myCard; // -1 means unrestricted cardinality
    std::vector<long> myCurrSubset;
  private: // impl details
    void MarkAsEnded();
  };

  inline bool IsEnded(const SubsetIter& it)  { return it.IamEnded(); }
  std::ostream& operator<<(std::ostream& out, const SubsetIter& it);


  // Iterator for tuples of {0,1,2,...,n-1} with cardinality k
  // 2024-02-08: Currently follows  "lex" order
  class TupleIter
  {
  public:
    TupleIter(long n, long card);
    TupleIter(const TupleIter&) = default;
    TupleIter& operator++();
    TupleIter operator++(int);
    bool IamEnded() const noexcept { return (myN==-1); }
    const std::vector<long>& operator*() const noexcept { return myCurrTuple; }  // returns empty subset when ended
    friend std::ostream& operator<<(std::ostream& out, const TupleIter& it);
  private: // data members
    long myN;    // -1 means iter has ended
    long myCard;
    std::vector<long> myCurrTuple;
  private: // impl details
    void MarkAsEnded();
  };

  inline bool IsEnded(const TupleIter& it)  { return it.IamEnded(); }
  std::ostream& operator<<(std::ostream& out, const TupleIter& it);

} // end of namespace CoCoA

#endif
