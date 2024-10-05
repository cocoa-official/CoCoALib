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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/combinatorics.H,v 1.11 2024/02/08 20:46:30 abbott Exp $
// $Log: combinatorics.H,v $
// Revision 1.11  2024/02/08 20:46:30  abbott
// Summary: Added new class TupleIter (redmine 379)
//
// Revision 1.10  2023/12/10 20:26:40  abbott
// Summary: Added SubsetIter with restricted cardinality (redmine 379)
//
// Revision 1.9  2023/01/01 11:36:06  abbott
// Summary: Added SubsetIter
//
// Revision 1.8  2022/08/05 12:02:41  abbott
// Summary: Added CatalanNumber
//
// Revision 1.7  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.5  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2020/09/22 18:11:22  abbott
// Summary: Added signature (of a perrm); made NumPartitions interruptible
//
// Revision 1.3  2020/01/26 14:16:52  abbott
// Summary: Moved NumPartitions from NumTheory to here
//
// Revision 1.2  2018/08/28 12:35:41  abbott
// Summary: Added new fn RandomPermutation
//
// Revision 1.1  2015/05/20 13:38:32  abbott
// Summary: New files for combinatorial fns: RandomSubset, RandomTuple
// Author: JAA
//
//

#endif
