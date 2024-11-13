#ifndef CoCoA_VectorOps_H
#define CoCoA_VectorOps_H

//   Copyright (c)  2005,2006,2009  John Abbott, and Anna M. Bigatti

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


// This contains sundry things which may eventually become stable parts of
// the CoCoA library.


#include "CoCoA/error.H"

#include <algorithm>
//using std::for_each;
//using std::find_if;
#include <iostream>
#include <list>
// using std::list
#include <vector>
// using std::vector;

namespace CoCoA
{

  ////////////////////////////////////////////////////////////
  // Basic pseudo-constructors

//   template <typename T>
//   inline std::vector<T> CoCoAVector(const T& a)
//   { return std::vector<T>(1,a); }

//   template <typename T>
//   inline std::vector<T> CoCoAVector(const T& a, const T& b)
//   { std::vector<T> v(1,a); v.push_back(b); return v; }

//   template <typename T>
//   inline std::vector<T> CoCoAVector(const T& a, const T& b, const T& c)
//   { std::vector<T> v(1,a); v.push_back(b); v.push_back(c); return v; }

//   template <typename T>
//   inline std::vector<T> CoCoAVector(const T& a,const T& b,const T& c,const T& d)
//   { std::vector<T> v(1,a);v.push_back(b);v.push_back(c);v.push_back(d); return v;}


  template <typename T>
  void/*std::vector<T>*/ concat_copy(std::vector<T> dest, std::vector<T> src)
  {
    dest.insert( dest.end(), src.begin(), src.end() );
  }


  template <typename T>
  void/*std::vector<T>*/ concat_move(std::vector<T> dest, std::vector<T> src) // destroys src!!!
  {
    dest.insert(dest.end(),
                std::make_move_iterator(src.begin()),
                std::make_move_iterator(src.end())
                );
  }


  ///////////////////////////////////////////////////////
  // Functions for printing vectors/lists/etc.

  // LLVM 6.0 (MacOS 10.10) needs this fwd decl
  template <typename FwdIterator_t>
  void OutputRange(std::ostream& out, FwdIterator_t start, const FwdIterator_t& end);


  // Use cout << v syntax to print out a std::vector<T> of printable values.
  template <typename T>
  inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
  {
    OutputRange(out, v.begin(), v.end());
    return out;
  }

  // Use cout << l syntax to print out a std::list<T> of printable values.
  template <typename T>
  inline std::ostream& operator<<(std::ostream& out, const std::list<T>& v)
  {
    OutputRange(out, v.begin(), v.end());
    return out;
  }

  // This prints out the elements (of a list or vector) in a given range.
  // The empty list appears as []; otherwise ",  " separated elements
  // between "[" and "]".
  // NOTE: This template fn is not intended for public use!
  template <typename FwdIterator_t>
  void OutputRange(std::ostream& out, FwdIterator_t start, const FwdIterator_t& end)
  {
    if (!out) return;
    if (start == end) { out << "[]"; return; }
    out << "[" << *start;
    while (++start != end)
      out << ",  " << *start;
    out << "]";
  }


  ////////////////////////////////////////////////////////////
  // Functions for computing products of lists and vectors.

  // This class is simply to permit use of std::for_each.
  template <typename arith_t>
  class SeqProd
  {
  public:
    SeqProd(const arith_t& InitialValue): myProdSoFar(InitialValue) {}
    // Default copy ctor, assignment, and dtor are all OK.
    void operator()(const arith_t& NextValue) { myProdSoFar *= NextValue; }
    const arith_t& myCurrentValue() const { return myProdSoFar; }
  private:
    arith_t myProdSoFar;
  };

  template <typename InputIterator, typename arith_t>
  arith_t product(InputIterator begin, InputIterator end, const arith_t& InitialValue)
  {
    return std::for_each(begin, end, SeqProd<arith_t>(InitialValue)).myCurrentValue();
  }

  template <typename ListOrVector>
  typename ListOrVector::value_type product(const ListOrVector& V)
  {
    if (V.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    return product(++V.begin(), V.end(), V.front());
  }


  ////////////////////////////////////////////////////////////
  // Functions for computing sums of lists and vectors.

  // This class is simply to permit use of for_each.
  template <typename arith_t>
  class SeqSum
  {
  public:
    SeqSum(const arith_t& InitialValue): mySumSoFar(InitialValue) {}
    // Default copy ctor, assignment, and dtor are all OK.
    void operator()(const arith_t& NextValue) { mySumSoFar += NextValue; }
    const arith_t& myCurrentValue() const { return mySumSoFar; }
  private:
    arith_t mySumSoFar;
  };

  template <typename InputIterator, typename arith_t>
  arith_t sum(InputIterator begin, InputIterator end, const arith_t& InitialValue)
  {
    return std::for_each(begin, end, SeqSum<arith_t>(InitialValue)).myCurrentValue();
  }

  template <typename ListOrVector>
  typename ListOrVector::value_type sum(const ListOrVector& V)
  {
    if (V.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    return sum(++V.begin(), V.end(), V.front());
  }


////////////////////////////////////////////////////////////
// Function for checking all elements have the same owner

  template <typename container>
  bool HasUniqueOwner(const container& v) noexcept
  {
    using elem_t = typename container::value_type;
    if (v.empty()) return true;  // ??? or error ???
    const auto& R = owner(v.front());
    auto OwnerIsR = [&R](const elem_t& x){ return (owner(x) == R);};
    return v.end() ==
      find_if_not(++v.begin(), v.end(), OwnerIsR);
  }


//////////////////////////////////////////////////////////////////
//  3-way comparison

//   // Generic 3-way comparison function: assumes operator< is defined.
//   // Result is <0 =0 or >0 according as a<b, a==b, a>b.
//   template <typename T1, typename T2>
//   inline int cmp(const T1& a, const T2& b)
//   {
//     if (a < b) return -1;
//     if (b < a) return 1;
//     return 0;
//   }

  // LexCmp3 3-way lex comparison.  I think this ought to be in the STL.
  // Result is <0, =0, or >0 according as begin1..end1 < = > begin2..end2
  // Assumes existence of fn cmp for single elements.
  template <class InputIterator1, class InputIterator2>
  inline int LexCmp3(InputIterator1 begin1, InputIterator1 end1,
                     InputIterator2 begin2, InputIterator2 end2) noexcept
  {
    while (begin1 != end1 && begin2 != end2)
    {
      const int sign = cmp(*begin1, *begin2);
      if (sign != 0) return sign;
      ++begin1;
      ++begin2;
    }
    if (begin1 != end1) return 1;
    if (begin2 != end2) return -1;
    return 0;
  }


  // // This function object class taken from Meyers's "Effective STL", page 38 (first ed.)
  // // It can be used to delete all the elements in an STL container of pointers.
  // struct DeleteObject
  // {
  //   template<typename T> void operator()(const T* ptr) const
  //   {
  //     delete ptr;
  //   }
  // };


} // end of namespace CoCoA

#endif
