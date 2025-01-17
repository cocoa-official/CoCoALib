#ifndef CoCoA_SparsePolyIter_H
#define CoCoA_SparsePolyIter_H

//   Copyright (c)  2005-2007,2009,2012  John Abbott, and Anna M. Bigatti

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


// Header file for the abstract class SparsePolyRingBase etc.

#include "CoCoA/PPMonoid.H"
#include "CoCoA/ring.H"

#include <iosfwd>
#include <memory>
//using std::uniqe_ptr;

namespace CoCoA
{

  class symbol;         // fwd decl -- defined in symbol.H
  class matrix;         // fwd decl -- defined in matrix.H

  class SparsePolyIter; // fwd decl -- defined later in this file

  class SparsePolyIterBase
  {
  public:
    virtual ~SparsePolyIterBase() {};
    virtual SparsePolyIterBase* myClone() const = 0;

    virtual bool IamEnded() const = 0;
    virtual void myForward() = 0;
//    virtual SparsePolyIter& operator++() = 0;
//    virtual SparsePolyIter operator++(int) = 0;
    virtual RingElemAlias myCoeff() const = 0;
    virtual ConstRefPPMonoidElem myPP() const = 0;
//    virtual bool operator==(const iter& i2) const = 0;
//    virtual bool operator!=(const iter& i2) const = 0;
  protected:
    virtual const void* myPolyPtr() const = 0;
    virtual const void* myTermPtr() const = 0;
    friend class SparsePolyIter;
    friend bool operator==(const SparsePolyIter& i1, const SparsePolyIter& i2);
  };

  // This is a "smart pointer class" -- it owns the object pointed to.
  class SparsePolyIter
  {
  public:
    SparsePolyIter(SparsePolyIterBase* ptr);
    SparsePolyIter(const SparsePolyIter& copy);
    // Default dtor works OK
    SparsePolyIter& operator=(const SparsePolyIter& rhs);
    SparsePolyIterBase* operator->();
    const SparsePolyIterBase* operator->() const;
    SparsePolyIter& operator++();
    SparsePolyIter operator++(int);
  private:
    std::unique_ptr<SparsePolyIterBase> myIterPtr;
    friend bool IsEnded(const SparsePolyIter& i);
    friend RingElemAlias coeff(const SparsePolyIter& i);
    friend ConstRefPPMonoidElem PP(const SparsePolyIter& i);
    friend bool operator==(const SparsePolyIter& i1, const SparsePolyIter& i2);
  };



  //--- iterator -------------------------------------
  inline SparsePolyIter::SparsePolyIter(SparsePolyIterBase* ptr):
      myIterPtr(ptr)
  {}


  inline SparsePolyIter::SparsePolyIter(const SparsePolyIter& copy):
      myIterPtr(copy.myIterPtr->myClone())
  {}


  inline SparsePolyIter& SparsePolyIter::operator=(const SparsePolyIter& rhs)
  {
//    if (&copy == this.get()) return *this;
    myIterPtr.reset(rhs.myIterPtr->myClone());
    return *this;
  }


  inline SparsePolyIterBase* SparsePolyIter::operator->()
  { return myIterPtr.get(); }


  inline const SparsePolyIterBase* SparsePolyIter::operator->() const
  { return myIterPtr.get(); }


  inline bool IsEnded(const SparsePolyIter& i)
  { return i.myIterPtr->IamEnded(); }


  inline SparsePolyIter& SparsePolyIter::operator++()
  {
    if (myIterPtr->IamEnded())
      CoCoA_THROW_ERROR1(ERR::PolyIterEnded);
    myIterPtr->myForward();
    return *this;
  }

  inline SparsePolyIter SparsePolyIter::operator++(int)
  {
    if (myIterPtr->IamEnded())
      CoCoA_THROW_ERROR1(ERR::PolyIterEnded);
    SparsePolyIter ans(*this);
    myIterPtr->myForward();
    return ans;
  }


  inline RingElemAlias coeff(const SparsePolyIter& i)
  {
    if (IsEnded(i))
      CoCoA_THROW_ERROR1(ERR::PolyIterEnded);
    return i.myIterPtr->myCoeff();
  }


  inline ConstRefPPMonoidElem PP(const SparsePolyIter& i)
  {
    if (IsEnded(i))
      CoCoA_THROW_ERROR1(ERR::PolyIterEnded);
    return i.myIterPtr->myPP();
  }


  inline bool operator==(const SparsePolyIter& i1, const SparsePolyIter& i2)
  {
    if (i1.myIterPtr->myPolyPtr() != i2.myIterPtr->myPolyPtr())
      CoCoA_THROW_ERROR1(ERR::MixedPolyIters);
    return i1.myIterPtr->myTermPtr() == i2.myIterPtr->myTermPtr();
  }


  inline bool operator!=(const SparsePolyIter& i1, const SparsePolyIter& i2)
  { return !(i1==i2); }


} // end of namespace CoCoA

#endif
