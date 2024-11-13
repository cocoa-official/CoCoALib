#ifndef CoCoA_SmartPtrIRC_H
#define CoCoA_SmartPtrIRC_H

//   Copyright (c)  2006-2009  John Abbott, and Anna M. Bigatti

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


#include "CoCoA/PREPROCESSOR_DEFNS.H"

#include <cstddef>
//using std::size_t;
#include <atomic>

namespace CoCoA
{

  // This template class assumes that T derives from the class IntrusiveReferenceCount (below).
  template <typename T> class SmartPtrIRC
  {
  public:
    explicit SmartPtrIRC(T* RawPtr): myRawPtrValue(RawPtr) { myRawPtrValue->myRefCountInc(); }
//    template<typename T1> explicit SmartPtrIRC(T1* RawPtr); //???
//     template <typename T2>
//       SmartPtrIRC(T2* RawPtr): myRawPtrValue(RawPtr) { myRawPtrValue->myRefCountInc(); }
    SmartPtrIRC(const SmartPtrIRC<T>& ToBeCopied): myRawPtrValue(ToBeCopied.myRawPtrValue) { myRawPtrValue->myRefCountInc(); }
//     template <typename T2>
//       SmartPtrIRC(const SmartPtrIRC<T2>& rhs): myRawPtrValue(rhs.myRawPtrValue) { myRawPtrValue->myRefCountInc(); }
    ~SmartPtrIRC() { myRawPtrValue->myRefCountDec(); }
    void mySwap(SmartPtrIRC& other) { T* SwapTmp = myRawPtrValue; myRawPtrValue = other.myRawPtrValue; other.myRawPtrValue = SwapTmp; }
    SmartPtrIRC<T>& operator=(const SmartPtrIRC<T>& rhs) { rhs.myRawPtrValue->myRefCountInc(); myRawPtrValue->myRefCountDec(); myRawPtrValue = rhs.myRawPtrValue; return *this; };
//     template <typename T2>
//       SmartPtrIRC<T>& operator=(const SmartPtrIRC<T2>& rhs) { rhs.myRawPtrValue->myRefCountInc(); myRawPtrValue->myRefCountDec(); myRawPtrValue = rhs.myRawPtrValue; return *this; }; // NEVER DEFINED -- assignment disabled
  public:
    T* myRawPtr() const { return myRawPtrValue; }    ///< to avoid too many friend stmts
    T* operator->() const { return myRawPtrValue; }  ///< allow const member fns to be called

  private: // Data members
    T* myRawPtrValue;
  };

  // For  ==  and  !=  both args must be pointers to exactly the same type.
  template <typename T>
  inline bool operator==(const SmartPtrIRC<T>& ptr1, const SmartPtrIRC<T>& ptr2)
  {
    return ptr1.myRawPtr() == ptr2.myRawPtr();
  }

  template <typename T>
  inline bool operator!=(const SmartPtrIRC<T>& ptr1, const SmartPtrIRC<T>& ptr2)
  {
    return !(ptr1 == ptr2);
  }


  // This class is intended to be used solely as a (protected) base class.
  class IntrusiveReferenceCount
  {
    // Note: this class has no public members.
  protected:
    IntrusiveReferenceCount();          ///< ctor just sets ref count to 0
    virtual ~IntrusiveReferenceCount(); ///< dtor MUST BE VIRTUAL (see doc)
    void myRefCountInc() const;         ///< increments ref count
    void myRefCountDec() const;         ///< decrements ref count, and self destructs if count becomes 0
    void myRefCountZero() const;        ///< sets myRefCount to 0
    std::size_t myRefCount() const;
  private: // data member
    mutable std::atomic<std::size_t> myRefCountValue;     ///< intrusive reference count value
  };


  inline IntrusiveReferenceCount::IntrusiveReferenceCount():
      myRefCountValue(0)
  {}

  inline IntrusiveReferenceCount::~IntrusiveReferenceCount()
  {}


  inline void IntrusiveReferenceCount::myRefCountInc() const
  {
#if !defined(CoCoA_THREADSAFE_HACK)
    ++myRefCountValue;
#endif
  }

  inline void IntrusiveReferenceCount::myRefCountDec() const
  {
#if !defined(CoCoA_THREADSAFE_HACK)
    --myRefCountValue;
    if (myRefCountValue == 0) delete this;
#endif
  }

  inline void IntrusiveReferenceCount::myRefCountZero() const
  {
    myRefCountValue = 0;
  }


  inline std::size_t IntrusiveReferenceCount::myRefCount() const 
  {
#if !defined(CoCoA_THREADSAFE_HACK)
    return myRefCountValue;
#else
    return 1; // just to keep compiler quiet
#endif
  }

} // end of namespace CoCoA

#endif
