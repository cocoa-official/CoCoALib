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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SmartPtrIRC.H,v 1.15 2022/02/18 14:11:49 abbott Exp $
// $Log: SmartPtrIRC.H,v $
// Revision 1.15  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.14  2021/10/30 17:21:42  abbott
// Summary: Removed cruft (and misleading comment)
//
// Revision 1.13  2021/06/21 11:19:02  abbott
// Summary: Now uses atomic counter -- surely I made this change ages ago!?!
//
// Revision 1.12  2017/04/18 15:52:02  abbott
// Summary: Added missing include for PREPROCESSOR_DEFNS.H
//
// Revision 1.11  2012/05/04 08:32:52  abbott
// Added the threadsafe hack (via several #ifdefs).
//
// Revision 1.10  2010/10/22 13:38:02  bigatti
// -- removed include of old documentation
//
// Revision 1.9  2010/09/30 14:29:54  abbott
// Added memfn so that classes derived from IntrusiveReferenceCount can
// read their own ref count values -- this can be useful for debugging.
//
// Revision 1.8  2010/09/06 09:03:29  abbott
// Added 2 commented-out ctors (which accept pointers to another type <T2>).
//
// Revision 1.7  2010/09/03 15:16:53  abbott
// Added an implementation of assignment for SmartPtrIRC.
//
// Revision 1.6  2010/07/14 08:40:25  bigatti
// -- added IhaveRefCountOne()
//
// Revision 1.5  2009/11/26 16:13:34  bigatti
// -- including stddef instead of stdlib
//
// Revision 1.4  2009/09/24 12:30:44  abbott
// Removed dependency on std::swap which forced inclusion of system header <algorithm>.
// Replaced call to std::swap by explicit code.  Corrected copyright year period.
//
// Revision 1.3  2008/04/10 15:02:43  abbott
// Added mySwap to SmartPtrIRC as it is needed for the new matrices.
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2006/12/06 17:25:43  cocoa
// -- removed #include "config.H"
//
// Revision 1.6  2006/11/20 15:55:03  cocoa
// ring is now a class again.  Improved definitions of operator-> in derived classes.
//
// Revision 1.5  2006/11/11 14:58:28  cocoa
// Documented SmartPtrIRC and IntrusiveReferenceCount, especially explaining
// about the virtual destructor.
//
// Revision 1.4  2006/11/03 16:16:13  cocoa
// changed data member from protected to private
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.3  2006/04/28 16:33:51  cocoa
// Used SmartPtrIRC for PPOrderings.
//
// Revision 1.2  2006/01/20 17:44:20  cocoa
// Friday afternoon check-in: you have been warned!
// Fixed ref counting bug in SmartPtrIRC.
// Added sum/product functions in tmp.H.
// Minor tidying here and there.
//
// Revision 1.1  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//

#endif
