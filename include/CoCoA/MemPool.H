#ifndef CoCoA_MemPool_H
#define CoCoA_MemPool_H

//   Copyright (c)  2005,2006,2010  John Abbott,  Anna M. Bigatti

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
#include <memory>
//using std::unique_ptr;
#include <string>
//using std::string;

namespace CoCoA
{

  class MemPoolFast;  // fwd decl

  std::ostream& MemPoolSetLogStream(std::ostream& out); // default logging stream is clog
  std::ostream& MemPoolSetErrStream(std::ostream& out); // default error stream is cerr


  // STILL EXPERIMENTAL: supposed to help you write exception clean code.
  class AutoPtrSlice
  {
  public:
    AutoPtrSlice(MemPoolFast* MemMgr, void* ptr);
    AutoPtrSlice(const AutoPtrSlice& copy);            // empties copy
    AutoPtrSlice& operator=(const AutoPtrSlice& rhs);  // empties rhs
    ~AutoPtrSlice();
    void* myPtr() const;
    void* myRelease();
  private:
    MemPoolFast* myMemMgr; // deliberately not mutable
    mutable void* mySlicePtr;
  };




  class loaf; // fwd decl, defined in MemPool.C

  /*-----------------------------------------------------------------*/
  /** \include MemPool.txt  */
  /*-----------------------------------------------------------------*/
  // This is the normal non-debugging version of MemPool; contrast with MemPoolDebug below.
  class MemPoolFast
  {
  public:
    enum FillNewLoaf_t {DontFillNewLoaf, FillNewLoaf}; // enum to avoid passing a bool argument in the ctor
  public:
    MemPoolFast(std::size_t sz, const std::string& name="Unnamed-MemPool", FillNewLoaf_t FillFlag = DontFillNewLoaf);
    ~MemPoolFast();
    void* alloc();
    void* alloc(std::size_t);
    void free(void*);
    void free(void*, std::size_t);
    AutoPtrSlice myAlloc(); // change name!!!
    bool IamOriginator(void* ptr) const;
    void SetVerbosityLevel(unsigned int lev);
    void myOutputStatus() const;
    std::size_t mySliceSize() const { return mySliceBytes; }

    typedef void** slice_t;
    static unsigned int ourInitialVerbosityLevel;
  public: // disable copy construction and assignment
    MemPoolFast(const MemPoolFast&) = delete;
    MemPoolFast& operator=(const MemPoolFast&) = delete;


  private: // member functions which are implementation details
    slice_t MakeNewLoaf();
  private: // data members
    const std::size_t mySliceSizeReq; ///< Size in bytes originally requested in ctor.
    const std::string myName;         ///< identification string, handy for debugging.
    const std::size_t mySliceWords;   ///< Size actually used as number of "words" (a word is a slice_t)
    const std::size_t mySliceBytes;   ///< Size actually used; always have mySliceBytes >= mySliceSizeReq.
    const bool myFillNewLoaf;         ///< if set, each new loaf is filled with copies of a known value.
    std::size_t mySlicesPerLoaf;      ///< approx number of slices new loaf should have (increases with each new loaf)
    slice_t myHeadOfFreeList;
    std::unique_ptr<loaf> myLoaves; // intrusive linked list of loaves
    unsigned int myVerbosityLevel; ///< if >=2 then a log message is printed each time a new loaf is made
  };


  inline AutoPtrSlice::AutoPtrSlice(MemPoolFast* MemMgr, void* ptr):
      myMemMgr(MemMgr),
      mySlicePtr(ptr)
  {}


  inline AutoPtrSlice::AutoPtrSlice(const AutoPtrSlice& copy):
      myMemMgr(copy.myMemMgr),
      mySlicePtr(copy.mySlicePtr)
  {
    copy.mySlicePtr = nullptr;
  }


  inline AutoPtrSlice::~AutoPtrSlice()
  {
    if (mySlicePtr) myMemMgr->free(mySlicePtr);
  }


  inline void* AutoPtrSlice::myPtr() const
  {
    return mySlicePtr;
  }


  inline void* AutoPtrSlice::myRelease()
  {
    void* ans = mySlicePtr;
    mySlicePtr = nullptr;
    return ans;
  }


  inline AutoPtrSlice MemPoolFast::myAlloc()
  {
    return AutoPtrSlice(this, alloc());
  }


  //---------------------------------------------------------------------------
  // This is the debugging version of MemPool; contrast with MemPoolFast above.
  class MemPoolDebug
  {
  public:
    MemPoolDebug(std::size_t sz, const std::string& name="Unnamed-MemPool", std::size_t debug_margin=ourDefaultMarginSize);
    // debug_margin is the number of words before and after each ptr
    // Copy ctor and assignment disabled -- see below.
    ~MemPoolDebug();
    void* alloc();
    void* alloc(std::size_t);
    void free(void*);
    void free(void*, std::size_t);

    void InterceptAlloc(std::size_t nth);
    void InterceptFree(std::size_t nth);
    void SetDebugLevel(unsigned int lev);
    void SetVerbosityLevel(unsigned int lev);
    std::size_t mySliceSize() const { return myMemMgr.mySliceSize() - 2*myMarginWords*sizeof(slice_t); } // may be slightly larger than the size requested

    // Globals used to control debugging level and verbosity.
    static unsigned int ourInitialVerbosityLevel;
    static unsigned int ourInitialDebugLevel;
    static unsigned int ourDefaultMarginSize; // size in words: 1 word = sizeof(slice_t)
    static double ourOutputStatusInterval; // print some info every so often
  public: // disable copy ctor and assignment
    MemPoolDebug(const MemPoolDebug&) = delete;
    MemPoolDebug& operator=(const MemPoolDebug&) = delete;

  private: // data members
    typedef MemPoolFast::slice_t slice_t;
    static constexpr long AliveMark = ~0xbadc0c0a;
    long myAliveOrDead;         ///< Heuristic trick to recognise destroyed MemPools when debugging.

    std::string myName;         ///< Identifier useful during debugging.
    std::size_t mySliceSizeReq; ///< Size in bytes originally requested in ctor.
    std::size_t myMarginWords;  ///< number of margin words before and after each block
    std::size_t mySliceWords;   ///< Size of each slice in words; 1 word = sizeof(slice_t).
    std::size_t mySliceBytes;   ///< Size actually used; always have mySliceBytes >= mySliceSizeReq + 2*myMarginWords*sizeof(slice_t)
    MemPoolFast myMemMgr;       ///< myMemMgr actually manages the alloc/free requests
    slice_t myHeadOfUsedList;   ///< Used only in the highest debug level

    unsigned int myDebugLevel;  ///< debug level should be one of 0, 1, or 2 (>2 is same as 2)
    unsigned int myVerbosityLevel; ///< verbosity level should be one of 0, 1, 2 or 3 (>3 is same as 3)
    std::size_t myAllocCount;
    std::size_t myAllocWatchPoint;
    std::size_t myFreeCount;
    std::size_t myFreeWatchPoint;
    std::size_t myInUseCount;
    std::size_t myInUseMax;
    double myNextOutputStatusTime;
  private: // implementation details.
    static std::size_t ourCheckCtorSizeArg(std::size_t sz);
    void myAliveCheck();

    void AllocMark(slice_t p);
    void FreeMark(slice_t p);

    void FullOverwriteFreeCheck() const;
    void OverwriteErrorMesg(slice_t, void*) const;
    void OverwriteFreeCheck(slice_t) const;
    bool FreeError(void* ptr) const;
    void FreeErrorMesg(void* ptr, const std::string& reason) const;
    bool AlreadyFreed(void* ptr) const;
    bool WrittenOutOfRange(void* ptr) const;

    void PrintSlice(void* ptr) const;

    void DtorErrorMesg() const;
    void PrintStats() const;
    void AllocMesg(void* ptr) const;
    void AllocWrongSizeMesg(std::size_t sz, void* ptr) const;
    void AllocLargeSliceMesg(void* ptr) const;
    void FreeMesg(void* ptr) const;
    void FreeWrongSizeMesg(std::size_t sz, void* ptr) const;
    void FreeLargeSliceMesg(void* ptr) const;
    void FreeZeroPtrMesg() const;
    void intercepted() const;
  };


  //------------------------------------------------------------
  // This is the fake version of MemPool -- it forwards everything to operator new
  // This is threadsafe if operator new is threadsafe.
  class MemPoolFake
  {
  public:
    MemPoolFake(std::size_t sz, const std::string& name="Unnamed-MemPool");
    ~MemPoolFake();
    void* alloc();
    void* alloc(std::size_t);
    void free(void*);
    void free(void*, std::size_t);
    void myOutputStatus() const;
    std::size_t mySliceSize() const { return mySliceSizeReq; }

  public: // disable copy construction and assignment
    MemPoolFake(const MemPoolFake&) = delete;
    MemPoolFake& operator=(const MemPoolFake&) = delete;

  private: // data members
    const std::size_t mySliceSizeReq; ///< Size in bytes originally requested in ctor.
    const std::string myName;         ///< identification string, handy for debugging.
  };



/////////////////////////////////////////////////////////////////////////////
// Define the type MemPool according to the setting of CoCoA_MEMPOOL_DEBUG & CoCoA_MEMPOOL_DISABLE.

  typedef int junk;
  
///#if defined!CoCoA_MEMPOOL_DISABLE) || !defined(CoCoA_ALLOW_NON_THREASDSAFE)
#if defined(CoCoA_MEMPOOL_DISABLE)
  typedef MemPoolFake MemPool;
#elif defined(CoCoA_MEMPOOL_DEBUG)
  typedef MemPoolDebug MemPool;
#else
  typedef MemPoolFast MemPool;
#endif

} // end of namespace CoCoA



// RCS Header/Log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/MemPool.H,v 1.14 2022/02/18 14:11:46 abbott Exp $
// $Log: MemPool.H,v $
// Revision 1.14  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.13  2021/10/30 17:07:23  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.12  2021/02/11 23:03:52  abbott
// Summary: Constexpr instead of just const
//
// Revision 1.11  2021/01/07 14:17:10  abbott
// Summary: Corrected copyright
//
// Revision 1.10  2020/12/04 10:13:21  abbott
// Summary: Updated: big slices now handled by sys mem mgr (redmine 1522)
//
// Revision 1.9  2020/10/27 09:51:14  abbott
// Summary: Added flag to disable MemPools (now indep of threadsafe flag)
//
// Revision 1.8  2019/03/18 17:56:49  abbott
// Summary: Changed 0 into nullptr where appropriate
//
// Revision 1.7  2019/03/04 10:18:39  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.6  2017/04/18 15:52:02  abbott
// Summary: Added missing include for PREPROCESSOR_DEFNS.H
//
// Revision 1.5  2012/04/23 10:29:25  abbott
// Added MemPoolFake for use when CoCoA_THREADSAFE_HACK is set.
//
// Revision 1.4  2010/12/20 11:47:15  abbott
// Updated MemPool so user can specify on which streams to print
// logging info and errors (when in verbose mode).
//
// Revision 1.3  2010/10/22 13:35:33  abbott
// Added mySliceSize member fn.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.3  2006/12/06 17:14:12  cocoa
// -- removed #include "config.H"
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
// Revision 1.7  2006/05/15 16:28:49  cocoa
// Fixed bug which appeared only when there was a double free.
//
// Revision 1.6  2006/04/27 13:54:19  cocoa
// Improved MemPools with an experimental way of handling raw memory in an
// exception clean manner.
//
// Revision 1.5  2006/04/11 16:34:29  cocoa
// Completed overhaul of MemPool (incl documentation).
// Modified GMPAllocator so that you can specify the slice size
// in the ctor -- useful for experimentation.
//
// Revision 1.4  2006/04/07 16:44:52  cocoa
// Considerably updated MemPool design -- it works, and I'm about to test
// its efficiency against the old one.
//
// Revision 1.3  2006/03/31 13:41:00  cocoa
// Temporary check in with partially updated MemPool code.
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.5  2004/11/19 15:14:09  cocoa
// (a) Added new check to MemPool so that it can signal an
//     error if one tries to use a MemPool after it has been
//     destroyed.
// (b) Improved makefile in TEST/ so that it checks output,
//     and prints useful messages if the test fails.
// (c) Tidied ring.txt a bit (still more to do).
//
// Revision 1.4  2004/11/08 14:03:32  cocoa
// -- small changes for doxygen support
//
// Revision 1.3  2004/11/04 18:49:03  cocoa
// (1) MemPool code cleaned up; now adheres better to the coding
//     conventions.
// (2) Log messages printed by MemPool (in debugging mode) are
//     more uniform.
//
// Revision 1.2  2004/01/28 15:17:58  cocoa
// Minor tidying: better indentation, cleaner code for non-debugging case.
// [backed off an abortive attempt to protect against mixing code compiled
// with and without CoCoA_MEMPOOL_DEBUG]
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.8  2003/06/23 16:19:47  abbott
// Minor cleaning prior to public release.
//
// Revision 1.7  2003/05/08 08:54:28  abbott
// Minor cleaning: replaced size_t by std::size_t.
// (apparently these changes were made in Nov 2002, but never checked in).
//
// Revision 1.6  2002/10/25 12:57:08  abbott
// Added new alloc/free functions which do not require a size parameter;
// the size defaults to the "natural" size for the MemPool.
//
// Removed the member MEMPOOL_ALLOCATED_WORD_INSIDE from the debugging
// version since now different values are used for different slices.
//
// Revision 1.5  2002/03/05 17:30:42  abbott
// Added two new debugging functions: InterceptAlloc and InterceptFree
// which allow one easily to intercept the nth call to alloc/free
// using a debugger (see details in the documentation).  To implement
// this I added two new data fields to MemPool (in the debug case).
//
// Added a new field containing the number of loaves created by each
// MemPool object: this is used to determine the growth of the size
// of each loaf (doubles every second time).
//
// Revision 1.4  2002/02/07 11:47:43  abbott
// Corrected two member function declarations, and added another which
// was missing.
//
// Revision 1.3  2001/11/28 21:31:17  abbott
// Radical reorganization of the code:
//  * the separation into two types of file MemPool* and MemPool_DEBUG*
//    has been eliminated; everything is now in MemPool.*
//  * compatible with gcc-3.0.2 regarding use of namespace std;
//  * many names changed (in accordance with the new coding conventions,
//    and for better comprehensibility);
//  * new more powerful checks added (e.g. write to a freed block),
//    and more levels of debugging supported;
//  * improved error messages, and progress messages.
//
// Revision 1.2  2001/11/16 19:25:02  bigatti
// added:  std::
// for compatibility with gcc-3
//
// Revision 1.1  2001/10/04 14:54:28  abbott
// Initial revision
//

#endif
