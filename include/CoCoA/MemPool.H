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

#endif
