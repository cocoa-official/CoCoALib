//   Copyright (c)  2005,2006,2010  John Abbott and Anna M. Bigatti

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


#include "CoCoA/MemPool.H"
#include "CoCoA/time.H"
#include "CoCoA/error.H"
#include "CoCoA/assert.H"

//#include <string>
using std::string;
//#include <cstddef>
using std::size_t;
#include <cstdint>
using std::uintptr_t;

#include <algorithm>
using std::fill;
#include <cmath>
using std::ceil;
#include <cstddef>
using std::ptrdiff_t;
#include <iomanip>
using std::setw;
#include <iostream>
using std::ostream;
using std::endl;
//#include <memory>
using std::unique_ptr;

namespace CoCoA
{

  // Code for managing the ostream used for logging and for reporting errors
  namespace // anonymous
  {
    ostream* LogStreamPtr = &std::clog; // default logging stream

    inline ostream& LogStream()
    {
      return *LogStreamPtr;
    }

    ostream* ErrStreamPtr = &std::cerr; // default error stream

    inline ostream& ErrStream()
    {
      return *ErrStreamPtr;
    }

  } // end of anonymous namespace

  ostream& MemPoolSetLogStream(ostream& out)
  {
    ostream& ans = *LogStreamPtr;
    LogStreamPtr = &out;
    return ans;
  }

  ostream& MemPoolSetErrStream(ostream& out)
  {
    ostream& ans = *ErrStreamPtr;
    ErrStreamPtr = &out;
    return ans;
  }

  //----------------------------------------------------------------------

  // If you call this fn, you don't deserve to have it inline >-}
  AutoPtrSlice& AutoPtrSlice::operator=(const AutoPtrSlice& rhs)
  {
    if (myMemMgr != rhs.myMemMgr) CoCoA_THROW_ERROR("MemPools are different", "AutoPtrSlice assignment");
    if (mySlicePtr == rhs.mySlicePtr) return *this; // ??? or give error???
    if (mySlicePtr) myMemMgr->free(mySlicePtr);
    mySlicePtr = rhs.mySlicePtr;
    rhs.mySlicePtr = nullptr;
    return *this;
  }


  typedef MemPoolFast::slice_t slice_t; // bring typedef of slice_t out to this level to save typing

  namespace // anonymous namespace for file local values, functions etc.
  {

    // Some constants -- implementation details.
    constexpr size_t MaxSliceSize = 256;  // slices larger than this are handled directly by ::operator new (delete)

    constexpr size_t WordSize = sizeof(slice_t);
    constexpr size_t allowance = 16;              // Aim to have loafsize+allowance just less than a power of 2.
    constexpr size_t MinLoafBytes = 1024;    // First loaf is about this size; later ones are bigger.
    constexpr size_t MaxLoafBytes = 64*1024; // We do not make any single loaf bigger than this.

    // Constants for MemPoolDebug (values are pretty arbitrary)
    const slice_t MEMPOOL_ALLOCATED_WORD_BEFORE = reinterpret_cast<slice_t>(13ul<<8);
    const slice_t MEMPOOL_ALLOCATED_WORD_AFTER = reinterpret_cast<slice_t>(37ul<<8);
    const slice_t MEMPOOL_FREE_WORD = reinterpret_cast<slice_t>(99ul<<8);


    // Initial loaf has size 2^k*MinLoafBytes where k is smallest such
    // that we can fit at least NumSlices into the loaf.
    size_t InitialLoafSlices(size_t sz, size_t NumSlices)
    {
      if (sz > MaxLoafBytes/2) return 2;    // ludicrously large slice

      size_t TargetLoafSize = MinLoafBytes;
      // NB Below we don't compute sz*NumSlices as it may overflow.
      while (TargetLoafSize <= MaxLoafBytes/2 && (TargetLoafSize-allowance)/sz < NumSlices)
        TargetLoafSize *= 2;
      return (TargetLoafSize-allowance)/sz;
    }


    // Double size up to a limit of MaxLoafBytes.
    // Choose number of slices so that loaf size is roughly a power of 2 times MinLoafBytes
    size_t IncrLoafSlices(size_t SlicesPerLoaf, size_t SliceBytes)
    {
      const size_t CurrLoafBytes = SlicesPerLoaf*SliceBytes; // cannot overflow.
      if (CurrLoafBytes > MaxLoafBytes/2) return SlicesPerLoaf; // do not increase further
      size_t TargetLoafSize = MinLoafBytes;
      while (TargetLoafSize <= CurrLoafBytes)
        TargetLoafSize *= 2;
      TargetLoafSize *= 2;
      return (TargetLoafSize-allowance)/SliceBytes;
    }

  } // end of anonymous namespace

  // Raw memory comes in "loaves" which are composed of "slices" of the desired size.
  class loaf
  {
  public:
    loaf(size_t NumSlices, size_t SliceBytes, bool FillBeforeSlicing);
    ~loaf();
    void myAppend(loaf* NewLoaf);
    slice_t myFirstSlice() const;
    bool IamOriginator(void* ptr) const;
    void myFreeCounterReset();
    void myCountFreeSlice(void* ptr);
    void myOutputStatus() const;
  public: // disable copy ctor and assignment
    loaf(const loaf&) = delete;
    loaf& operator=(const loaf&) = delete;
  private: // data members of loaf
    unique_ptr<loaf> myNext;
    const size_t mySliceWords;
    const size_t myNumSlices;
    slice_t* const myBegin;
    slice_t* const myEnd;
    size_t myFreeCounter; // only used when printing stats
    bool myRangeCheck(void* ptr) const;
    bool myAlignmentCheck(void* ptr) const;
  };


  /////////////////////////////////////////////////////////////////////////////
  // ---------------------- loaf functions ----------------------

  loaf::loaf(size_t NumSlices, size_t SliceBytes, bool FillBeforeSlicing):
      myNext(),
      mySliceWords(SliceBytes/WordSize),
      myNumSlices(NumSlices),
      myBegin(static_cast<slice_t*>(::operator new(NumSlices*SliceBytes))),
      myEnd(myBegin + NumSlices*mySliceWords),
      myFreeCounter(0)
  {
    CoCoA_ASSERT(SliceBytes%WordSize == 0);
    // First word in each slice must point to the start of the next slice.  Last ptr is 0.

    if (FillBeforeSlicing) // cond true if MemPoolFast was created by a MemPoolDebug
      fill(myBegin, myEnd, MEMPOOL_FREE_WORD);

    // The loop below slices up the newly allocated loaf into a linked list of free slices;
    // each slice (except the last) points to the slice immediately following.
    slice_t* next = nullptr;
    for (slice_t* curr = myBegin + (myNumSlices-1)*mySliceWords; curr >= myBegin; curr -= mySliceWords)
    {
      *curr = reinterpret_cast<slice_t>(next);
      next = curr;
    }
  }


  loaf::~loaf()
  { ::delete myBegin; }


  void loaf::myAppend(loaf* NewLoaf)
  {
    CoCoA_ASSERT(myNext.get() == nullptr);
    myNext.reset(NewLoaf); // assumes ownership!
  }


  slice_t loaf::myFirstSlice() const
  {
    return reinterpret_cast<slice_t>(myBegin);
  }


  inline bool loaf::myRangeCheck(void* ptr) const
  {
    return (ptr >= myBegin && ptr < myEnd);
  }

  inline bool loaf::myAlignmentCheck(void* ptr) const
  {
    // ASSUME ptr is in my range...
    CoCoA_ASSERT(myRangeCheck(ptr));

    // ...now check whether it is correctly aligned.
    const size_t SliceBytes = mySliceWords*WordSize;
    ptrdiff_t d = reinterpret_cast<char*>(ptr) - reinterpret_cast<char*>(myBegin);
    return (d%SliceBytes == 0);
  }


  // Check whether ptr belongs to some loaf in the chain.
  bool loaf::IamOriginator(void* ptr) const
  {
    if (myRangeCheck(ptr))
      return myAlignmentCheck(ptr);

    if (myNext.get()) return myNext->IamOriginator(ptr);
    return false;
  }


  // The next three functions are closely related; any ideas for a better impl?
  void loaf::myFreeCounterReset()
  {
    myFreeCounter = 0;
    if (myNext.get()) myNext->myFreeCounterReset();
  }

  void loaf::myCountFreeSlice(void* ptr)
  {
    if (myRangeCheck(ptr))
    {
      ++myFreeCounter;
      return;
    }
    CoCoA_ASSERT(myNext.get() != nullptr);
    myNext->myCountFreeSlice(ptr);
  }

  void loaf::myOutputStatus() const
  {
    if (myNext.get()) myNext->myOutputStatus();
    const double full = 1-double(myFreeCounter)/double(myNumSlices);
    void* LastByte = reinterpret_cast<char*>(myEnd)-1;
    LogStream() << "[Log] loaf=[" << myBegin << "--" << LastByte << "]\t  slices=" << myNumSlices << "\t  full=" << full << endl;
  }



  /////////////////////////////////////////////////////////////////////////////
  // Implementations for MemPoolFast


  // non-constant static member initializations
  unsigned int MemPoolFast::ourInitialVerbosityLevel = 0;

  //-------------------- constructor & destructor --------------------//

  MemPoolFast::MemPoolFast(size_t sz, const string& name, FillNewLoaf_t FillFlag):
      mySliceSizeReq(sz),
      myName(name),
      mySliceWords(1+(sz-1)/WordSize),
      mySliceBytes(WordSize*mySliceWords),
      myFillNewLoaf(FillFlag == FillNewLoaf)
  {
    if (sz == 0) CoCoA_THROW_ERROR(ERR::MemPoolZero, "MemPoolFast ctor");

    mySlicesPerLoaf = InitialLoafSlices(mySliceBytes, 16); // first loaf should have at least 16 slices
    myHeadOfFreeList = nullptr;
    myVerbosityLevel = ourInitialVerbosityLevel;
  }


  MemPoolFast::~MemPoolFast()
  {}


//-------------------- alloc & free --------------------//

// No overall benefit was observed from making this function inline.
  void* MemPoolFast::alloc()
  {
    if (mySliceSizeReq > MaxSliceSize) return ::operator new(mySliceSizeReq);
    if (myHeadOfFreeList == nullptr)
      myHeadOfFreeList = MakeNewLoaf();
    slice_t p = myHeadOfFreeList;
    myHeadOfFreeList = reinterpret_cast<slice_t>(*p);

    return p;
  }


  void* MemPoolFast::alloc(size_t sz)
  {
    if (sz != mySliceSizeReq) return ::operator new(sz);
    return alloc();
  }


// No overall benefit was observed from making this function inline.
  void MemPoolFast::free(void* ptr)
  {
    if (ptr == nullptr) return;
    if (mySliceSizeReq > MaxSliceSize) { ::operator delete(ptr); return; }
    slice_t old_head = myHeadOfFreeList;
    myHeadOfFreeList = static_cast<slice_t>(ptr);
    *myHeadOfFreeList = old_head;
  }


  void MemPoolFast::free(void* ptr, size_t sz)
  {
    if (sz != mySliceSizeReq)  { ::operator delete(ptr);  return; }
    free(ptr);
  }


  bool MemPoolFast::IamOriginator(void* ptr) const
  {
    return myLoaves->IamOriginator(ptr);
  }


  void MemPoolFast::SetVerbosityLevel(unsigned int lev)
  {
    myVerbosityLevel = lev;
  }


  void MemPoolFast::myOutputStatus() const
  {
    LogStream() << "[Log] -------------------------------------------------------" << endl;
    LogStream() << "[Log] Status for MemPool(\"" << myName << "\")   CpuTime=" << CpuTime() << endl;
    if (mySliceSizeReq > MaxSliceSize)
      LogStream() << "[Log] --- LargeSlice(" << mySliceSizeReq << " > " << MaxSliceSize << "): all calls forwarded to system mem mgr" << endl;
    if (!myLoaves.get())
      LogStream() << "[Log] --- This MemPool created no loaves ---" << endl;
    else
    {
      myLoaves->myFreeCounterReset(); // Clear free slice counters in each loaf.
      // For each slice in the free list increment the counter of the loaf owning that slice...
      for (slice_t ptr=myHeadOfFreeList; ptr != nullptr; ptr = *reinterpret_cast<slice_t*>(ptr))
        myLoaves->myCountFreeSlice(ptr);
      myLoaves->myOutputStatus(); // each loaf prints the proportion of its slices still in use.
    }
    LogStream() << "[Log] -------------------------------------------------------" << endl;
  }


  //-------------------- private functions --------------------//

  slice_t MemPoolFast::MakeNewLoaf()
  {
    if (myLoaves.get() != nullptr)
      mySlicesPerLoaf = IncrLoafSlices(mySlicesPerLoaf, mySliceBytes);

    // Create a new loaf and insert it at the front of the list of loaves.
    // Use of raw ptr here is safe as we cannot throw exception once loaf has been built.
    loaf* NewLoaf = new loaf(mySlicesPerLoaf, mySliceBytes, myFillNewLoaf);
    NewLoaf->myAppend(myLoaves.release());
    myLoaves.reset(NewLoaf);

    if (myVerbosityLevel > 1)
    {
      LogStream() << "[Log]"
                  << " MemPoolName=\"" << myName << "\""
                  << " fn=MakeNewLoaf"
                  << " LoafBytes=" << mySlicesPerLoaf*mySliceBytes
                  << " SliceBytes=" << mySliceBytes
                  << " LoafSlices=" << mySlicesPerLoaf
                  << " cputime=" << CpuTime()
                  << endl;
    }
    return myLoaves->myFirstSlice();
  }



  /////////////////////////////////////////////////////////////////////////////
  // Implementation for MemPoolDebug



  //-------------------- stats & debug functions --------------------//

  // non-constant static member initializations
  unsigned int MemPoolDebug::ourInitialVerbosityLevel = 0;
  unsigned int MemPoolDebug::ourInitialDebugLevel = 0;
  unsigned int MemPoolDebug::ourDefaultMarginSize = 4;
  double MemPoolDebug::ourOutputStatusInterval = 1.0E16; // print interval in seconds; default is "almost never".

  //---------------------Functions to fill freed/allocated slices---------//

  void MemPoolDebug::AllocMark(slice_t p)
  {
    // We shall fill the visible part of the slice with a value which depends on p.
    const uintptr_t word_inside = ~(reinterpret_cast<uintptr_t>(p));
    const slice_t MEMPOOL_ALLOCATED_WORD_INSIDE = reinterpret_cast<slice_t>(word_inside);

    size_t i;
    for (i = 0; i < myMarginWords; ++i)
      p[i] = MEMPOOL_ALLOCATED_WORD_BEFORE;
    for ( ; i < mySliceWords - myMarginWords; ++i)
      p[i] = MEMPOOL_ALLOCATED_WORD_INSIDE;
    for ( ; i < mySliceWords; ++i)
      p[i] = MEMPOOL_ALLOCATED_WORD_AFTER;
  }


  void MemPoolDebug::FreeMark(slice_t p)
  {
    // Skip first word as it is used as a "next' pointer.
    for (size_t i = 1; i < mySliceWords ; ++i)
      p[i] = MEMPOOL_FREE_WORD;
  }


  //-----------------------Integrity test functions---------------//

  void MemPoolDebug::FullOverwriteFreeCheck() const
  {
    for (slice_t ptr = myHeadOfUsedList; ptr != nullptr; ptr = reinterpret_cast<slice_t>(*ptr))
      OverwriteFreeCheck(ptr);
  }


  void MemPoolDebug::OverwriteFreeCheck(slice_t p) const
  {
    // Do not check index 0 as it is used for a "next" pointer.
    for (size_t i = 1; i < mySliceWords; ++i)
      if (p[i] != MEMPOOL_FREE_WORD)
      {
        OverwriteErrorMesg(p+myMarginWords, &p[i]);
        p[i] = MEMPOOL_FREE_WORD;  // reset to expected value to avoid repeated error mesgs.
      }
  }


  void MemPoolDebug::OverwriteErrorMesg(slice_t slice_addr, void* overwritten_addr) const
  {
    ErrStream() << "[ERR] MemPoolDebug(\"" << myName
                << "\") ERROR: OVERWRITTEN freed slice, slice_addr=" << slice_addr
                << ", overwritten at addr=" << overwritten_addr << endl;
  }


  //-------------------------------------------------------------------------//

  void MemPoolDebug::FreeErrorMesg(void* ptr, const string& reason) const
  {
    ErrStream() << "[ERR] MemPoolDebug(\"" << myName << "\") ERROR:  free, seq=" << setw(4) << myFreeCount
                << ", addr=" << ptr  << ": " << reason << endl;
  }


  // ptr should be a value returned by alloc (rather than the start
  // of the first margin).
  bool MemPoolDebug::FreeError(void* ptr) const
  {
    if (!myMemMgr.IamOriginator(static_cast<slice_t>(ptr)-myMarginWords))
    {
      FreeErrorMesg(ptr, "not allocated by me (or misaligned)");
      return true;
    }

    if (AlreadyFreed(ptr))
    {
      FreeErrorMesg(ptr, "already freed");
      return true;
    }

    if (WrittenOutOfRange(ptr))
    {
      FreeErrorMesg(ptr, "written out of range");
      PrintSlice(ptr);
      return true;
    }
    return false;
  }


  // Heuristic test to see whether a slice has already been freed;
  // heuristic fails if the margins have been overwritten after being freed.
  // ptr should be a value returned by alloc (rather than the start
  // of the first margin).
  bool MemPoolDebug::AlreadyFreed(void* ptr) const
  {
    if (myMarginWords == 0) return false; // disable test if myMarginWords==0

    slice_t p = static_cast<slice_t>(ptr) - myMarginWords;

    // Check only the margins, in case the data area has been overwritten.
    for (size_t i = 1; i < myMarginWords; ++i)
      if (p[i] != MEMPOOL_FREE_WORD) return false;
    for (size_t i = mySliceWords - myMarginWords; i < mySliceWords; ++i)
      if (p[i] != MEMPOOL_FREE_WORD) return false;

    // Pretty sure the block has already been freed, check for overwriting.
    OverwriteFreeCheck(p);
    return true;
  }


  // Check that the margins (if any) are uncorrupted.
  // ptr should be a value returned by alloc (rather than the start
  // of the first margin).
  bool MemPoolDebug::WrittenOutOfRange(void* ptr) const
  {
    slice_t p = static_cast<slice_t>(ptr) - myMarginWords;

    for (size_t i = 0; i < myMarginWords; ++i)
      if (p[i] != MEMPOOL_ALLOCATED_WORD_BEFORE) return true;
    for (size_t i = mySliceWords - myMarginWords; i < mySliceWords; ++i)
      if (p[i] != MEMPOOL_ALLOCATED_WORD_AFTER) return true;
    return false;
  }


  //------------------------Error message functions---------------//


  void MemPoolDebug::PrintSlice(void* ptr) const
  {
    slice_t p = static_cast<slice_t>(ptr) - myMarginWords;
    size_t i=0;
    ErrStream() << "[ERR] Nature        addr  :    value      (decimal)" << endl;
    ErrStream() << "[ERR] =============================================" << endl;
    for (; i < myMarginWords; ++i)
    {
      ErrStream() << "[ERR] margin  " << setw(12) << static_cast<void*>(&p[i]) << ":"
                  << setw(12) << static_cast<void*>(p[i])
                  << setw(12) << reinterpret_cast<uintptr_t>(p[i]);
      if (p[i] != MEMPOOL_ALLOCATED_WORD_BEFORE)
        ErrStream() << " <--- WRONG!  Should have been "
                    << static_cast<void*>(MEMPOOL_ALLOCATED_WORD_BEFORE);
      ErrStream() << endl;
    }
    ErrStream() << "[ERR] ---------------------------------------------" << endl;
    for (; i < mySliceWords - myMarginWords; ++i)
    {
      ErrStream() << "[ERR] data    " << setw(12) << static_cast<void*>(&p[i]) << ":"
                  << setw(12) << static_cast<void*>(p[i])
                  << setw(12) << reinterpret_cast<uintptr_t>(p[i]) << endl;
    }
    ErrStream() << "[ERR] ---------------------------------------------" << endl;
    for (i = mySliceWords - myMarginWords; i < mySliceWords; ++i)
    {
      ErrStream() << "[ERR] margin  " << setw(12) << static_cast<void*>(&p[i]) << ":"
                  << setw(12) << static_cast<void*>(p[i])
                  << setw(12) << reinterpret_cast<uintptr_t>(p[i]);
      if (p[i] != MEMPOOL_ALLOCATED_WORD_AFTER)
        ErrStream() << " <--- WRONG!  Should have been "
                    << static_cast<void*>(MEMPOOL_ALLOCATED_WORD_AFTER);
      ErrStream() << endl;
    }
    ErrStream() << "[ERR] =============================================" << endl;
    ErrStream() << endl;
  }


  //--------------------------- Message functions ----------------------------//

  void MemPoolDebug::DtorErrorMesg(void) const
  {
    ErrStream() << "[ERR] MemPoolDebug(\"" << myName << "\") ERROR:  dtor"
                << ", unfreed slices: NumSlices=" << myInUseCount << endl;
  }


  void MemPoolDebug::PrintStats() const
  {
    LogStream() << "[Log] MemPoolDebug(\"" << myName << "\"):   STATS:  "
                << "  NumAlloc=" << myAllocCount
                << "  NumFree=" << myFreeCount;
    if (mySliceSizeReq <= MaxSliceSize)
    {
      LogStream() << "  NumUsed=" << myInUseCount
                  << "  MaxUsed=" << myInUseMax << endl;
    }
    else
    {
      LogStream() << "  LargeSlice(" << mySliceSizeReq << " > " << MaxSliceSize << ")" << endl;
    }
    myMemMgr.myOutputStatus();
  }


  void MemPoolDebug::AllocMesg(void* ptr) const
  {
    LogStream() << "[Log] MemPoolDebug(\"" << myName << "\"):   ALLOC " << setw(10) << ptr
                << ", seq=" << setw(4) << myAllocCount
                << ", #alloc=" << myAllocCount << ", #free=" << myFreeCount
                << ", in_use=" << myInUseCount << ", in_use_max=" << myInUseMax << endl;
  }


  void MemPoolDebug::AllocWrongSizeMesg(size_t sz, void* ptr) const
  {
    LogStream() << "[Log] MemPoolDebug(\"" << myName << "\") WARNING: ALLOC " << setw(10) << ptr
                << ", seq=" << setw(4) << myAllocCount
                << "WRONG SIZE sz=" << sz << ", but mySliceSize=" << mySliceSizeReq << endl;
  }

  void MemPoolDebug::AllocLargeSliceMesg(void* ptr) const
  {
    LogStream() << "[Log] MemPoolDebug(\"" << myName << "\")  LARGE SLICE ALLOC " << setw(10) << ptr
                << ", seq=" << setw(4) << myAllocCount << endl;
  }


  void MemPoolDebug::FreeMesg(void* ptr) const
  {
    LogStream() << "[Log] MemPoolDebug(\"" << myName << "\"):   FREED " << setw(10) << ptr
                << ", seq=" << setw(4) << myFreeCount
                << ", #alloc=" << myAllocCount << ", #free=" << myFreeCount
                << ", in_use=" << myInUseCount << ", in_use_max=" << myInUseMax << endl;
  }


  void MemPoolDebug::FreeWrongSizeMesg(size_t sz, void* ptr) const
  {
    LogStream() << "[Log] MemPoolDebug(\"" << myName << "\") WARNING: FREED " << setw(12) << ptr
                << ", seq=" << setw(4) << myFreeCount
                << "WRONG SIZE sz=" << sz << ", but mySliceSize=" << mySliceSizeReq << endl;
  }


  void MemPoolDebug::FreeLargeSliceMesg(void* ptr) const
  {
    LogStream() << "[Log] MemPoolDebug(\"" << myName << "\")  LARGE SLICE FREED " << setw(12) << ptr
                << ", seq=" << setw(4) << myFreeCount << endl;
  }


  void MemPoolDebug::FreeZeroPtrMesg() const
  {
    LogStream() << "[Log] MemPoolDebug(\"" << myName << "\") WARNING: FREED " << 0/*nullptr*/
                << ", seq=" << setw(4) << myFreeCount
                << " ZERO PTR" << endl;
  }


  void MemPoolDebug::intercepted() const
  {
    LogStream() << "[Log] MemPoolDebug(\"" << myName << "\"): INTERCEPTED" << endl;
  }


  size_t MemPoolDebug::ourCheckCtorSizeArg(size_t sz)
  {
    if (sz == 0) CoCoA_THROW_ERROR(ERR::MemPoolZero, "MemPoolDebug ctor");
    return sz;
  }

  //-------------------- constructor & destructor --------------------//

  MemPoolDebug::MemPoolDebug(size_t sz, const string& name, size_t DebugMargin):
      myAliveOrDead(AliveMark),
      myName(name),
      mySliceSizeReq(ourCheckCtorSizeArg(sz)),
      myMarginWords(DebugMargin),
      mySliceWords(2*myMarginWords+1+(sz-1)/WordSize),
      mySliceBytes(mySliceWords*WordSize),
      myMemMgr(mySliceBytes, name, MemPoolFast::FillNewLoaf), // last arg forces filling of newly allocated loaves
      myHeadOfUsedList(nullptr)
  {
    myDebugLevel = ourInitialDebugLevel;
    myVerbosityLevel = ourInitialVerbosityLevel; // default initial verbosity
    myMemMgr.SetVerbosityLevel(myVerbosityLevel);

    myAllocCount = 0;
    myAllocWatchPoint = 0;
    myFreeCount = 0;
    myFreeWatchPoint = 0;
    myInUseCount = 0;
    myInUseMax = 0;
    myNextOutputStatusTime = ourOutputStatusInterval*ceil((0.1+CpuTime())/ourOutputStatusInterval); // make sure value is strictly positive.
    if (myVerbosityLevel > 0)
      LogStream() << "[Log] MemPoolDebug(\"" << myName << "\"):   ctor completed,"
                  << "  RequestedSize=" << mySliceSizeReq
                  << "  ActualSize=" << mySliceBytes << endl;
  }


  inline void MemPoolDebug::myAliveCheck()
  {
    // This if block is an ugly hack to produce the regular reports Abshoff wants.
    if (myAliveOrDead == AliveMark && CpuTime() > myNextOutputStatusTime)
    {
      myMemMgr.myOutputStatus();
      double t = CpuTime();
      while (myNextOutputStatusTime <= t)
        myNextOutputStatusTime += ourOutputStatusInterval;
    }
    // End of ugly hack.
    if (myAliveOrDead == AliveMark) return;
    CoCoA_THROW_ERROR(ERR::DeadMemPool, "MemPoolDebug::myAliveCheck");
  }


  MemPoolDebug::~MemPoolDebug()
  {
    myAliveCheck();
    if (myVerbosityLevel > 0) LogStream() << "[Log] MemPoolDebug(\"" << myName << "\"):   dtor commencing" << endl;
    if (myInUseCount != 0) DtorErrorMesg();
    if (myVerbosityLevel > 0) PrintStats();
    if (myVerbosityLevel > 0) LogStream() << "[Log] MemPoolDebug(\"" << myName << "\"):   dtor completed" << endl;
    myAliveOrDead = ~AliveMark;  // Useful telltale for debugging
  }



  //-------------------- alloc & free --------------------//

  void* MemPoolDebug::alloc()
  {
    return alloc(mySliceSizeReq);
  }


  void* MemPoolDebug::alloc(size_t sz)
  {
    myAliveCheck();
    if (myDebugLevel > 1) FullOverwriteFreeCheck();
    if (++myAllocCount == myAllocWatchPoint) intercepted();
    if (sz != mySliceSizeReq)
    {
      void* ptr = ::operator new(sz);
      AllocWrongSizeMesg(sz, ptr); // mesg always printed in debug mode
      return ptr;
    }
    if (mySliceSizeReq > MaxSliceSize)
    {
      // No margins if handled by ::op new; let valgrind deal with any problems ;-)
      void* ptr = ::operator new(mySliceSizeReq);
      if (myVerbosityLevel > 2) AllocLargeSliceMesg(ptr);
      return ptr;
/////      slice_t ptr = myMarginWords+static_cast<slice_t>(::operator new(mySliceBytes));
/////      AllocLargeSliceMesg(ptr);
/////      AllocMark(ptr);
/////      return ptr;
    }

    slice_t ptr = static_cast<slice_t>(myMemMgr.alloc());
    if (myDebugLevel > 0)
    {
      OverwriteFreeCheck(ptr);
      AllocMark(ptr);
    }

    ++myInUseCount;
    if (myInUseCount > myInUseMax) myInUseMax = myInUseCount;
    if (myVerbosityLevel > 2) AllocMesg(ptr+myMarginWords);
    return ptr + myMarginWords;
  }


  void MemPoolDebug::free(void* ptr)
  {
    free(ptr, mySliceSizeReq);
  }


  void MemPoolDebug::free(void* ptr, size_t sz)
  {
    myAliveCheck();
    if (myDebugLevel > 1) FullOverwriteFreeCheck();
    if (++myFreeCount == myFreeWatchPoint) intercepted();
    if (ptr == nullptr)
    {
      if (myDebugLevel > 0 || myVerbosityLevel > 1) FreeZeroPtrMesg();
      return;
    }
    if (sz != mySliceSizeReq)
    {
      FreeWrongSizeMesg(sz, ptr); // mesg always printed in debug mode
      ::operator delete(ptr);
      return;
    }
    if (mySliceSizeReq > MaxSliceSize)
    {
      // No margins were allocated; we let valgrind check memory handled by the system
      if (myVerbosityLevel > 2) FreeLargeSliceMesg(ptr);
      ::operator delete(ptr /*, mySliceSizeReq*/);
      return;
/////      slice_t true_ptr = static_cast<slice_t>(ptr) - myMarginWords;
/////      FreeLargeSlice(true_ptr);
/////      FreeMark(true_ptr);
/////      ::operator delete(true_ptr /*, mySliceBytes*/);
/////      return;
    }

    if (myDebugLevel > 0 && FreeError(ptr)) return;
    --myInUseCount;
    if (myVerbosityLevel > 2) FreeMesg(ptr);

    // p points to the start of the first margin.
    slice_t p = static_cast<slice_t>(ptr) - myMarginWords;
    if (myDebugLevel > 0) FreeMark(p);

    if (myDebugLevel > 1)
    {
      // high debug level -- place freed block on special list of freed blocks
      slice_t old_head = myHeadOfUsedList;
      myHeadOfUsedList = p;
      *myHeadOfUsedList = old_head;
    }
    else
    {
      // low debug level -- return block to underlying MemPoolFast
      myMemMgr.free(p);
    }
  }


  void MemPoolDebug::InterceptAlloc(size_t nth)
  {
    myAliveCheck();
    myAllocWatchPoint = nth;
  }


  void MemPoolDebug::InterceptFree(size_t nth)
  {
    myAliveCheck();
    myFreeWatchPoint = nth;
  }


  void MemPoolDebug::SetDebugLevel(unsigned int lev)
  {
    myAliveCheck();
    if (myAllocCount == 0) // can change debug level only before allocating any space
      myDebugLevel = lev;
  }


  void MemPoolDebug::SetVerbosityLevel(unsigned int lev)
  {
    myAliveCheck();
    myVerbosityLevel = lev;
    myMemMgr.SetVerbosityLevel(lev);
  }


  // ------------------------------------------------------------------
  // Fake MemPool -- temporary hack to allow threadsafety.

  MemPoolFake::MemPoolFake(std::size_t sz, const std::string& name):
      mySliceSizeReq(sz),
      myName(name)
  {}
  MemPoolFake::~MemPoolFake() {}
  void* MemPoolFake::alloc() { return ::operator new(mySliceSizeReq); }
  void* MemPoolFake::alloc(std::size_t sz) { return ::operator new(sz); }
  void MemPoolFake::free(void* ptr) { ::operator delete(ptr); }
  void MemPoolFake::free(void* ptr, std::size_t) { ::operator delete(ptr); }
//    void MemPoolFake::myOutputStatus() const { LogStream() << "FAKE MemPool" << endl; }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MemPool.C,v 1.18 2022/02/18 14:11:55 abbott Exp $
// $Log: MemPool.C,v $
// Revision 1.18  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.17  2022/02/13 09:58:06  abbott
// Summary: Use uintptr_t instead of unsigned long (after email from Bruns 2022-02-12)
//
// Revision 1.16  2021/10/30 16:51:46  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.15  2021/02/17 17:53:22  abbott
// Summary: Increased max slice size (see redmine 1522)
//
// Revision 1.14  2021/01/07 15:16:51  abbott
// Summary: Corrected copyright
//
// Revision 1.13  2020/12/04 10:13:21  abbott
// Summary: Updated: big slices now handled by sys mem mgr (redmine 1522)
//
// Revision 1.12  2020/06/17 15:49:24  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.11  2019/10/15 11:54:08  abbott
// Summary: Changed 0 into nullptr (where appropriate)
//
// Revision 1.10  2019/03/04 10:19:05  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.9  2017/12/01 17:30:30  bigatti
// -- reduced loaf sizes (was tested a long time ago...)
//
// Revision 1.8  2015/11/04 10:43:55  abbott
// Summary: Added a space in a warning message (about ZERO PTR)
//
// Revision 1.7  2012/04/23 12:18:26  abbott
// Now use std library fill fn inside loaf ctor.
// Code is shorter & clearer (& maybe faster).
//
// Revision 1.6  2012/04/23 10:29:25  abbott
// Added MemPoolFake for use when CoCoA_THREADSAFE_HACK is set.
//
// Revision 1.5  2010/12/20 11:47:15  abbott
// Updated MemPool so user can specify on which streams to print
// logging info and errors (when in verbose mode).
//
// Revision 1.4  2010/12/17 16:13:24  abbott
// Minor cosmetic change to the format of logging and error messages.
//
// Revision 1.3  2010/09/15 21:30:29  abbott
// Added extra print stmt when printing status of an empty MemPool.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.4  2007/03/07 18:08:02  cocoa
// Removed magic constant; added some comments.
//
// Revision 1.3  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/10/06 09:48:52  cocoa
// Changed data member order in a loaf; ctor now initializes all data members
// before ctor body.  Fixed problem with value of myNextOutputStatusTime: now
// it is set to value strictly in the future.  Changed some longs into unsigned
// longs.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.8  2006/05/15 16:28:49  cocoa
// Fixed bug which appeared only when there was a double free.
//
// Revision 1.7  2006/04/27 13:54:19  cocoa
// Improved MemPools with an experimental way of handling raw memory in an
// exception clean manner.
//
// Revision 1.6  2006/04/11 16:34:29  cocoa
// Completed overhaul of MemPool (incl documentation).
// Modified GMPAllocator so that you can specify the slice size
// in the ctor -- useful for experimentation.
//
// Revision 1.5  2006/04/07 16:44:52  cocoa
// Considerably updated MemPool design -- it works, and I'm about to test
// its efficiency against the old one.
//
// Revision 1.4  2006/03/31 16:03:03  cocoa
// Further refinement to MemPool.  Probably more to come (on Monday).
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
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.5  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from ZZ).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.4  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/03/29 17:36:47  cocoa
// Just checking in before going home -- test-matrix1 does not yet link!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.6  2004/11/19 15:14:09  cocoa
// (a) Added new check to MemPool so that it can signal an
//     error if one tries to use a MemPool after it has been
//     destroyed.
// (b) Improved makefile in TEST/ so that it checks output,
//     and prints useful messages if the test fails.
// (c) Tidied ring.txt a bit (still more to do).
//
// Revision 1.5  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.4  2004/11/05 15:32:47  cocoa
// Minor finishing touches.
//   Better spacing in log messages printed out by MemPools.
//   Added an idea for a future improvement.
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
// Revision 1.9  2003/06/23 16:19:59  abbott
// Minor cleaning prior to public release.
//
// Revision 1.8  2003/05/08 08:59:49  abbott
// Now handles the case of a MemPool on blocks of size 0 without crashing.
// Strictly I believe this should provoke an error, but currently the size
// 0 is silently replaced by 1.  In debugging mode a "rude message" is
// printed out when a MemPool for block of size 0 is constructed.
//
// The debugging messages have slightly improved formatting now.
//
// Revision 1.7  2002/10/25 13:06:53  abbott
// Added implementations of the new alloc/free member functions which
// require no size parameter (two implementations of each, for the normal
// and the debugging versions).
//
// Cleaned up the use of "using": previously it was a heavy-handed
// "using namespace std;" now only those names needed are listed.
//
// Removed the definition of MemPool::MEMPOOL_ALLOCATED_WORD_INSIDE since
// now a value depending on the slice is used: the function AllocMark has
// been modified accordingly -- the public area of each slice is filled
// with the ones-complement of the pointer to the "before" margin.
//
// The PrintSlice procedure now prints on "cerr" rather than "cout".
//
// Revision 1.6  2002/03/05 17:35:45  abbott
// Added two new debugging functions: InterceptAlloc and InterceptFree
// which allow one easily to intercept the nth call to alloc/free.
//
// Added a new field containing the number of loaves created by each
// MemPool object: this is used to determine the growth of the size
// of each loaf (doubles every second time).
//
// Modified slightly some of the printed messages: better uniformity,
// and slightly shorter lines of output.
//
// OverwriteFreeCheck now resets any corrupted blocks back to their
// expected state (i.e. full of MEMPOOL_FREE_WORD values).
//
// AlreadyFreed now checks only the margins to decide whether the
// block has already been freed.  If it looks like an already-freed
// block then it performs an OverwriteFreeCheck on the block before
// reporting it as a multiply-freed block (so you may get two messages
// for the same block).
//
// Revision 1.5  2002/01/11 11:57:23  abbott
// Modified the progress messages printed by AllocMesg and FreeMesg
// so that it should be possible to use the new_delete program to
// find which allocations were never freed.  Also modified the
// messages warning about "wrong size" allocations; these now print
// out the pointer concerned.  Added a check for freeing the zero
// pointer which is probably an error for a MemPool: prints out a
// warning message if the debug level is at least 1.
//
// Revision 1.4  2001/12/03 13:36:08  abbott
// Slight changes to allow "level 0" debugging:
//  * almost as fast as without debugging;
//  * computes stats;
//  * checks that the MemPool is empty upon destruction.
//
// Revision 1.3  2001/11/28 21:31:53  abbott
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
// Revision 1.2  2001/11/16 19:26:08  bigatti
// added:  using namespace std;
// for compatibility with gcc-3
//
// Revision 1.1  2001/10/04 14:55:00  abbott
// Initial revision
//
