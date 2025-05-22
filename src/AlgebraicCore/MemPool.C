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
      myBegin(static_cast<slice_t*>(std::malloc(NumSlices*SliceBytes))),
      myEnd(myBegin + NumSlices*mySliceWords),
      myFreeCounter(0)
  {
    if (myBegin == nullptr)  throw std::bad_alloc(); //?? CoCoA_THROW_ERROR
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
  { std::free(myBegin); } // myBegin is pointer to raw memory


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
