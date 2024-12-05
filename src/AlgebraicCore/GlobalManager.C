//   Copyright (c)  2007,2010,2011,2016  John Abbott and Anna M. Bigatti

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


#include "CoCoA/GlobalManager.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/DenseUPolyRing.H" // for Hilbert-Poincare' series
#include "CoCoA/utils.H"

#include "CoCoA/PREPROCESSOR_DEFNS.H"
#ifdef CoCoA_WITH_GFAN
#include "gfanlib/gfanlib.h"
#endif

#include "CoCoA/error.H"
#include "TmpHilbertDir/TmpPoincareCPP.H" // for Hilbert-Poincare' series

#include "gmp.h"

#include <algorithm>
using std::min;
using std::max;
#include <cstdlib>
using std::malloc;
using std::realloc;
using std::free;
#include <iostream>
// using std::cerr & std::endl for serious warning in GlobalManager dtor
#include <cstring>
using std::memcpy;


// These 3 fns must have C linkage to work with GMP's mem mgr setter.
extern "C"
{
  void* CoCoA_GMP_alloc(size_t sz);
  void* CoCoA_GMP_realloc(void* ptr, size_t oldsz, size_t newsz);
  void CoCoA_GMP_free(void* ptr, size_t sz);
}

void* CoCoA_GMP_alloc(size_t sz)
{
  CoCoA_ASSERT(CoCoA::GlobalGMPSliceSize() != 0);
  if (sz <= CoCoA::GlobalGMPSliceSize())
    return CoCoA::GlobalGMPPoolPtr()->alloc();
  return malloc(sz);
}

void* CoCoA_GMP_realloc(void* ptr, size_t oldsz, size_t newsz)
{
  CoCoA_ASSERT(CoCoA::GlobalGMPSliceSize() != 0);
  if (oldsz <= CoCoA::GlobalGMPSliceSize() &&
      newsz <= CoCoA::GlobalGMPSliceSize())
    return ptr;

  if (oldsz > CoCoA::GlobalGMPSliceSize() &&
      newsz > CoCoA::GlobalGMPSliceSize())
    return realloc(ptr, newsz);

  const size_t n = min(oldsz, newsz);
  void* dest = CoCoA_GMP_alloc(newsz);
  memcpy(dest, ptr, n);
  CoCoA_GMP_free(ptr, oldsz);
  return dest;
}

void CoCoA_GMP_free(void* ptr, size_t sz)
{
  CoCoA_ASSERT(CoCoA::GlobalGMPSliceSize() != 0);
  if (sz <= CoCoA::GlobalGMPSliceSize())
    CoCoA::GlobalGMPPoolPtr()->free(ptr);
  else
    free(ptr);
}



namespace CoCoA
{

  // Pseudo-ctors for RingZZ and RingQ.
  ring MakeUniqueInstanceOfRingZZ(); // Defined in RingZZ.C.
  FractionField MakeUniqueInstanceOfRingQQ(const ring&); // Defined in RingQQ.C.

  // Checking fns to be called immediately before calling dtors for RingQ and RingZZ
  bool RingZZStillInUse(const ring& ZZ);  // Defined in RingZZ.C
  bool RingQQStillInUse(const FractionField& Q);  // Defined in RingQ.C


  // The static members of GlobalManager -- effectively global variables.
  bool GlobalManager::DtorFailed = false;
  GlobalManager* GlobalManager::ourGlobalDataPtr = nullptr;
  std::size_t GlobalManager::GMPSliceSize = 0; // size in bytes of slices in the MemPool (compile-time constant)
  MemPool* GlobalManager::GMPPoolPtr = nullptr;
  long GlobalManager::ourHPMaxPower = 100;  // for Hilbert-Poincare' series
  bool GlobalManager::ourAllowObsolescentFnsFlag = false;

  GlobalManager* GlobalManager::ptr(const ErrorContext& ErrCtx)
  {
    if (GlobalManager::ourGlobalDataPtr == nullptr)
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::NoGlobalMgr, ErrCtx);
    return GlobalManager::ourGlobalDataPtr;
  }


  GlobalManager::GMPMemMgr::GMPMemMgr(GlobalSettings::GMPAllocatorType choice, std::size_t SliceSize)
  {
    if (choice == GlobalSettings::GMPAllocatorType::SystemDefault) return;

    myPoolPtr.reset(new MemPool(SliceSize, "Global GMP MemPool")); // must do this first to be exception safe
    GlobalManager::GMPPoolPtr = myPoolPtr.get();
    GlobalManager::GMPSliceSize = GlobalManager::GMPPoolPtr->mySliceSize();
    mp_get_memory_functions(&myPrevAlloc, &myPrevRealloc, &myPrevFree);
    mp_set_memory_functions(&CoCoA_GMP_alloc, &CoCoA_GMP_realloc, &CoCoA_GMP_free);
  }


  GlobalManager::GMPMemMgr::~GMPMemMgr()
  {
    if (myPoolPtr.get() == nullptr) return;

    mp_set_memory_functions(myPrevAlloc, myPrevRealloc, myPrevFree);
    GlobalManager::GMPSliceSize = 0;
    GlobalManager::GMPPoolPtr = nullptr;
  }


  // ----------------------------------------------------------------------

  const std::size_t GlobalSettings::ourDefaultSliceSize = 2*sizeof(long);
  const GlobalSettings::ResidueRepr GlobalSettings::ourDefaultResidueRepr = GlobalSettings::ResidueRepr::symmetric;
  const GlobalSettings::GMPAllocatorType GlobalSettings::ourDefaultGMPAllocatorType = GlobalSettings::GMPAllocatorType::SystemDefault;


  GlobalSettings::GlobalSettings():
      myResidueReprHasBeenSet(false),
      myGMPAllocatorTypeHasBeenSet(false),
      mySliceSizeHasBeenSet(false),
      myObsolescentFnPolicyHasBeenSet(false),
      myResidueRepr(ourDefaultResidueRepr),
      myGMPAllocatorType(ourDefaultGMPAllocatorType),
      mySliceSize(ourDefaultSliceSize),
      myObsolescentFnPolicy(ObsolescentFnPolicy::forbid)
  {}

  GlobalSettings& GlobalSettings::mySetResidueRepr(ResidueRepr r)
  {
    CoCoA_ASSERT(!myResidueReprHasBeenSet);
    myResidueReprHasBeenSet = true;
    myResidueRepr = r;
    return *this;
  }

  GlobalSettings& GlobalSettings::mySetGMPAllocatorType(GMPAllocatorType a)
  {
    CoCoA_ASSERT(!myGMPAllocatorTypeHasBeenSet);
    myGMPAllocatorTypeHasBeenSet = true;
    myGMPAllocatorType = a;
    return *this;
  }

  GlobalSettings& GlobalSettings::mySetSliceSize(std::size_t SliceSize)
  {
    CoCoA_ASSERT(!mySliceSizeHasBeenSet);
    mySliceSizeHasBeenSet = true;
    mySliceSize = SliceSize;
    return *this;
  }

  GlobalSettings& GlobalSettings::mySetObsolescentFnsPolicy(ObsolescentFnPolicy ObsFn)
  {
    CoCoA_ASSERT(!myObsolescentFnPolicyHasBeenSet);
    myObsolescentFnPolicyHasBeenSet = true;
    myObsolescentFnPolicy = ObsFn;
    return *this;
  }

  GlobalSettings GlobalSettings::operator()(std::size_t SliceSize) const
  {
    CoCoA_ASSERT(!mySliceSizeHasBeenSet && myGMPAllocatorType != GMPAllocatorType::SystemDefault);
    GlobalSettings ans(*this);
    return ans.mySetSliceSize(SliceSize);
  }


  GlobalSettings operator+(const GlobalSettings& arg1, const GlobalSettings& arg2)
  {
    GlobalSettings ans;
    if (arg1.myResidueReprHasBeenSet && arg2.myResidueReprHasBeenSet)
      CoCoA_THROW_ERROR(ERR::BadGlobalSettings, "residue setting");
    if (arg1.myResidueReprHasBeenSet) ans.mySetResidueRepr(arg1.myResidueRepr);
    if (arg2.myResidueReprHasBeenSet) ans.mySetResidueRepr(arg2.myResidueRepr);

    if (arg1.myGMPAllocatorTypeHasBeenSet && arg2.myGMPAllocatorTypeHasBeenSet)
      CoCoA_THROW_ERROR(ERR::BadGlobalSettings, "GMP allocator type");
    if (arg1.myGMPAllocatorTypeHasBeenSet) ans.mySetGMPAllocatorType(arg1.myGMPAllocatorType);
    if (arg2.myGMPAllocatorTypeHasBeenSet) ans.mySetGMPAllocatorType(arg2.myGMPAllocatorType);

    if (arg1.mySliceSizeHasBeenSet && arg2.mySliceSizeHasBeenSet)
      CoCoA_THROW_ERROR(ERR::BadGlobalSettings, "GMPAllocator slice size");
    if (arg1.mySliceSizeHasBeenSet) ans.mySetSliceSize(arg1.mySliceSize);
    if (arg2.mySliceSizeHasBeenSet) ans.mySetSliceSize(arg2.mySliceSize);

    if (arg1.myObsolescentFnPolicyHasBeenSet && arg2.myObsolescentFnPolicyHasBeenSet)
      CoCoA_THROW_ERROR(ERR::BadGlobalSettings, "obsolescent fns policy");
    if (arg1.myObsolescentFnPolicyHasBeenSet) ans.mySetObsolescentFnsPolicy(arg1.myObsolescentFnPolicy);
    if (arg2.myObsolescentFnPolicyHasBeenSet) ans.mySetObsolescentFnsPolicy(arg2.myObsolescentFnPolicy);

    return ans;
  }


  const GlobalSettings UseSymmResidues(GlobalSettings().mySetResidueRepr(GlobalSettings::ResidueRepr::symmetric));
  const GlobalSettings UseNonNegResidues(GlobalSettings().mySetResidueRepr(GlobalSettings::ResidueRepr::NonNegative));
  const GlobalSettings UseSystemAllocatorForGMP(GlobalSettings().mySetGMPAllocatorType(GlobalSettings::GMPAllocatorType::SystemDefault));
  const GlobalSettings UseGMPAllocator(GlobalSettings().mySetGMPAllocatorType(GlobalSettings::GMPAllocatorType::cocoa));
  const GlobalSettings ForbidObsolescentFns(GlobalSettings().mySetObsolescentFnsPolicy(GlobalSettings::ObsolescentFnPolicy::forbid));
  const GlobalSettings AllowObsolescentFns(GlobalSettings().mySetObsolescentFnsPolicy(GlobalSettings::ObsolescentFnPolicy::allow));


  // ----------------------------------------------------------------------


  GlobalManager::ZZQQMgr::ZZQQMgr():
      myRingZZ(MakeUniqueInstanceOfRingZZ()),
      myRingQQ(MakeUniqueInstanceOfRingQQ(myRingZZ))
  {}

  GlobalManager::ZZQQMgr::~ZZQQMgr()
  {
    if (RingZZStillInUse(myRingZZ) || RingQQStillInUse(myRingQQ))
      DtorError(); // *IMPORTANT* cannot throw here -- inside a dtor!
  }


  // ----------------------------------------------------------------------

  GlobalManager::GlobalManager(const GlobalSettings& settings):
      myResidueRepr(settings.myResidueRepr),
      myGMPMemMgr(settings.myGMPAllocatorType, settings.mySliceSize),
      myZZQQMgr()
  {
// !!!***NOT THREAD SAFE***!!!  Must make next 3 lines atomic.
    // Complain if a GlobalManager object has already been created
    if (ourGlobalDataPtr != nullptr)
      CoCoA_THROW_ERROR(ERR::GlobalManager2, "GlobalManager ctor");
    ourAllowObsolescentFnsFlag = (settings.myObsolescentFnPolicy == GlobalSettings::ObsolescentFnPolicy::allow);
#ifdef CoCoA_WITH_GFAN
    gfan::initializeCddlibIfRequired();
#endif
    ourGlobalDataPtr = this;  // this line MUST BE LAST!
  }


  GlobalManager::~GlobalManager()
  {
    // Delete registered globals in reverse order
    while (!myDtorStack.empty())
    {
      myDtorStack.top().RunDtor(); // try...catch????
      myDtorStack.pop();
    }
#ifdef CoCoA_WITH_GFAN
    gfan::deinitializeCddlibIfRequired();
#endif
    ourGlobalDataPtr = nullptr; // "deregister" the global data
  }


  void GlobalManager::DtorError()
  {
    DtorFailed = true;
    std::cerr << std::endl
              << "============================================" << std::endl
              << ">>> CoCoA: PROBLEM DURING FINAL CLEAN-UP <<<" << std::endl
              << "============================================" << std::endl
              << std::endl
              << "--------------------------------------------------------------" << std::endl
              << ">>>  CoCoA::GlobalManager dtor: CoCoA objects still live!  <<<" << std::endl
              << "--------------------------------------------------------------" << std::endl
              << std::endl;
  }


  GlobalSettings::ResidueRepr DefaultResidueRepr()
  {
    return GlobalManager::ptr(CoCoA_ERROR_CONTEXT)->myResidueRepr;
  }

  
  //----------------------------------------------------------------------
  // pre-computed power list for univariate Hilbert-Poincare Series
  void MakeGlobalHPPowerList(const DenseUPolyRing& P)
  {  
    MakeHPPowerList(GlobalManager::ptr(CoCoA_ERROR_CONTEXT)->myHPPowerList,
                    P,
                    GlobalManager::ourHPMaxPower);
  }


  long HPPowerListMaxDeg()
  {  
    return GlobalManager::ptr(CoCoA_ERROR_CONTEXT)->ourHPMaxPower;
  }

  
  ConstRefRingElem HPPowerList(int exp)
  {
    if (exp>HPPowerListMaxDeg())
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "HPPowerList");
    return GlobalManager::ptr(CoCoA_ERROR_CONTEXT)->myHPPowerList[exp];
  }
  
    
  void CopyHPPower(RingElem& res, int exp)
  {
    if (exp<=HPPowerListMaxDeg())
      res = HPPowerList(exp);
    else
    {
      res = HPPowerList(HPPowerListMaxDeg());
      const DenseUPolyRing HSRing = owner(res);
      for (long i=HPPowerListMaxDeg(); i<exp; ++i)
        HSRing->myMulBy1MinusXExp(raw(res), 1);
    }    
  }
  
  RandomSource& GlobalRandomSource()
  {
    return GlobalManager::ptr(CoCoA_ERROR_CONTEXT)->myRandomSource;
  }

  //------------------------------------------------------------------
  // Things related to registration of pseudo-dtors for globals.

  GlobalManager::PseudoDtor::PseudoDtor(void (*dtor)()):
      Dtor0arg(dtor),
      Dtor1arg(nullptr),
      ObjPtr(nullptr)
  {}

  GlobalManager::PseudoDtor::PseudoDtor(void (*dtor)(void*), void* ptr):
      Dtor0arg(nullptr),
      Dtor1arg(dtor),
      ObjPtr(ptr)
  {}


  void GlobalManager::PseudoDtor::RunDtor()
  {
    CoCoA_ASSERT((Dtor0arg == nullptr)^(Dtor1arg == nullptr));
    CoCoA_ASSERT(Dtor1arg == nullptr || ObjPtr != nullptr);
    if (Dtor0arg) Dtor0arg();
    else Dtor1arg(ObjPtr);
    // Clear all pointers -- unnecessary, but may help debugging?
    Dtor0arg = nullptr;
    Dtor1arg = nullptr;
    ObjPtr = nullptr;
  }


  void RegisterDtorForGlobal(void (*dtor)())
  {
    GlobalManager::ptr(CoCoA_ERROR_CONTEXT)->myDtorStack.push(GlobalManager::PseudoDtor(dtor));
  }

  void RegisterDtorForGlobal(void (*dtor)(void*), void* ptr)
  {
    GlobalManager::ptr(CoCoA_ERROR_CONTEXT)->myDtorStack.push(GlobalManager::PseudoDtor(dtor,ptr));
  }

} // end of namespace CoCoA
