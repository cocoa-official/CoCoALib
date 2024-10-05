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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/GlobalManager.C,v 1.44 2024/07/22 15:43:57 abbott Exp $
// $Log: GlobalManager.C,v $
// Revision 1.44  2024/07/22 15:43:57  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.43  2024/03/16 10:15:12  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.42  2022/02/25 10:37:40  abbott
// Summary: Added conditional include (after having removed it from ExtLib-GFan.H)
//
// Revision 1.41  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.40  2021/03/04 21:03:45  abbott
// Summary: enum revision and renaming (redmine 894)
//
// Revision 1.39  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.38  2021/01/07 15:07:03  abbott
// Summary: Corrected copyright
//
// Revision 1.37  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.36  2020/01/24 21:38:13  abbott
// Summary: Added init and de-init for GFan (if necessary)
//
// Revision 1.35  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.34  2017/07/08 19:05:51  abbott
// Summary: major revision to interrupt mechanism
//
// Revision 1.33  2017/06/22 11:02:04  abbott
// Summary: Minor change to format of "imminent disaster" message
//
// Revision 1.32  2017/03/13 12:20:32  abbott
// Summary: Removed CPPFlags_check; function subsumed by PREPROCESSOR_DEFNS.H
//
// Revision 1.31  2017/03/01 17:16:23  abbott
// Summary: Added automatic check for some CPP flags (THREADSAFE_HACK)
//
// Revision 1.30  2016/11/18 18:10:35  abbott
// Summary: Renamed InterruptFlag to InterruptSignalReceived
//
// Revision 1.29  2016/11/04 20:41:14  abbott
// Summary: Added stuff to allow user to enable/disable calling obsolescent fns
//
// Revision 1.28  2016/11/03 12:29:58  abbott
// Summary: Added file for obsolescent fns; also there is a global flag saying whether to give error if calling one.
//
// Revision 1.27  2016/09/21 14:24:39  abbott
// Summary: Added GlobalManagerDtorFailed (and global var GlobalManager::DtorFailed)
//
// Revision 1.26  2016/03/21 17:03:18  abbott
// Summary: Changed "imminent disaster" into a more helpful message
//
// Revision 1.25  2015/11/30 21:56:43  abbott
// Summary: Moved "imminent disaster" mesg into separate fn DtorError
//
// Revision 1.24  2015/11/04 10:08:05  abbott
// Summary: Added RegisterDtorForGlobal
//
// Revision 1.23  2015/09/02 11:40:27  abbott
// Summary: Added useful comment
//
// Revision 1.22  2015/06/29 10:24:29  abbott
// Summary: Added GlobalManager::ourInterruptFlag; cleaner impl
// Author: JAA
//
// Revision 1.21  2015/05/21 12:22:41  abbott
// Summary: Added initializer to NULL ptr in ctor for GlobalManager
// Author: JAA
//
// Revision 1.20  2015/05/20 15:37:00  abbott
// Summary: Removed the interrupt flag specifier from GlobalSettings
// Author: JAA
//
// Revision 1.19  2015/05/20 14:49:12  abbott
// Summary: Added fns for specifying the interrupt flag to monitor
// Author: JAA
//
// Revision 1.18  2014/07/09 13:01:17  abbott
// Summary: Removed AsDenseUPolyRing
// Author: JAA
//
// Revision 1.17  2014/07/01 12:40:53  bigatti
// -- added CopyHPPower and argument check in HPPowerList
//
// Revision 1.16  2014/04/30 16:07:40  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.15  2013/06/17 08:54:02  abbott
// Added RegisterDtorForGlobal.
//
// Revision 1.14  2012/10/15 12:35:24  abbott
// Added  std::  prefix.
//
// Revision 1.13  2012/02/08 13:47:16  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.12  2011/05/19 13:54:48  abbott
// Replaced DefaultResiduesAreSymm by DefaultResidueSetting.
//
// Revision 1.11  2011/05/03 10:03:32  abbott
// Added GlobalRandomSource.
// Internally added GlobalManager::ptr to allow neater implementations.
//
// Revision 1.10  2010/11/17 15:52:33  abbott
// Removed out-of-date include of GMPAllocator.H.
//
// Revision 1.9  2010/11/11 17:45:08  abbott
// Moved GMPMemMgr so that it is a nested class inside GlobalManager.
//
// Revision 1.8  2010/10/29 12:06:41  bigatti
// -- added globals for Hilbert-Poincare' series
//
// Revision 1.7  2010/10/27 20:58:45  abbott
// Major reorganization of GlobalManager and GMPAllocator.
//
// Revision 1.6  2010/10/22 14:03:04  abbott
// Major change to GMPAllocator -- it is now set/activated by the GlobalManager.
// This is a Friday afternoon check-in... hope to check in cleaner code in the
// next few days.
//
// Revision 1.5  2010/09/30 14:28:23  abbott
// Replaced auto_ptrs to RingZ and RingQ by direct values; ctor changed accordingly.
//
// Dtor now checks that ref counts in RingZ and RingQ are correct; if not, a rude
// message is printed on cerr (and the program will probably crash after the
// GlobalManager has been destroyed).
//
// Revision 1.4  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.3  2009/05/14 09:39:29  abbott
// Added possibility to specify "symmetric" or "non-negative" residues
// in quotients of ZZ.  Affects printing of elements in quotients of ZZ
// (also changed printing of elements in general quotient rings).
// Consequent changes in several tests.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2007/03/05 21:33:13  cocoa
// Improved/cleaned GlobalManager; added doc too.
//
// Revision 1.1  2007/03/03 14:02:11  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.1  2007/03/02 16:46:28  cocoa
// New foundations object which calls ctors and dtors of global objects.
//
