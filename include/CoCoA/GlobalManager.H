#ifndef CoCoA_GlobalManager_H
#define CoCoA_GlobalManager_H

//   Copyright (c)  2007,2009-2011,2016,2017  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/FractionField.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/random.H"

#include <memory>
// using std::unique_ptr;
#include <stack>
// using std::stack;

namespace CoCoA
{

  class DenseUPolyRing; // fwd decl -- defined in DenseUPolyRing.H

  // This class simply manages the parameters for the GlobalManager ctor.
  class GlobalSettings
  {
  public:
//    enum class ResidueSetting {SymmResidues, NonNegResidues};
    enum class ResidueRepr {symmetric, NonNegative};
//    static constexpr ResidueRepr SymmResidues = ResidueRepr::symmetric;
//    static constexpr ResidueRepr NonNegResidues = ResidueRepr::NonNegative;
//    enum class AllocatorSetting {GMPAllocator, SystemAllocator};
    enum class GMPAllocatorType {cocoa, SystemDefault};
//    enum class ObsolescentFnSetting {ObsolescentFnsForbidden, ObsolescentFnsAllowed};
    enum class ObsolescentFnPolicy {forbid, allow};

  public:
    GlobalSettings();
    GlobalSettings& mySetResidueRepr(ResidueRepr r);
    GlobalSettings& mySetGMPAllocatorType(GMPAllocatorType a);
    GlobalSettings& mySetSliceSize(std::size_t SliceSize);
    GlobalSettings& mySetObsolescentFnsPolicy(ObsolescentFnPolicy a);
    GlobalSettings operator()(std::size_t SliceSize) const; // NB creates a new object!!

  private: // data members
    bool myResidueReprHasBeenSet;
    bool myGMPAllocatorTypeHasBeenSet;
    bool mySliceSizeHasBeenSet;
    bool myObsolescentFnPolicyHasBeenSet;
    ResidueRepr myResidueRepr;
    GMPAllocatorType myGMPAllocatorType;
    std::size_t mySliceSize;
    ObsolescentFnPolicy myObsolescentFnPolicy;
    static const std::size_t ourDefaultSliceSize;
    static const ResidueRepr ourDefaultResidueRepr;
    static const GMPAllocatorType ourDefaultGMPAllocatorType;

    friend GlobalSettings operator+(const GlobalSettings& arg1, const GlobalSettings& arg2);

    friend class GlobalManager;
  };

  GlobalSettings operator+(const GlobalSettings& arg1, const GlobalSettings& arg2);

  extern const GlobalSettings UseSymmResidues;
  extern const GlobalSettings UseNonNegResidues;
  extern const GlobalSettings UseSystemAllocatorForGMP;
  extern const GlobalSettings UseGMPAllocator;
  extern const GlobalSettings ForbidObsolescentFns;
  extern const GlobalSettings AllowObsolescentFns;


  // At most one instance of this class may exist at any given time.
  class GlobalManager // : private BOOST::noncopyable ???
  {
  public:
    GlobalManager(const GlobalSettings& settings = GlobalSettings());
    ~GlobalManager();
    const ring& myRingZZ() const { return myZZQQMgr.myRingZZ; }
    const FractionField& myRingQQ() const { return myZZQQMgr.myRingQQ; }
    static GlobalManager* ptr(const ErrorContext& ErrCtx);

  public: // disable copy ctor and assignment
    GlobalManager(const GlobalManager&) = delete;
    GlobalManager& operator=(const GlobalManager&) = delete;
  private:
    static void DtorError();
    static bool DtorFailed; // set iff DtorError was called
    friend bool GlobalManagerDtorFailed();
  private: // the true global variables
    static std::size_t GMPSliceSize;          // size in bytes of slices in the GMP MemPool
    static MemPool* GMPPoolPtr;               // raw ptr to GMP MemPool (or nullptr if there is none)
    static GlobalManager* ourGlobalDataPtr; ///< nullptr or addr of unique existing GlobalManager
    friend bool IsInitialized();
    friend bool IsAllowedObsolescentFnCall();

    friend class GMPMemMgr;  // so it can set GMPPoolPtr and GMPSliceSize.
    friend std::size_t GlobalGMPSliceSize();  // accessor fn
    friend MemPool* GlobalGMPPoolPtr();       // accessor fn
  private:
    friend const ring& RingZZ();              ///< defined in RingZZ.C
    friend const FractionField& RingQQ();     ///< defined in RingQQ.C
    friend GlobalSettings::ResidueRepr DefaultResidueRepr();     ///< called by ctors for rings of the form Z/(nZ)
    friend RandomSource& GlobalRandomSource();    

    //-- for Hilbert-Poincare' series --
    friend void MakeGlobalHPPowerList(const DenseUPolyRing& HSRing);
    friend long HPPowerListMaxDeg();
    friend ConstRefRingElem HPPowerList(int exp);
    void CopyHPPower(RingElem& res, int exp);
  private:
    class GMPMemMgr
    {
    public:
      GMPMemMgr(GlobalSettings::GMPAllocatorType choice, std::size_t SliceSize);
      ~GMPMemMgr();
    private: // data members of GMPMemMgr
      std::unique_ptr<MemPool> myPoolPtr;
      void *(*myPrevAlloc) (std::size_t);
      void *(*myPrevRealloc) (void *, std::size_t, std::size_t);
      void (*myPrevFree) (void *, std::size_t);
    };
  private:
    class ZZQQMgr
    {
    public:
      ZZQQMgr();
      ~ZZQQMgr(); // prints rude message if ZZ or QQ still have external references
    public: // data members of ZZQQMgr
      ring myRingZZ;                 // Must come before RingQQ
      FractionField myRingQQ;
    };
  private: // data members of GlobalManager
    const GlobalSettings::ResidueRepr myResidueRepr; ///< default value used for creating rings Z/(nZ)
    GMPMemMgr myGMPMemMgr;        // Must come before myZZQQMgr
    ZZQQMgr myZZQQMgr;            // Must come IMMEDIATELY after myGMPMemMgr
    RandomSource myRandomSource;  // for global random source (must come after myZZQQMgr)
    static bool ourAllowObsolescentFnsFlag; // flag to allow obsolescent fns to be called
    // Parts related to registration of pseudo-dtors for globals
    class PseudoDtor
    {
    public:
      PseudoDtor(void (*dtor)());
      PseudoDtor(void (*dtor)(void*), void* ptr);  // OBSOLESCENT?
      void RunDtor();
    private: // data members
      void (*Dtor0arg)();
      void (*Dtor1arg)(void*);
      void* ObjPtr;
    };
    std::stack<PseudoDtor> myDtorStack;
    friend void RegisterDtorForGlobal(void (*dtor)());
    friend void RegisterDtorForGlobal(void (*dtor)(void*), void* ptr);  // OBSOLESCENT?

    //-- for Hilbert-Poincare' series --
    static long ourHPMaxPower; // ??? static?
    std::vector<RingElem> myHPPowerList;
  };


  inline bool IsInitialized()
  {
    return (GlobalManager::ourGlobalDataPtr != nullptr);
  }


  inline bool IsAllowedObsolescentFnCall()
  {
    return GlobalManager::ptr(CoCoA_ERROR_CONTEXT)->ourAllowObsolescentFnsFlag;
  }


  GlobalSettings::ResidueRepr DefaultResidueRepr();
  //-- for Hilbert-Poincare' series --
  void RegisterDtorForGlobal(void (*dtor)(void*), void* ptr);


  inline bool GlobalManagerDtorFailed()
  {
    return GlobalManager::DtorFailed;
  }


  void MakeGlobalHPPowerList(const DenseUPolyRing& HSRing);
  long HPPowerListMaxDeg();
  ConstRefRingElem HPPowerList(int exp);
  void CopyHPPower(RingElem& res, int exp);
  inline std::size_t GlobalGMPSliceSize() { return GlobalManager::GMPSliceSize; }
  inline MemPool* GlobalGMPPoolPtr() { return GlobalManager::GMPPoolPtr; }
  RandomSource& GlobalRandomSource();
  void RegisterDtorForGlobal(void (*dtor)());

} // end of namespace CoCoA



// RCS header/log in the next few lines.
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/GlobalManager.H,v 1.37 2024/03/16 10:15:11 abbott Exp $
// $Log: GlobalManager.H,v $
// Revision 1.37  2024/03/16 10:15:11  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.36  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.35  2021/10/30 17:06:58  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.34  2021/03/04 21:03:45  abbott
// Summary: enum revision and renaming (redmine 894)
//
// Revision 1.33  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.32  2021/01/07 14:14:42  abbott
// Summary: Corrected copyright
//
// Revision 1.31  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.30  2019/03/18 17:56:49  abbott
// Summary: Changed 0 into nullptr where appropriate
//
// Revision 1.29  2019/03/04 10:19:28  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.28  2017/07/08 19:05:51  abbott
// Summary: major revision to interrupt mechanism
//
// Revision 1.27  2017/03/13 12:17:58  abbott
// Summary: Removed CPPFLAGS_check stuff; function subsumed by PREPROCESSOR_DEFNS.H
//
// Revision 1.26  2017/03/01 17:16:23  abbott
// Summary: Added automatic check for some CPP flags (THREADSAFE_HACK)
//
// Revision 1.25  2016/11/18 18:07:15  abbott
// Summary: Changed name of InterruptFlag to InterruptSignalReceived (since it is not a flag)
//
// Revision 1.24  2016/11/04 20:38:56  abbott
// Summary: Added stuff to allow user to enable/disable calling of obsolescent fns
//
// Revision 1.23  2016/11/03 12:29:58  abbott
// Summary: Added file for obsolescent fns; also there is a global flag saying whether to give error if calling one.
//
// Revision 1.22  2016/09/21 14:24:39  abbott
// Summary: Added GlobalManagerDtorFailed (and global var GlobalManager::DtorFailed)
//
// Revision 1.21  2015/11/30 21:57:15  abbott
// Summary: New mem gn DtorError to print out "imminent disaster"
//
// Revision 1.20  2015/11/04 10:06:52  abbott
// Summary: Added RegisterDtorForGlobal
//
// Revision 1.19  2015/09/02 11:39:25  abbott
// Summary: Added IsInitialized
//
// Revision 1.18  2015/06/29 10:24:29  abbott
// Summary: Added GlobalManager::ourInterruptFlag; cleaner impl
// Author: JAA
//
// Revision 1.17  2015/05/20 15:36:59  abbott
// Summary: Removed the interrupt flag specifier from GlobalSettings
// Author: JAA
//
// Revision 1.16  2015/05/20 14:49:12  abbott
// Summary: Added fns for specifying the interrupt flag to monitor
// Author: JAA
//
// Revision 1.15  2014/07/01 12:40:14  bigatti
// -- added CopyHPPower
//
// Revision 1.14  2013/06/17 08:54:02  abbott
// Added RegisterDtorForGlobal.
//
// Revision 1.13  2012/12/04 20:04:36  abbott
// Includes new unified random header file.
//
// Revision 1.12  2012/02/08 13:37:35  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.11  2011/05/19 13:54:48  abbott
// Replaced DefaultResiduesAreSymm by DefaultResidueSetting.
//
// Revision 1.10  2011/05/03 10:03:32  abbott
// Added GlobalRandomSource.
// Internally added GlobalManager::ptr to allow neater implementations.
//
// Revision 1.9  2010/11/11 17:45:08  abbott
// Moved GMPMemMgr so that it is a nested class inside GlobalManager.
//
// Revision 1.8  2010/11/05 17:39:14  bigatti
// -- fixed name for HPPowerList function declaration
//
// Revision 1.7  2010/10/29 09:36:04  bigatti
// -- added globals for Hilbert-Poincare' series
//
// Revision 1.6  2010/10/27 20:58:45  abbott
// Major reorganization of GlobalManager and GMPAllocator.
//
// Revision 1.5  2010/10/22 14:03:04  abbott
// Major change to GMPAllocator -- it is now set/activated by the GlobalManager.
// This is a Friday afternoon check-in... hope to check in cleaner code in the
// next few days.
//
// Revision 1.4  2010/09/30 14:28:23  abbott
// Replaced auto_ptrs to RingZ and RingQ by direct values; ctor changed accordingly.
//
// Dtor now checks that ref counts in RingZ and RingQ are correct; if not, a rude
// message is printed on cerr (and the program will probably crash after the
// GlobalManager has been destroyed).
//
// Revision 1.3  2009/05/14 09:39:29  abbott
// Added possibility to specify "symmetric" or "non-negative" residues
// in quotients of ZZ.  Affects printing of elements in quotients of ZZ
// (also changed printing of elements in general quotient rings).
// Consequent changes in several tests.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
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

#endif
