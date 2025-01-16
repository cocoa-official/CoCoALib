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

#endif
