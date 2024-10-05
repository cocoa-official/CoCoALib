//   Copyright (c)  2005,2021  John Abbott and Anna M. Bigatti

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


// Implementation file for the class FreeModule


#include "CoCoA/FreeModule.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/ModuleOrdering.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/assert.H"
#include "CoCoA/degree.H"
#include "CoCoA/module.H"

#include <vector>
using std::vector;
#include <algorithm>
using std::max;
//using std::swap;
#include <functional>
using std::ptr_fun;
#include <iostream>
using std::ostream;
#include<iterator>
using std::back_inserter;
#include <memory>
using std::unique_ptr;
#include <new>
//for placement new

namespace CoCoA
{

  const FreeModuleBase* FreeModulePtr(const module& M)
  { return dynamic_cast<const FreeModuleBase*>(M.myRawPtr()); }

  
  const FreeModuleBase* FreeModulePtr(const module& M, const ErrorContext& ErrCtx)
  {
    const FreeModuleBase* ptr = FreeModulePtr(M);
    if (ptr == nullptr)
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqFreeModule, ErrCtx);
    return ptr;
  }


  namespace // anonymous namespace
  {
    inline const RingElem* import(ModuleBase::ConstRawPtr rawx)
    {
      return static_cast<const RingElem*>(rawx.ptr);
    }

    inline RingElem* import(ModuleBase::RawPtr& rawx)
    {
      return static_cast<RingElem*>(rawx.ptr);
    }


    // This is an exception clean "constructor" for module elements
    void CreateZeroVector(void* ptr, long NumCompts, ring R)  // we know that NumCompts > 0
    {
      RingElem* array = static_cast<RingElem*>(ptr);
      long i;
      try
      {
        for (i=0; i < NumCompts; ++i)
          new(&array[i]) RingElem(R);  // placement new
      }
      catch (...)
      {
        while (i-- > 0)
          array[i].~RingElem();
        throw;
      }
    }

    // This is an exception clean "constructor" for module elements
    void CopyVector(void* ptr, long NumCompts, const void* copyme)  // we know that NumCompts > 0
    {
      RingElem* array = static_cast<RingElem*>(ptr);
      const RingElem* rhs = static_cast<const RingElem*>(copyme);
      long i;
      try
      {
        for (i=0; i < NumCompts; ++i)
          new(&array[i]) RingElem(rhs[i]);  // placement new
      }
      catch (...)
      {
        while (i-- > 0)
          array[i].~RingElem();
      }
    }

    void DeleteVector(void* ptr, long NumCompts)  // we know that NumCompts > 0
    {
      RingElem* array = static_cast<RingElem*>(ptr);
      for (long i=NumCompts; i-- > 0;)
        array[i].~RingElem();
    }

  } // end of anonymous namespace


  class FreeModuleImpl: public FreeModuleBase
  {
    // Two typedefs to save typing.
    typedef FGModuleBase::RawPtr RawPtr;
    typedef FGModuleBase::ConstRawPtr ConstRawPtr;

  protected:
    friend FreeModule NewFreeModule(const ring& R, long NumCompts);
    FreeModuleImpl(const ring& R, long NumCompts);
//???    ~FreeModuleImpl();  // can be handy for debugging
  public:
    long myNumCompts() const override;
    const ring& myRing() const override;
    const FreeModule& myAmbientFreeModule() const override;
    const std::vector<ModuleElem>& myGens() const override;
    const std::vector<ModuleElem>& myMinGens(const CpuTimeLimit& /*CheckForTimeout*/) const override;
    const std::vector<ModuleElem>& myTidyGens(const CpuTimeLimit& /*CheckForTimeout*/) const override;

    const ModuleElem& myZero() const override;
    void myNew(RawPtr& rawv) const override;
    void myNew(RawPtr& rawv, ConstRawPtr rawt) const override;
    void myDelete(RawPtr& rawv) const override;                                            ///< destroys v (incl all resources)
    void mySwap(RawPtr& rawv, RawPtr& raww) const override;                                ///< swap(v, w)
    void myAssign(RawPtr& rawlhs, ConstRawPtr rawv) const override;                        ///< lhs = v
    ConstRefRingElem myCompt(const RawPtr& rawv, long pos) const override;                 ///< v[pos]
    void myNegate(RawPtr& rawlhs, ConstRawPtr rawv) const override;                        ///< lhs = -v
    void myAdd(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const override;         ///< lhs = v+w
    void mySub(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const override;         ///< lhs = v-w

    void myMul(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const override; ///< lhs = r*v
    void myDiv(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const override; ///< lhs = (1/r)*v
    void myOutput(std::ostream& out, ConstRawPtr rawv) const override;                     ///< out << v
    void myOutputSelf(std::ostream& out) const override;                                   ///< out << M
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawv) const override;              ///< OMOut << v
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                            ///< OMOut << M
    bool myIsZero(ConstRawPtr rawv) const override;                                        ///< v == 0
//???    bool IsZeroAddMul(RawPtr& rawlhs, RingElemConstRawPtr rawy, ConstRawPtr rawz) const;  // lhs += y*z, result says whether lhs == 0.
    bool myIsEqual(ConstRawPtr, ConstRawPtr) const override;
  ///  void convert(string&, RawPtr) const;

  protected: // data members
    const ring myR;
    const long myNumComptsValue;            // always > 0
    mutable MemPool myMemMgr;               // This must come before myZeroValue.
    std::unique_ptr<ModuleElem> myZeroValue;  // This must come after myMemMgr.
    std::vector<ModuleElem> myGensArray;
    const FreeModule myM; // so that myAmbientFreeModule can work
  };



  FreeModuleImpl::FreeModuleImpl(const ring& R, long n):
      myR(R),
      myNumComptsValue(n),
      myMemMgr(n*sizeof(RingElem), "FreeModuleImpl.myMemMgr"), // MUST HAVE n != 0
      myZeroValue(new ModuleElem(module(this))),
      myM(this)
  {
    CoCoA_ASSERT("FreeModule must have at least one compt" && n > 0);
    CoCoA_ASSERT("FreeModule dim ludicrous" && n < 65536);
    for (long i=0; i < myNumComptsValue; ++i)
    {
      myGensArray.push_back(myZero());
      import(raw(myGensArray[i]))[i] = 1;
    }
    myRefCountZero();
  }


//    FreeModuleImpl::~FreeModuleImpl()
//    {}


  long FreeModuleImpl::myNumCompts() const
  {
    return myNumComptsValue;
  }


  const ring& FreeModuleImpl::myRing() const
  {
    return myR;
  }


  const FreeModule& FreeModuleImpl::myAmbientFreeModule() const
  {
    return myM;
  }


  const std::vector<ModuleElem>& FreeModuleImpl::myGens() const
  {
    return myGensArray;
  }


  const std::vector<ModuleElem>& FreeModuleImpl::myMinGens(const CpuTimeLimit& /*CheckForTimeout*/) const
  {
    return myGensArray;
  }


  const std::vector<ModuleElem>& FreeModuleImpl::myTidyGens(const CpuTimeLimit& /*CheckForTimeout*/) const
  {
    return myGensArray;
  }


  const ModuleElem& FreeModuleImpl::myZero() const
  {
    return *myZeroValue;
  }


  void FreeModuleImpl::myNew(RawPtr& rawv) const
  {
    rawv.ptr = myMemMgr.alloc();
    CreateZeroVector(rawv.ptr, myNumComptsValue, myR);
  }


  void FreeModuleImpl::myNew(RawPtr& rawv, ConstRawPtr rawcopy) const
  {
    rawv.ptr = myMemMgr.alloc();
    CopyVector(rawv.ptr, myNumComptsValue, rawcopy.ptr);
  }


  void FreeModuleImpl::myDelete(RawPtr& rawv) const
  {
    // Kill elements in reverse order (might scramble the MemPool less?)
    DeleteVector(rawv.ptr, myNumComptsValue);
    myMemMgr.free(rawv.ptr);
  }


  void FreeModuleImpl::mySwap(RawPtr& rawv, RawPtr& raww) const
  {
    std::swap(rawv.ptr, raww.ptr);
  }


  ConstRefRingElem FreeModuleImpl::myCompt(const RawPtr& rawv, long pos) const
  {
    CoCoA_ASSERT(0 <= pos && pos < myNumComptsValue);
    const RingElem* vv = import(rawv);
    return vv[pos];
  }


  void FreeModuleImpl::myAssign(RawPtr& rawlhs, ConstRawPtr rawv) const
  {
    RingElem* lvec = import(rawlhs);
    const RingElem* vec = import(rawv);
    for (long i=0; i < myNumComptsValue; ++i)
      lvec[i] = vec[i];
  }


  void FreeModuleImpl::myNegate(RawPtr& rawlhs, ConstRawPtr rawv) const
  {
    RingElem* lvec= import(rawlhs);
    const RingElem* vec = import(rawv);
    for (long i=0; i < myNumComptsValue; ++i)
      lvec[i] = -vec[i];
  }


  void FreeModuleImpl::myAdd(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const
  {
    RingElem* lvec = import(rawlhs);
    const RingElem* vec1 = import(rawv);
    const RingElem* vec2 = import(raww);
    for (long i=0; i < myNumComptsValue; ++i)
      lvec[i] = vec1[i] + vec2[i];
  }


  void FreeModuleImpl::mySub(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const
  {
    RingElem* lvec = import(rawlhs);
    const RingElem* vec1 = import(rawv);
    const RingElem* vec2 = import(raww);
    for (long i=0; i < myNumComptsValue; ++i)
      lvec[i] = vec1[i] - vec2[i];
  }


  void FreeModuleImpl::myMul(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const
  {
    RingElem* lvec = import(rawlhs);
    const RingElem* vec = import(rawv);
    for (long i=0; i < myNumComptsValue; ++i)
      myR->myMul(raw(lvec[i]), rawx, raw(vec[i]));
  }


  void FreeModuleImpl::myDiv(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const
  {
    RingElem* lvec = import(rawlhs);
    const RingElem* vec = import(rawv);
    for (long i=0; i < myNumComptsValue; ++i)
      myR->myDiv(raw(lvec[i]), rawx, raw(vec[i]));
  }


  void FreeModuleImpl::myOutput(std::ostream& out, ConstRawPtr rawv) const
  {
    if (!out) return;  // short-cut for bad ostreams

    // Guaranteed that myNumComptsValue > 0
    const RingElem* vec = import(rawv);
    out << "[" << vec[0];
    for (long i=1; i < myNumComptsValue; ++i)
    {
      out << ", " << vec[i];
    }
    out << "]";
  }


  void FreeModuleImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "FreeModule(" << myR << ", " << myNumComptsValue << ")";
  }


  void FreeModuleImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawv) const
  {
    // Guaranteed that myNumComptsValue > 0
    const RingElem* vec = import(rawv);
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "FreeModuleElement"); // BUG: what should this OMSymbol be???
    myOutputSelf_OM(OMOut);
    for (long i=0; i < myNumComptsValue; ++i)
    {
      OMOut << vec[i];
    }
    OMOut->mySendApplyEnd();
  }


  void FreeModuleImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "FreeModule"); // BUG: what should this OMSymbol be???
    OMOut << myR << myNumComptsValue;
    OMOut->mySendApplyEnd();
  }


  bool FreeModuleImpl::myIsZero(ConstRawPtr rawv) const
  {
    const RingElem* vec = import(rawv);
    for (long i=0; i < myNumComptsValue; ++i)
      if (!IsZero(vec[i])) return false;
    return true;
  }


///???    bool IsZeroAddMul(RawPtr& rawlhs, RingElemConstRawPtr rawy, ConstRawPtr rawz) const;  // lhs += y*z, result says whether lhs == 0.


  bool FreeModuleImpl::myIsEqual(ConstRawPtr rawv, ConstRawPtr raww) const
  {
    const RingElem* vec1 = import(rawv);
    const RingElem* vec2 = import(raww);
    for (long i=0; i < myNumComptsValue; ++i)
      if (vec1[i] != vec2[i]) return false;
    return true;
  }


  //----------------------------------------------------------------------
  // non-member functions

  // fwd decl -- define later in this file
  FreeModule NewFreeModuleSpPR(const SparsePolyRing& P, const ModuleOrdering& O);

  
  FreeModule NewFreeModule(const ring& R, long NumCompts)
  {
    if (NumCompts < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    if (IsSparsePolyRing(R))
      return NewFreeModuleSpPR(R, NewOrdPosn(ordering(PPM(R)), NumCompts));
    return FreeModule(new FreeModuleImpl(R, NumCompts));
  }


  const std::vector<ModuleElem>& CanonicalBasis(const FreeModule& F)
  { return gens(F); }


  long FirstNonZeroPosn(const ModuleElem& v)
  {
    const long n = NumCompts(owner(v));
    for (long i=0; i < n; ++i )  if ( !IsZero(v[i]) )  return i;
    CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    return 0; // just to keep the compiler quiet
  }


  RingElem FirstNonZero(const ModuleElem& v)
  {
    const long n = NumCompts(owner(v));
    for (long i=0; i < n; ++i )  if ( !IsZero(v[i]) )  return v[i];
    CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    return v[0]; // just to keep the compiler quiet
  }


  class FreeModuleSpPRImpl: public FreeModuleImpl
  {
    // should I use RawPtrs? it would make the code uglier
    //    typedef FGModuleBase::RawPtr RawPtr;
    //    typedef FGModuleBase::ConstRawPtr ConstRawPtr;

    class CmpBase
    {
    public:
      CmpBase(const std::vector<degree>& shifts);
      virtual ~CmpBase() {};
    public:
      virtual long myLPosn(const ModuleElem& v) const =0;
      virtual int myCmpLT(const ModuleElem& v1, const ModuleElem& v2) const =0;
      void myWDeg(degree& d, const ModuleElem& v) const;
      int myCmpWDeg(const ModuleElem& v1, const ModuleElem& v2) const;
    protected:
      std::vector<degree> myShiftsValue;
    };

    // ???ANNA: this is NOT degree compatible
    class PosnOrdImpl: public CmpBase
    {
    public:
      PosnOrdImpl(const std::vector<degree>& shifts);
      virtual ~PosnOrdImpl() {};
      long myLPosn(const ModuleElem& v) const override;
      int myCmpLT(const ModuleElem& v1, const ModuleElem& v2) const override;
      int myCmpPosn(const long pos1, const long pos2) const;
    private:
      std::vector<long> myComptsOrd;
    };

    class OrdPosnImpl: public CmpBase
    {
    public:
      OrdPosnImpl(const std::vector<degree>& shifts);
      virtual ~OrdPosnImpl() {};
      long myLPosn(const ModuleElem& v) const override;
      int myCmpLT(const ModuleElem& v1, const ModuleElem& v2) const override;
      int myCmpPosn(const long pos1, const long pos2) const;
    private:
      std::vector<long> myComptsOrd;
    };


    class WDegPosnOrdImpl: public CmpBase
    {
    public:
      WDegPosnOrdImpl(const std::vector<degree>& shifts);
      virtual ~WDegPosnOrdImpl() {};
      long myLPosn(const ModuleElem& v) const override;
      int myCmpLT(const ModuleElem& v1, const ModuleElem& v2) const override;
      int myCmpPosn(const long pos1, const long pos2) const;
    private:
      std::vector<long> myComptsOrd;
    };

    //--- end of CmpBase ----------------------------------------------------

  private:
    FreeModuleSpPRImpl(const SparsePolyRing& P, const ModuleOrdering& O);
//???    ~FreeModuleSpPRImpl();  // can be handy for debugging
    friend FreeModule NewFreeModuleSpPR(const SparsePolyRing& P, long NumCompts);
    friend FreeModule NewFreeModuleSpPR(const SparsePolyRing& P, const ModuleOrdering& O);
    friend FreeModule NewFreeModule(const ring& P, long NumCompts);
    friend FreeModule NewFreeModule(const ring& P, const std::vector<degree>& shifts);
    friend FreeModule NewFreeModule(const ring& P, const std::vector<degree>& shifts, const ModuleOrderingCtor& OrdCtor);
    friend FreeModule NewFreeModule(const ring& P, long NumCompts, const ModuleOrderingCtor& OrdCtor);

  public:
    const ModuleOrdering& myModuleOrdering() const;
    long myLPosn(const ModuleElem&  v) const;
    int myCmpLT(const ModuleElem&  v1, const ModuleElem&  v2) const;
    void myWDeg(degree& d, const ModuleElem&  v) const;
    int myCmpWDeg(const ModuleElem&  v1, const ModuleElem&  v2) const;

  private: // data members
    const ModuleOrdering myModuleOrd; ///< abstract description of the ordering
    std::unique_ptr<CmpBase> myOrdPtr; ///< actual implementation of the ordering [should be const???]
  };


  class FreeModuleSpPR: public FreeModule  // copied from SparsePolyRing
  {
  public:
    FreeModuleSpPR(const module& M);
    explicit FreeModuleSpPR(const FreeModuleSpPRImpl* MPtr);
    // Default copy ctor works fine.
    // Default dtor works fine.
  public: // disable assignment
    FreeModuleSpPR& operator=(const FreeModuleSpPR& rhs) = delete;
  public:
    const FreeModuleSpPRImpl* operator->() const; // allow member fns to be called
  };


  const FreeModuleSpPRImpl* FreeModuleSpPRPtr(const module& M)
  { return dynamic_cast<const FreeModuleSpPRImpl*>(M.myRawPtr()); }
  

  const FreeModuleSpPRImpl* FreeModuleSpPRPtr(const module& M, const ErrorContext& ErrCtx)
  {
    const FreeModuleSpPRImpl* ptr = FreeModuleSpPRPtr(M);
    if (ptr == nullptr) CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqFreeModule, ErrCtx);
    return ptr;
  }


  // const FreeModuleSpPRImpl* FreeModuleSpPRPtr(const module& M, const char* const FnName)
  // {
  //   const FreeModuleSpPRImpl* ptr = FreeModuleSpPRPtr(M);
  //   if (ptr == nullptr) CoCoA_THROW_ERROR(ERR::ReqFreeModule, FnName);
  //   return ptr;
  // }


  inline bool IsFreeModuleSpPR(const module& M)
  { return FreeModuleSpPRPtr(M) != nullptr; }


  inline FreeModuleSpPR::FreeModuleSpPR(const module& M):
      FreeModule(FreeModuleSpPRPtr(M, CoCoA_ERROR_CONTEXT))
  {}


  inline FreeModuleSpPR::FreeModuleSpPR(const FreeModuleSpPRImpl* MPtr):
    FreeModule(MPtr)
  {}


  inline const FreeModuleSpPRImpl* FreeModuleSpPR::operator->() const
  { return static_cast<const FreeModuleSpPRImpl*>(myRawPtr()); }


  //--- FreeModuleSpPRImpl

  FreeModuleSpPRImpl::FreeModuleSpPRImpl(const SparsePolyRing& P, const ModuleOrdering& O):
      FreeModuleImpl(P, NumComponents(O)),
      myModuleOrd(O),
      myOrdPtr()
  {
    if ( IsWDegPosnOrd(O) )  myOrdPtr.reset(new WDegPosnOrdImpl(shifts(O)));
    else if ( IsOrdPosn(O) ) myOrdPtr.reset(new OrdPosnImpl(shifts(O)));
    else  CoCoA_THROW_ERROR2(ERR::NYI, "this ModuleOrdering");
  }


  FreeModule NewFreeModuleSpPR(const SparsePolyRing& P, long NumCompts)
  {
    if (NumCompts < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    return FreeModule(new FreeModuleSpPRImpl(P, NewOrdPosn(ordering(PPM(P)),NumCompts)));
  }


  FreeModule NewFreeModuleSpPR(const SparsePolyRing& P, const ModuleOrdering& O)
  {return FreeModule(new FreeModuleSpPRImpl(P, O));}


  FreeModule NewFreeModule(const ring& P, const std::vector<degree>& shifts)
  {
    if (!IsSparsePolyRing(P))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    return FreeModule(new FreeModuleSpPRImpl(P, NewOrdPosn(ordering(PPM(P)), shifts)));
  }


  FreeModule NewFreeModule(const ring& P, const std::vector<degree>& shifts, const ModuleOrderingCtor& OrdCtor)
  {
    if (!IsSparsePolyRing(P))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    return FreeModule(new FreeModuleSpPRImpl(P, OrdCtor(ordering(PPM(P)), shifts)));
  }


  FreeModule NewFreeModule(const ring& P, long NumComponents, const ModuleOrderingCtor& OrdCtor)
  {
    if (!IsSparsePolyRing(P))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    return FreeModule(new FreeModuleSpPRImpl(P, OrdCtor(ordering(PPM(P)), NumComponents)));
  }


  const ModuleOrdering& FreeModuleSpPRImpl::myModuleOrdering() const
  { return myModuleOrd; }


  long FreeModuleSpPRImpl::myLPosn(const ModuleElem&  v) const
  { return myOrdPtr->myLPosn(v); }


  int FreeModuleSpPRImpl::myCmpLT(const ModuleElem&  v1, const ModuleElem&  v2) const
  { return myOrdPtr->myCmpLT(v1, v2); }


  void FreeModuleSpPRImpl::myWDeg(degree& d, const ModuleElem&  v) const
  { return myOrdPtr->myWDeg(d, v); }


  int FreeModuleSpPRImpl::myCmpWDeg(const ModuleElem&  v1, const ModuleElem&  v2) const
  { return myOrdPtr->myCmpWDeg(v1, v2); }


  //--- non-member fns

  const ModuleOrdering& ordering(const FreeModule& M)
  {
    if (!IsFreeModuleSpPR(M))  CoCoA_THROW_ERROR1(ERR::ReqModuleSpPR);
   return dynamic_cast<const FreeModuleSpPRImpl*>(M.myRawPtr())->myModuleOrdering();
  }


  const std::vector<degree>& shifts(const FreeModule& M)
  {
    if (!IsFreeModuleSpPR(M))  CoCoA_THROW_ERROR1(ERR::ReqModuleSpPR);
    return shifts(ordering(M));
  }


  //----------------------------------------------------------------------
  //--- CmpBase


  FreeModuleSpPRImpl::CmpBase::CmpBase(const std::vector<degree>& shifts):
    myShiftsValue(shifts)
  {}


  void FreeModuleSpPRImpl::CmpBase::myWDeg(degree& d, const ModuleElem& v) const
  {
    const long posn = myLPosn(v);
    d = wdeg(v[posn]) + myShiftsValue[posn];
  }


  int FreeModuleSpPRImpl::CmpBase::myCmpWDeg(const ModuleElem& v1, const ModuleElem& v2) const
  {
    return CmpWDeg(v1[LPosn(v1)], v2[LPosn(v2)]);
  }


  //----------------------------------------------------------------------
  //--- OrdPosnImpl

  FreeModuleSpPRImpl::OrdPosnImpl::OrdPosnImpl(const std::vector<degree>& shifts):
    CmpBase(shifts)
  {
    const long n = len(shifts);
    myComptsOrd = LongRange(0, n-1);
  }
  

  long FreeModuleSpPRImpl::OrdPosnImpl::myLPosn(const ModuleElem& v) const
  {
    long i;
    for ( i=0 ; i < len(myComptsOrd) ; ++i )
      if (!IsZero(v[myComptsOrd[i]]))
        break;
    CoCoA_ASSERT( i != len(myComptsOrd) );
    long MaxPos = myComptsOrd[i];
    degree MaxWDeg(wdeg(v[MaxPos])+myShiftsValue[MaxPos]);
    for ( ++i ; i < len(myComptsOrd) ; ++i )
      if ( !IsZero(v[myComptsOrd[i]]))
      {
        const long pos = myComptsOrd[i];
        int CmpWDegree = cmp(MaxWDeg, wdeg(v[pos])+myShiftsValue[pos]);
        if ( CmpWDegree==-1 || (CmpWDegree==0 && LPP(v[MaxPos])<LPP(v[pos])) )
        {
          MaxPos = pos;
          MaxWDeg = wdeg(v[MaxPos]) + myShiftsValue[MaxPos];
        }
      }
    return MaxPos;
  }


  int FreeModuleSpPRImpl::OrdPosnImpl::myCmpLT(const ModuleElem& v1, const ModuleElem& v2) const
  {
    const long LPosn1 = myLPosn(v1);
    const long LPosn2 = myLPosn(v2);
    int CmpFlag = cmp( wdeg(v1[LPosn1]) + myShiftsValue[LPosn1],
                       wdeg(v2[LPosn2]) + myShiftsValue[LPosn2] );
    if ( CmpFlag !=0 ) return CmpFlag;
    CmpFlag = cmp(LPP(v1[LPosn1]), LPP(v2[LPosn2]) );
    if ( CmpFlag !=0 ) return CmpFlag;
    return myCmpPosn(LPosn1, LPosn2);
  }


  int FreeModuleSpPRImpl::OrdPosnImpl::myCmpPosn(const long pos1, const long pos2) const
  {
    CoCoA_ASSERT( 0 <= pos1 && pos1 < len(myComptsOrd) );
    CoCoA_ASSERT( 0 <= pos2 && pos2 < len(myComptsOrd) );
    if ( pos1 == pos2 ) return 0;
    for ( long i=0 ; i < len(myComptsOrd) ; ++i )
    {
      if ( myComptsOrd[i]==pos1 ) return 1;
      if ( myComptsOrd[i]==pos2 ) return -1;
    }
    return 0;  // just to keep the compiler quiet
  }

  //----------------------------------------------------------------------
  //--- PosnOrdImpl

  FreeModuleSpPRImpl::PosnOrdImpl::PosnOrdImpl(const std::vector<degree>& shifts):
    CmpBase(shifts)
  {
    const long n = len(shifts);
    for ( long i=0 ; i<n ; ++i )
      myComptsOrd.push_back(i);
    // myComptsOrd = LongRange(0, len(shifts)-1);
  }


  long FreeModuleSpPRImpl::PosnOrdImpl::myLPosn(const ModuleElem& v) const
  {
    long i;
    for ( i=0 ; i < len(myComptsOrd) ; ++i )
      if (!IsZero(v[myComptsOrd[i]]))
        break;
    return i;
  }


  int FreeModuleSpPRImpl::PosnOrdImpl::myCmpLT(const ModuleElem& v1, const ModuleElem& v2) const
  {
    const long LPos1 = myLPosn(v1);
    const long LPos2 = myLPosn(v2);
    int CmpFlag = cmp( wdeg(v1[LPos1]) + myShiftsValue[LPos1],
                       wdeg(v2[LPos2]) + myShiftsValue[LPos2] );
    if ( CmpFlag !=0 ) return CmpFlag;
    CmpFlag = cmp(LPP(v1[LPos1]), LPP(v2[LPos2]) );
    if ( CmpFlag !=0 ) return CmpFlag;
    return myCmpPosn(LPos1, LPos2);
  }


  int FreeModuleSpPRImpl::PosnOrdImpl::myCmpPosn(const long pos1, const long pos2) const
  {
    CoCoA_ASSERT( 0 <= pos1 && pos1 < len(myComptsOrd) );
    CoCoA_ASSERT( 0 <= pos2 && pos2 < len(myComptsOrd) );
    if ( pos1 == pos2 ) return 0;
    for ( long i=0 ; i < len(myComptsOrd) ; ++i )
    {
      if ( myComptsOrd[i]==pos1 ) return 1;
      if ( myComptsOrd[i]==pos2 ) return -1;
    }
    return 0;  // just to keep the compiler quiet
  }

  //----------------------------------------------------------------------
  //--- WDegPosnOrdImpl

  FreeModuleSpPRImpl::WDegPosnOrdImpl::WDegPosnOrdImpl(const std::vector<degree>& shifts):
    CmpBase(shifts)
  {
    for ( long i=0 ; i < len(shifts) ; ++i )
      myComptsOrd.push_back(i);
    // myComptsOrd = LongRange(0, len(shifts)-1);
  }


  long FreeModuleSpPRImpl::WDegPosnOrdImpl::myLPosn(const ModuleElem& v) const
  {
    long i;
    for ( i=0 ; i < len(myComptsOrd) ; ++i )
      if (!IsZero(v[myComptsOrd[i]]))
        break;
    CoCoA_ASSERT( i != len(myComptsOrd) );
    long MaxPos = myComptsOrd[i];

    degree MaxWDeg(wdeg(v[MaxPos]) + myShiftsValue[MaxPos]);
    for ( ++i ; i < len(myComptsOrd) ; ++i )
      if ( !IsZero(v[myComptsOrd[i]]))
      {
        int CmpWDegree = cmp(MaxWDeg,
                             wdeg(v[myComptsOrd[i]])+myShiftsValue[myComptsOrd[i]]);
        if ( CmpWDegree==-1 )
        {
          MaxPos = myComptsOrd[i];
          MaxWDeg = wdeg(v[MaxPos]) + myShiftsValue[MaxPos];
        }
      }
    return MaxPos;
  }


  int FreeModuleSpPRImpl::WDegPosnOrdImpl::myCmpLT(const ModuleElem& v1, const ModuleElem& v2) const
  {
    const long LPos1 = myLPosn(v1);
    const long LPos2 = myLPosn(v2);
    int CmpFlag = cmp( wdeg(v1[LPos1]) + myShiftsValue[LPos1],
                       wdeg(v2[LPos2]) + myShiftsValue[LPos2] );
    if ( CmpFlag !=0 ) return CmpFlag;
    if ( LPos1 != LPos2 ) return ( LPos1>LPos2 ? 1 : -1 );
    return cmp(LPP(v1[LPos1]), LPP(v2[LPos2]));
  }


  int FreeModuleSpPRImpl::WDegPosnOrdImpl::myCmpPosn(const long pos1, const long pos2) const
  {
    CoCoA_ASSERT( 0 <= pos1 && pos1 < len(myComptsOrd) );
    CoCoA_ASSERT( 0 <= pos2 && pos2 < len(myComptsOrd) );
    if ( pos1 == pos2 ) return 0;
    for ( long i=0 ; i < len(myComptsOrd) ; ++i )
    {
      if ( myComptsOrd[i]==pos1 ) return 1;
      if ( myComptsOrd[i]==pos2 ) return -1;
    }
    return 0;  // just to keep the compiler quiet
  }

  //----------------------------------------------------------------------


  //  ConstRefPPMonoidElem LPP(const ModuleElem& v);

  long LPosn(const ModuleElem& v)
  {
    const module& M = owner(v);
    if ( !IsFreeModuleSpPR(M) )  CoCoA_THROW_ERROR1(ERR::ReqModuleSpPR);
    return FreeModuleSpPRPtr(M)->myLPosn(v);
  }


  ConstRefPPMonoidElem LPP(const ModuleElem& v)
  {
    const module& M = owner(v);
    if ( !IsFreeModuleSpPR(M) )  CoCoA_THROW_ERROR1(ERR::ReqModuleSpPR);
    return LPP(v[FreeModuleSpPRPtr(M)->myLPosn(v)]);
  }


  RingElemAlias LC(const ModuleElem& v)
  {
    const module& M = owner(v);
    if ( !IsFreeModuleSpPR(M) )  CoCoA_THROW_ERROR1(ERR::ReqModuleSpPR);
    return LC(v[FreeModuleSpPRPtr(M)->myLPosn(v)]);
  }


  degree wdeg(const ModuleElem& v)
  {
    const module& M = owner(v);
    if ( !IsFreeModuleSpPR(M) )  CoCoA_THROW_ERROR1(ERR::ReqModuleSpPR);
    degree ans(GradingDim(ordering(M)));
    FreeModuleSpPRPtr(M)->myWDeg(ans, v);
    return ans;
  }


  long StdDeg(const ModuleElem& v)
  {
    if ( !IsFreeModuleSpPR(owner(v)) )  CoCoA_THROW_ERROR1(ERR::ReqModuleSpPR);
    if (IsZero(v))  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    long PolyDegree = 0;
    const long n = NumCompts(v);
    for (long i=0; i<n; ++i)
      PolyDegree = max(PolyDegree, StdDeg(v[i]));
    return PolyDegree;
  }


  long deg(const ModuleElem& v)
  { return StdDeg(v); }
  

  int CmpWDeg(const ModuleElem& v1, const ModuleElem& v2)
  {
    const module& M1 = owner(v1);
    if ( M1 != owner(v2) )  CoCoA_THROW_ERROR1(ERR::MixedModules);
    if ( !IsFreeModuleSpPR(M1) )  CoCoA_THROW_ERROR1(ERR::ReqModuleSpPR);
    return FreeModuleSpPRPtr(M1)->myCmpWDeg(v1, v2);
  }


  ModuleElem homog(const ModuleElem& v, ConstRefRingElem h)
  {
    const FreeModule M(owner(v), CoCoA_ERROR_CONTEXT);
//    const module& M = owner(v);
//    if (!IsFreeModuleSpPR(M))
//      CoCoA_THROW_ERROR(ERR::ReqModuleSpPR, "ModuleElem(v, h)");
    if ( GradingDim(RingOf(M)) != 1 )  CoCoA_THROW_ERROR2(ERR::NYI, "GradingDim!=1");
    CoCoA_ASSERT( IsIndet(h) );
    CoCoA_ASSERT( IsOne(wdeg(h)[0]) );
    const BigInt d = wdeg(v)[0];
    const vector<ModuleElem>& e(gens(M));
    ModuleElem resv(M);
    for (long i=0; i<NumCompts(v); ++i)
      resv += e[i] * (homog(v[i],h) * power(h,d-wdeg(v)[0])); // and shifts ???
    return resv;
  }


  bool IsHomog(const ModuleElem& v)
  {
    if ( !IsFreeModuleSpPR(owner(v)) )  CoCoA_THROW_ERROR1(ERR::ReqModuleSpPR);
    if (IsZero(v)) return true;
    const FreeModule FM = owner(v);
    const std::vector<degree>& s = shifts(FM);
    const long FNZP = FirstNonZeroPosn(v);
    degree d(wdeg(v[FNZP])+s[FNZP]);
    for (long i=FNZP ; i<NumCompts(FM) ; ++i )
    {
      if (!IsZero(v[i]))
        if (!IsHomog(v[i]) || wdeg(v[i])+s[i] != d)
          return false;
    }
    return true;
  }


  FreeModule NewFreeModuleForSyz(const std::vector<RingElem>& L, const ErrorContext& ErrCtx)
  {
    if (L.empty()) CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqNonEmpty, ErrCtx);
    if (!IsSparsePolyRing(owner(L[0])))
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqElemPolyRing, ErrCtx);
    const SparsePolyRing P = owner(L[0]);
    std::vector<degree> sh;
    // transform(L.begin(), L.end(), back_inserter(sh),
    //           [](const RingElem& arg) { return wdeg(arg); });
    for (const auto& f:L)
      if (!IsZero(f)) sh.push_back(wdeg(f));
      else CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqNonZero, ErrCtx);
    return NewFreeModuleSpPR(P, NewWDegPosnOrd(ordering(PPM(P)), sh));
  }


  FreeModule NewFreeModuleForSyz(const std::vector<ModuleElem>& L, const ErrorContext& ErrCtx)
  {
    if (L.empty()) CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqNonEmpty, ErrCtx);
    if (!IsSparsePolyRing(RingOf(owner(L[0]))))
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqElemPolyRing, ErrCtx);
    const SparsePolyRing P = RingOf(owner(L[0]));
    std::vector<degree> sh;
    for (const auto& f:L)
      if (!IsZero(f)) sh.push_back(wdeg(f));
      else CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqNonZero, ErrCtx);
    return NewFreeModuleSpPR(P, NewWDegPosnOrd(ordering(PPM(P)), sh));
    //    return NewFreeModule(P, InputShifts, WDegPosnOrd);
  }


  FreeModule NewFreeModuleForSyz(const std::vector<RingElem>& L)
  { return NewFreeModuleForSyz(L, CoCoA_ERROR_CONTEXT); }

  FreeModule NewFreeModuleForSyz(const std::vector<ModuleElem>& L)
  { return NewFreeModuleForSyz(L, CoCoA_ERROR_CONTEXT); }



  //---------------------------------------------------------------------------

//   ConstRefRingElem FreeModuleElem::operator[](long pos) const
//   {
// //     if (!IsFGModule(owner(*this)))
// //       CoCoA_THROW_ERROR(ERR::NotFGModule, "FreeModuleElem[pos]");
//     if (pos < 0 || pos >= NumCompts(AsFGModule(owner(*this))))
//       CoCoA_THROW_ERROR(ERR::BadComptIndex, "FreeModuleElem[pos]");
//     return AsFGModule(owner(*this))->myCompt(raw(*this), pos);
//   }


//   RingElem& FreeModuleElem::operator[](long pos)
//   {
// //     if (!IsFGModule(owner(*this)))
// //       CoCoA_THROW_ERROR(ERR::NotFGModule, "FreeModuleElem[pos]");
//     if (pos < 0 || pos >= NumCompts(AsFGModule(owner(*this))))
//       CoCoA_THROW_ERROR(ERR::BadComptIndex, "FreeModuleElem[pos]");
//     return AsFGModule(owner(*this))->myCompt(raw(*this), pos);
//   }





}  // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/FreeModule.C,v 1.53 2024/08/01 16:55:28 bigatti Exp $
// $Log: FreeModule.C,v $
// Revision 1.53  2024/08/01 16:55:28  bigatti
// Summary: error codes: Req(Elem)SparsePolyRing;
// now using CoCoA_THROW_ERROR1/2, CoCoA_THROW_ERROR_WITH_CONTEXT2
//
// Revision 1.52  2024/07/22 15:43:57  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.51  2024/07/15 17:10:00  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.50  2024/07/03 14:46:49  bigatti
// Summary: error codes: Not... changed into ReqFGModule, ReqFreeModule, ReqModuleSpPR
//
// Revision 1.49  2024/07/02 15:31:12  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.48  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.47  2024/03/26 10:51:14  bigatti
// Summary: added NewFreeModuleForSyz with ErrorContext
//
// Revision 1.46  2024/03/15 18:57:32  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.45  2023/07/04 09:20:12  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.44  2023/03/13 19:50:26  abbott
// Summary: Deleted some cruft
//
// Revision 1.43  2023/03/09 22:31:58  abbott
// Summary: General tidying
//
// Revision 1.42  2022/02/18 15:38:24  abbott
// Summary: Using lambda fn to avoid clang "deprecated" warning
//
// Revision 1.41  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.40  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.39  2022/02/05 19:57:32  abbott
// Summary: Added CanonicalBasis (redmine 1235)
//
// Revision 1.38  2021/10/30 16:48:27  abbott
// Summary: Used keywords override and delete (redmine 1625 & 1627)
//
// Revision 1.37  2021/01/07 15:07:03  abbott
// Summary: Corrected copyright
//
// Revision 1.36  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.35  2020/02/11 16:56:40  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.34  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.33  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.32  2019/03/04 10:29:04  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.31  2018/08/01 15:06:55  bigatti
// -- removed commented code in NewFreeModuleForSyz (module)
//
// Revision 1.30  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.29  2018/05/18 16:40:27  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.28  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.27  2015/12/01 13:11:00  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.26  2014/07/30 14:05:02  abbott
// Summary: Changed BaseRing into RingOf, myBaseRing --> myRing
// Author: JAA
//
// Revision 1.25  2014/07/09 14:27:53  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.24  2014/07/07 12:22:41  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.23  2014/05/05 18:24:20  abbott
// Added missing include<iterator> for using std::back_inserter
//
// Revision 1.22  2014/04/30 16:07:21  abbott
// Summary: Removed pointless include
// Author: JAA
//
// Revision 1.21  2014/04/09 13:10:07  bigatti
// -- fixed shifts in SyzOfGens
//
// Revision 1.20  2014/04/02 05:46:55  bigatti
// -- added comment about use of LongRange (to be tested)
//
// Revision 1.19  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.18  2013/08/02 16:40:34  bigatti
// -- added NewFreeModule with shifts
//
// Revision 1.17  2013/08/02 14:41:14  bigatti
// -- added LC
// -- changed LPos --> LPosn
//
// Revision 1.16  2013/07/12 14:53:55  abbott
// Improved layout.
//
// Revision 1.15  2013/06/03 09:11:50  bigatti
// renamed ModuleTermOrdering into ModuleOrdering
//
// Revision 1.14  2013/05/28 07:00:32  bigatti
// -- fixed shifts
// -- more of pos --> posn
//
// Revision 1.13  2013/05/27 16:36:53  bigatti
// -- new constructor for FreeModule (with new implementation of orderings)
//
// Revision 1.12  2013/02/21 17:14:36  bigatti
// -- added NewFreeModuleForSyz
//
// Revision 1.11  2013/02/13 09:08:06  bigatti
// -- added homog(ModuleElem)
//
// Revision 1.10  2013/02/12 16:24:21  bigatti
// -- added deg, StdDeg
//
// Revision 1.9  2013/02/12 07:10:14  bigatti
// -- added FirstNonZero
// -- minor cleaning
//
// Revision 1.8  2012/10/12 12:38:18  abbott
// Removed element accessor (via operator[]) and non-const mem fn  ModuleBase::myCompt.
//
// Revision 1.7  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.6  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.5  2009/12/03 17:26:34  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.4  2008/04/21 12:57:04  abbott
// Changed a size_t arg into MachineInteger.  Added defn of myCompts memfn.
// Improved layout of some comments.
//
// Revision 1.3  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.15  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.14  2007/01/15 13:37:51  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.13  2006/12/07 17:27:29  cocoa
// -- for compilation with _Wextra: commented out useless comparisons
//
// Revision 1.12  2006/11/24 17:01:43  cocoa
// -- reorganized includes of header files
//
// Revision 1.11  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.10  2006/11/14 17:46:20  cocoa
// -- changed: reference counting in modules now uses SmartPtrIRC
//
// Revision 1.9  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.8  2006/10/19 13:56:24  cocoa
// Added #include<new> whenever needed (i.e. for placement new).
//
// Revision 1.7  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.6  2006/10/06 16:43:54  cocoa
// -- added PosnOrdImpl (warning: non-degree compatible ordering)
//
// Revision 1.5  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.4  2006/10/06 10:48:34  cocoa
// Removed #include references to GradedFreeModule.H
//
// Revision 1.3  2006/10/06 09:53:06  cocoa
// Moved definition of FreeModuleImpl out of the header file.
// Merged GradedFreeModule code into here so it can access FreeModuleImpl
// (which is no longer visible in the header file).
//
// Revision 1.2  2006/07/18 10:51:25  cocoa
// -- added: FirstNonZeroPos
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.8  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.7  2006/05/04 14:25:16  cocoa
// -- major cleaning of FreeModule: created GradedFreeModule and moved
//    some code around
//
// Revision 1.6  2006/04/21 14:56:33  cocoa
// Changed return type of myCompt member function: now it returns a
// ConstRefRingElem instead of a RingElem (i.e. a copy).
//
// Revision 1.5  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.4  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.3  2005/11/29 13:04:47  cocoa
// -- added "const" to myCompt argument
//
// Revision 1.2  2005/11/24 16:09:38  cocoa
// -- added operator[] for ModuleElem
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.4  2005/09/28 12:06:40  cocoa
// -- fixed "n undeclared" when compiled with Debug on (inside CoCoA_ASSERT)
//
// Revision 1.3  2005/09/28 11:50:34  cocoa
// -- new code for graded modules
//
// Revision 1.2  2005/06/22 14:42:16  cocoa
// Renamed MemPool data member to myMemMgr
// (seems more sensible than myMemory).
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
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
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.3  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.2  2004/01/28 15:56:49  cocoa
// "Old style" code, brought in alignment with the new coding conventions.
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
