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
