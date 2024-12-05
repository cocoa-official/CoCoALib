//   Copyright (c)  2005  John Abbott and Anna M. Bigatti

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


// Implementation file for the class submodule

#include "CoCoA/QuotientModule.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/FGModule.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/ring.H"

#include <vector>
using std::vector;
#include <iostream>
using std::ostream;


namespace CoCoA
{

  class QuotientModuleImpl: public FGModuleBase
  {
    // Two typedefs to save typing.
    typedef ModuleBase::RawPtr RawPtr;
    typedef const ModuleBase::RawPtr& ConstRawPtr;

  public:
    QuotientModuleImpl(const FGModule& Mnumer, const FGModule& Mdenom);
    long myNumCompts() const override;
    const ring& myRing() const override;
    const FreeModule& myAmbientFreeModule() const override;
    const std::vector<ModuleElem>& myGens() const override;
    const std::vector<ModuleElem>& myMinGens(const CpuTimeLimit& CheckForTimeout) const override;
    const std::vector<ModuleElem>& myTidyGens(const CpuTimeLimit& CheckForTimeout) const override;

    const ModuleElem& myZero() const override;
    void myNew(RawPtr& rawv) const override;
    void myNew(RawPtr& rawv, ConstRawPtr rawt) const override;
    void myDelete(RawPtr& rawv) const override;                                       // destroys v (incl all resources)
    void mySwap(RawPtr& rawv, RawPtr& raww) const override;                           // swap(v, w);
    void myAssign(RawPtr& rawlhs, ConstRawPtr rawv) const override;                   // lhs = v;
    ConstRefRingElem myCompt(const RawPtr& rawv, long pos) const override;            ///< v[pos]
    void myNegate(RawPtr& rawlhs, ConstRawPtr rawv) const override;                   // lhs = -v
    void myAdd(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const override;    // lhs = v+w;
    void mySub(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const override;    // lhs = v-w;

    void myMul(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const override; // lhs = r*v;
    void myDiv(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const override; // lhs = (1/r)*v;
    void myOutput(std::ostream& out, ConstRawPtr rawv) const override;                // out << v
    void myOutputSelf(std::ostream& out) const override;                              // out << M
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawv) const override;            // OMOut << v
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                          // OMOut << M
    bool myIsZero(ConstRawPtr rawv) const override;                                   // v == 0
//???    bool IsZeroAddMul(RawPtr& rawlhs, RingElemConstRawPtr rawy, ConstRawPtr rawz) const;  // lhs += y*z, result says whether lhs == 0.
    bool myIsEqual(ConstRawPtr, ConstRawPtr) const override;

  private: // data members
    const FreeModule myM;
    const FGModule myNumer;
    const FGModule myDenom;
    std::vector<ModuleElem> myGensArray;
    mutable bool myTidyGensIsValid;
    mutable std::vector<ModuleElem> myTidyGensArray;
//???    std::vector<ModuleElem>& ComputeTidyGens() const;
  };



  QuotientModuleImpl::QuotientModuleImpl(const FGModule& Mnumer, const FGModule& Mdenom):
      myM(AmbientFreeModule(Mnumer)),
      myNumer(Mnumer),
      myDenom(Mdenom)
  {
    if (AmbientFreeModule(Mnumer) != AmbientFreeModule(Mdenom))
      CoCoA_THROW_ERROR("Modules reside in different ambient free modules", "QuotientModule ctor");
    //??????? MISSING CODE  inherit gens of Mnumer directly??
    myRefCountZero();
  }


  long QuotientModuleImpl::myNumCompts() const
  {
    return NumCompts(myM);
  }


  const ring& QuotientModuleImpl::myRing() const
  {
    return RingOf(myM);
  }


  const FreeModule& QuotientModuleImpl::myAmbientFreeModule() const
  {
    return myM;
  }


  const vector<ModuleElem>& QuotientModuleImpl::myGens() const
  {
    return myGensArray;
  }


  const vector<ModuleElem>& QuotientModuleImpl::myMinGens(const CpuTimeLimit& /*CheckForTimeout*/) const
  {
//    if (!IsSparsePolyRing(myRing()))
//      CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, "QuotientModuleImpl::myMinGens()");
    CoCoA_THROW_ERROR(ERR::NYI, "QuotientModuleImpl::myMinGens");
    return myGensArray;
  }


  const vector<ModuleElem>& QuotientModuleImpl::myTidyGens(const CpuTimeLimit& /*CheckForTimeout*/) const
  {
    if (!myTidyGensIsValid)
    {
      CoCoA_THROW_ERROR(ERR::NYI, "QuotientModuleImpl::myTidyGens");
    }
    return myGensArray;
  }


  const ModuleElem& QuotientModuleImpl::myZero() const
  {
    return zero(myM);
  }


  void QuotientModuleImpl::myNew(RawPtr& rawv) const
  {
    myM->myNew(rawv);
  }


  void QuotientModuleImpl::myNew(RawPtr& rawv, ConstRawPtr rawcopy) const
  {
    myM->myNew(rawv, rawcopy);
  }


  void QuotientModuleImpl::myDelete(RawPtr& rawv) const
  {
    myM->myDelete(rawv);
  }


  void QuotientModuleImpl::mySwap(RawPtr& rawv, RawPtr& raww) const
  {
    myM->mySwap(rawv, raww);
  }


  void QuotientModuleImpl::myAssign(RawPtr& rawlhs, ConstRawPtr rawv) const
  {
    myM->myAssign(rawlhs, rawv);
  }


  ConstRefRingElem QuotientModuleImpl::myCompt(const RawPtr& rawv, long pos) const
  {
    CoCoA_ASSERT(0 <= pos && pos < myNumCompts());
    return myM->myCompt(rawv, pos);
  }

  void QuotientModuleImpl::myNegate(RawPtr& rawlhs, ConstRawPtr rawv) const
  {
    myM->myNegate(rawlhs, rawv);
  }


  void QuotientModuleImpl::myAdd(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const
  {
    myM->myAdd(rawlhs, rawv, raww);
    //??? reduce module Mdenom
  }


  void QuotientModuleImpl::mySub(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const
  {
    myM->mySub(rawlhs, rawv, raww);
    //??? reduce module Mdenom
  }


  void QuotientModuleImpl::myMul(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const
  {
    myM->myMul(rawlhs, rawx, rawv);
    //??? reduce module Mdenom
  }


  void QuotientModuleImpl::myDiv(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const
  {
    myM->myDiv(rawlhs, rawx, rawv);
    //??? reduce module Mdenom
  }


  void QuotientModuleImpl::myOutput(ostream& out, ConstRawPtr rawv) const
  {
    if (!out) return;  // short-cut for bad ostreams
    myM->myOutput(out, rawv);
  }


  void QuotientModuleImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "QuotientModule(" << myNumer << ", " << myDenom << ")";
  }


  void QuotientModuleImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawv) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "ModuleElement"); // BUG: what should this OMSymbol be???
    OMOut << myM;
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "list"); // BUG: what should this OMSymbol be???
    OMOut << myNumCompts();
    myM->myOutput_OM(OMOut, rawv); // BUG: this should be a "naked" output???
    OMOut->mySendApplyEnd();
    OMOut->mySendApplyEnd();
  }


  void QuotientModuleImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "QuotientModule"); // BUG: what should this OMSymbol be???
    OMOut << myNumer << myDenom;
    OMOut->mySendApplyEnd();
  }


  bool QuotientModuleImpl::myIsZero(ConstRawPtr rawv) const
  {
    return myM->myIsZero(rawv);
  }


//???    bool IsZeroAddMul(RawPtr& rawlhs, RingElemConstRawPtr rawy, ConstRawPtr rawz) const;  // lhs += y*z, result says whether lhs == 0.


  bool QuotientModuleImpl::myIsEqual(ConstRawPtr rawv, ConstRawPtr raww) const
  {
    return myM->myIsEqual(rawv, raww);
  }


  FGModule NewQuotientModule(const FGModule& Mnumer, const FGModule& Mdenom)
  {
    return FGModule(new QuotientModuleImpl(Mnumer, Mdenom));
  }


}  // end of namespace CoCoA
