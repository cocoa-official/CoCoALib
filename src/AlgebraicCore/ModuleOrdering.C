//   Copyright (c)  2005-2007,2011,2021  John Abbott, Anna Bigatti and Massimo Caboara

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

#include "CoCoA/ModuleOrdering.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/error.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/VectorOps.H"  // for template to output a vector

#include <iostream>
using std::ostream;
//#include <vector>
using std::vector;

namespace CoCoA
{

  ModuleOrderingBase::ModuleOrderingBase(const PPOrdering& PPO, const std::vector<degree>& shifts):
    IntrusiveReferenceCount(),
    myPPO(PPO)
  {
    const long n = len(shifts);
    for ( long i=0 ; i < n ; ++i )
      if ( GradingDim(PPO) != GradingDim(shifts[i]) )
        CoCoA_THROW_ERROR2(ERR::IncompatDims, "GradingDim PPOrd & shifts[i]");
    myShiftsValue = shifts;
  }


  ModuleOrderingBase::ModuleOrderingBase(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm):
    IntrusiveReferenceCount(),
    myPPO(PPO)
  {
    const long n = len(shifts);
    for ( long i=0 ; i < n ; ++i )
      if ( GradingDim(PPO) != GradingDim(shifts[i]) )
        CoCoA_THROW_ERROR2(ERR::IncompatDims, "GradingDim PPOrd & shifts[i]");
    myShiftsValue = shifts;
    if ( len(perm) != n )
      CoCoA_THROW_ERROR2(ERR::IncompatDims, "Shifts and permutation");
    CoCoA_THROW_ERROR2(ERR::NYI, "Checking entries of permutation");
    myPermutationValue = perm;
  }


  ModuleOrderingBase::~ModuleOrderingBase()
  {}


  //---------------------------------------------------------------------------


  const std::vector<degree>& ModuleOrderingBase::myShifts() const
  { return myShiftsValue; }


  const std::vector<long>& ModuleOrderingBase::myPerm() const
  { return myPermutationValue; }


  const PPOrdering& ModuleOrderingBase::myPPOrdering() const
  { return myPPO; }


  void ModuleOrderingBase::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    myOutputName(out);
    out << "(" << myPPO << ", " << myShiftsValue;
    if ( !myPermutationValue.empty() )
      out << ", " << myPermutationValue;
    out << ")";
  }


  void ModuleOrderingBase::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    // missing shifts and permutation
    OMOut->mySendApplyStart();
    myOutputName_OM(OMOut);
    OMOut << myPPO;
    OMOut->mySendApplyEnd();
  }


  //---------------------------------------------------------------------------
  //  (ex-inline) non-member functions

  std::ostream& operator<<(std::ostream& out, const ModuleOrdering& MTO)
  {
    if (!out) return out;  // short-cut for bad ostreams
    MTO->myOutputSelf(out);
    return out;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const ModuleOrdering& MTO)
  {
    MTO->myOutputSelf_OM(OMOut);
    return OMOut;
  }


  const std::vector<degree>& shifts(const ModuleOrdering& MTO)
  { return MTO->myShifts(); }


  const PPOrdering& ModPPOrdering(const ModuleOrdering& MTO)
  { return (MTO->myPPOrdering()); }


  long NumComponents(const ModuleOrdering& MTO)
  { return len(shifts(MTO)); }


  long GradingDim(const ModuleOrdering& MTO)
  { return GradingDim(ModPPOrdering(MTO)); }


  // ----------------- concrete classes ---------------------
  namespace ModuleOrd
  {

    class PosnOrdImpl: public ModuleOrderingBase
    {
    private: // pseudo-ctors
      friend ModuleOrdering CoCoA::NewPosnOrd(const PPOrdering& PPO, long NumComponents);
      friend ModuleOrdering CoCoA::NewPosnOrd(const PPOrdering& PPO, const std::vector<degree>& shifts);
    private:
      PosnOrdImpl(const PPOrdering& PPO, const std::vector<degree>& shifts);
    private: // disable copy ctor & assignment
      PosnOrdImpl(const PosnOrdImpl&) = delete;
      PosnOrdImpl& operator=(const PosnOrdImpl&) = delete;
    public:
      void myOutputName(std::ostream& out) const override;
      void myOutputName_OM(OpenMathOutput& OMOut) const override;
    };


    class OrdPosnImpl: public ModuleOrderingBase
    {
    private: // pseudo-ctors
      friend ModuleOrdering CoCoA::NewOrdPosn(const PPOrdering& PPO, long NumComponents);
      friend ModuleOrdering CoCoA::NewOrdPosn(const PPOrdering& PPO, const std::vector<degree>& shifts);
      friend ModuleOrdering CoCoA::NewOrdPosn(const PPOrdering& PPO, const std::vector<long>& perm);
      friend ModuleOrdering CoCoA::NewOrdPosn(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm);
    private:
      OrdPosnImpl(const PPOrdering& PPO, const std::vector<degree>& shifts);
      OrdPosnImpl(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm);
    private: // disable copy ctor & assignment
      OrdPosnImpl(const OrdPosnImpl&) = delete;
      OrdPosnImpl& operator=(const OrdPosnImpl&) = delete;
    public:
      void myOutputName(std::ostream& out) const override;
      void myOutputName_OM(OpenMathOutput& OMOut) const override;
    };


    class WDegPosnOrdImpl: public ModuleOrderingBase
    {
    private: // pseudo-ctors
      friend ModuleOrdering CoCoA::NewWDegPosnOrd(const PPOrdering& PPO, long NumComponents);
      friend ModuleOrdering CoCoA::NewWDegPosnOrd(const PPOrdering& PPO, const std::vector<degree>& shifts);
      friend ModuleOrdering CoCoA::NewWDegPosnOrd(const PPOrdering& PPO, const std::vector<long>& perm);
      friend ModuleOrdering CoCoA::NewWDegPosnOrd(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm);
    private:
      WDegPosnOrdImpl(const PPOrdering& PPO, const std::vector<degree>& shifts);
      WDegPosnOrdImpl(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm);
    private: // disable copy ctor & assignment
      WDegPosnOrdImpl(const WDegPosnOrdImpl&) = delete;
      WDegPosnOrdImpl& operator=(const WDegPosnOrdImpl&) = delete;
    public:
      void myOutputName(std::ostream& out) const override;
      void myOutputName_OM(OpenMathOutput& OMOut) const override;
    };


    //---------- PosnOrdImpl ----------------------------------------

    PosnOrdImpl::PosnOrdImpl(const PPOrdering& PPO, const std::vector<degree>& shifts):
      ModuleOrderingBase(PPO, shifts)
    {}


    void PosnOrdImpl::myOutputName(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "ModuleOrderingOrdPosn";
    }

    void PosnOrdImpl::myOutputName_OM(OpenMathOutput& OMOut) const
    {  OMOut << OpenMathSymbol("cocoa", "ModuleOrderingOrdPosn"); }


    //---------- OrdPosnImpl ----------------------------------------

    OrdPosnImpl::OrdPosnImpl(const PPOrdering& PPO, const std::vector<degree>& shifts):
      ModuleOrderingBase(PPO, shifts)
    {}


    OrdPosnImpl::OrdPosnImpl(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm):
      ModuleOrderingBase(PPO, shifts, perm)
    {}


    void OrdPosnImpl::myOutputName(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "ModuleOrderingOrdPosn";
    }

    void OrdPosnImpl::myOutputName_OM(OpenMathOutput& OMOut) const
    {  OMOut << OpenMathSymbol("cocoa", "ModuleOrderingOrdPosn"); }


    //---------- WDegPosnOrdImpl ----------------------------------------

    WDegPosnOrdImpl::WDegPosnOrdImpl(const PPOrdering& PPO, const std::vector<degree>& shifts):
      ModuleOrderingBase(PPO, shifts)
    {}


    WDegPosnOrdImpl::WDegPosnOrdImpl(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm):
      ModuleOrderingBase(PPO, shifts, perm)
    {}


    void WDegPosnOrdImpl::myOutputName(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "ModuleOrderingOrdPosn";
    }


    void WDegPosnOrdImpl::myOutputName_OM(OpenMathOutput& OMOut) const
    {  OMOut << OpenMathSymbol("cocoa", "ModuleOrderingOrdPosn"); }


    //------------------------------------------------------------//


  } // end of namespace ModuleOrd


  ModuleOrdering NewPosnOrd(const PPOrdering& PPO, long NumComponents)
  {
    if (NumComponents < 0)  CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "NumComponents");
    return ModuleOrdering(new ModuleOrd::PosnOrdImpl(PPO, vector<degree>(NumComponents, degree(GradingDim(PPO)))));
  }

  ModuleOrdering NewOrdPosn(const PPOrdering& PPO, long NumComponents)
  {
    if (NumComponents < 0)  CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "NumComponents");
    return ModuleOrdering(new ModuleOrd::OrdPosnImpl(PPO, vector<degree>(NumComponents, degree(GradingDim(PPO)))));
  }


  ModuleOrdering NewWDegPosnOrd(const PPOrdering& PPO, long NumComponents)
  {
    if (NumComponents < 0)  CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "NumComponents");
    return ModuleOrdering(new ModuleOrd::WDegPosnOrdImpl(PPO, vector<degree>(NumComponents, degree(GradingDim(PPO)))));
  }


  ModuleOrdering NewPosnOrd(const PPOrdering& PPO, const std::vector<degree>& shifts)
  { return ModuleOrdering(new ModuleOrd::PosnOrdImpl(PPO, shifts)); }


  ModuleOrdering NewOrdPosn(const PPOrdering& PPO, const std::vector<degree>& shifts)
  { return ModuleOrdering(new ModuleOrd::OrdPosnImpl(PPO, shifts)); }


  ModuleOrdering NewWDegPosnOrd(const PPOrdering& PPO, const std::vector<degree>& shifts)
  { return ModuleOrdering(new ModuleOrd::WDegPosnOrdImpl(PPO, shifts)); }


  ModuleOrdering NewOrdPosn(const PPOrdering& PPO, const std::vector<long>& perm)
  {
    return ModuleOrdering(new ModuleOrd::OrdPosnImpl(PPO, vector<degree>(len(perm), degree(GradingDim(PPO))), perm));
  }


  ModuleOrdering NewWDegPosnOrd(const PPOrdering& PPO, const std::vector<long>& perm)
  {
    return ModuleOrdering(new ModuleOrd::WDegPosnOrdImpl(PPO, vector<degree>(len(perm), degree(GradingDim(PPO))), perm));
  }


  ModuleOrdering NewOrdPosn(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm)
  { return ModuleOrdering(new ModuleOrd::OrdPosnImpl(PPO, shifts, perm)); }


  ModuleOrdering NewWDegPosnOrd(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm)
  { return ModuleOrdering(new ModuleOrd::WDegPosnOrdImpl(PPO, shifts, perm)); }


  bool IsPosnOrd(const ModuleOrdering& MOrd)
  {
    if (dynamic_cast<const ModuleOrd::PosnOrdImpl*>(MOrd.myRawPtr())) return true;
    // must decide whether the matrix is PosnWDeg..., possibly in disguise
    return false;
  }


  bool IsOrdPosn(const ModuleOrdering& MOrd)
  {
    if (dynamic_cast<const ModuleOrd::OrdPosnImpl*>(MOrd.myRawPtr())) return true;
    // must decide whether the matrix is WDeg..., possibly in disguise
    return false;
  }


  bool IsWDegPosnOrd(const ModuleOrdering& MOrd)
  {
    if (dynamic_cast<const ModuleOrd::WDegPosnOrdImpl*>(MOrd.myRawPtr())) return true;
    // must decide whether the matrix is WDeg..., possibly in disguise
    return false;
  }


  // STOPGAP Placeholder defn
  ModuleOrdering PosnOrdCtor::operator()(const PPOrdering& PPO, const std::vector<degree>& shifts) const
  { CoCoA_THROW_ERROR2(ERR::NYI, "PosnOrdCtor with shifts"); return operator()(PPO,0); }


  //----- declaration of ordering ctors ---------------------------
  OrdPosnCtor OrdPosn;
  PosnOrdCtor PosnOrd;
  WDegPosnOrdCtor WDegPosnOrd;
  //----- declaration of ordering ctors ---------------------------

} // end of namespace CoCoA
