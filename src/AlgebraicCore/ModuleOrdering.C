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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ModuleOrdering.C,v 1.15 2024/09/04 15:07:42 bigatti Exp $
// $Log: ModuleOrdering.C,v $
// Revision 1.15  2024/09/04 15:07:42  bigatti
// Summary: Now using CoCoA_THROW_ERROR1/2
//
// Revision 1.14  2024/07/02 15:31:12  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.13  2024/03/26 14:04:16  bigatti
// Summary: changed err::BadArg into NotNonNegative
//
// Revision 1.12  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.11  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.10  2021/10/30 11:53:48  abbott
// Summary: Used keyword override (redmine 1625); a little tidying too
//
// Revision 1.9  2020/06/17 15:49:24  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.8  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.7  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.6  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.5  2018/05/17 15:37:31  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.4  2016/11/11 14:15:32  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.3  2015/12/01 13:11:00  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.2  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.1  2013/06/03 08:51:58  bigatti
// -- was ModuleTermOrdering
//
// Revision 1.6  2013/05/31 10:31:09  abbott
// Moved NYI impl body of PosnOrdCtor::myCtor to *.C file to avoid multiple compiler
// warnings about unused parameter.
//
// Revision 1.5  2013/05/27 16:35:05  bigatti
// -- major reorganisation of the implementation, changed names
// ---- WDegPosTO --> WDegPosnOrd,  WDegTOPos --> OrdPosn,  PosWDegTO --> PosnOrd
//
// Revision 1.4  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.3  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.10  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.9  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.8  2007/01/18 15:34:04  cocoa
// -- changed namespace MTO into ModuleTermOrd
//
// Revision 1.7  2006/11/24 17:14:10  cocoa
// -- reorganized includes of header files
// -- SmartPtrIRC for reference counting
//
// Revision 1.6  2006/11/10 13:30:57  cocoa
// -- fixed: "const &" to PPOrdering arguments
// -- some more documentation
//
// Revision 1.5  2006/11/10 13:06:03  cocoa
// -- cleaned code
//
// Revision 1.4  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.3  2006/10/06 15:11:01  cocoa
// -- removed commented out code
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.8  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.7  2006/05/11 15:59:23  cocoa
// -- changed reference count is done using SmartPtrIRC
//
// Revision 1.6  2006/05/04 14:19:02  cocoa
// -- moved some code from .H to .C
//
// Revision 1.5  2006/04/28 11:32:16  cocoa
// -- moved concrete class definition from .H to .C
//
// Revision 1.4  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.3  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.2  2005/11/24 17:06:40  cocoa
// -- fixed IsWDegTOPos, IsWDegPosTO, IsWDegTOPermPos, IsWDegPermPosTO
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1  2005/09/28 11:50:34  cocoa
// -- new code for graded modules
//
