//   Copyright (c)  2005-2017 John Abbott, Anna M. Bigatti
//   Authors: 2005-2010 Massimo Caboara, 2010-2017 Anna M. Bigatti

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

#include "CoCoA/GBEnv.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DenseMatrix.H" // for DetermineIfMyGradingIsPosPlus
#include "CoCoA/FreeModule.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/RingQQ.H" // for DetermineIfMyGradingIsPosPlus
#include "CoCoA/VectorOps.H" // just for debugging and statistics
#include "CoCoA/assert.H"
#include "CoCoA/matrix.H" // for DetermineIfMyGradingIsPosPlus
#include "CoCoA/symbol.H"

using std::vector;
#include <limits>
using std::numeric_limits;
#include <algorithm>
using std::min;
using std::max;
#include <iostream>
using std::ostream;
using std::endl;
//using std::swap;



namespace CoCoA
{  

  // BUG BUG BUG remove 30000 from line below!!!
  /*static*/ const long GRingInfo::myMaxComponentIndex = min(30000ul,
                                                             min(static_cast<unsigned long>(numeric_limits<long>::max()-1),
                                                                 static_cast<unsigned long>(numeric_limits<SmallExponent_t>::max()-1))); // max num of compts -- depends on type SmallExponent_t



// Two utilities for GRingInfo ctors

  RingHom CreateNewP2OldPRingHom(const SparsePolyRing& theNewSPR,
                                 const SparsePolyRing& theOldSPR)
  {
    if (theNewSPR==theOldSPR)
      return  IdentityHom(theNewSPR);
    std::vector<RingElem> images;
    for (long i=0; i!=NumIndets(theOldSPR); ++i)
      images.push_back(indet(theOldSPR, i)); // x[i] |-> x[i]
    for (long i=0; i!=GradingDim(theNewSPR); ++i)
      images.push_back(one(theOldSPR)); // y[i] |-> 1
    images.push_back(one(theOldSPR)); // e |-> 1
    return  PolyRingHom(theNewSPR,theOldSPR,CoeffEmbeddingHom(theOldSPR),images);
  }//CreateNewP2OldPRingHom

  RingHom CreateOldP2NewPRingHom(const SparsePolyRing& theNewSPR,
                                 const SparsePolyRing& theOldSPR)
  {
     if (theNewSPR==theOldSPR)
       return  IdentityHom(theNewSPR);
     std::vector<RingElem> images;
     for (long i=0; i < NumIndets(theOldSPR); ++i)
       images.push_back(indet(theNewSPR, i));       // x[i] |-> x[i]
     return PolyRingHom(theOldSPR,theNewSPR,CoeffEmbeddingHom(theNewSPR),images);
  }//CreateOldP2NewPRingHom





  //-------------------------------------------------------
  //---------class GRingInfo-------------------------------
  //-------------------------------------------------------

  ComputationInputAndGradingType  DetermineComputationType(long GrDim,
                                                           const bool IsHomog,
                                                           const bool IsSatAlg)
  {
    if (IsSatAlg)
    { 
      //  if (!IsHomog) CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "DetermineComputationType");
      if (GrDim==0) return SaturatingAlgNoDRL;
      return SaturatingAlg;
    }
    if (GrDim==0) return NOGRADING;
    if (IsHomog) return HOMOG;
    return NONHOMOG_GRADING;
  }//DetermineComputationType
  

  // ----------------------------------------------------------------------
  // GRingInfo ctors

  void GRingInfo::myCtorAux(const SparsePolyRing& theNewSPR,
                            const bool IsHomog,
                            const bool IsSatAlg)
  {
    myInputAndGradingValue=DetermineComputationType(GradingDim(theNewSPR),
                                                    IsHomog,
                                                    IsSatAlg);
    myGradingPosPlusValue=DetermineIfMyGradingIsPosPlus(theNewSPR);
    mySetCoeffRingType(CoeffEncoding::Field);
  }
  

  GRingInfo::GRingInfo(const SparsePolyRing& theNewSPR,
                       const bool IsHomog,
                       const bool IsSatAlg,
                       const DivMaskRule& DivMaskR,
                       const CpuTimeLimit& CheckForTimeout):
    myNewSPRValue(theNewSPR),
    myOldSPRValue(theNewSPR),
    myPPMValue(NewPPMonoidEv(symbols(PPM(theNewSPR)), ordering(PPM(theNewSPR)))),
    myFreeModuleValue(NewFreeModule(myNewSPRValue,1)),
    myOutputFreeModuleValue(NewFreeModule(myNewSPRValue,1)),
    myNewP2OldPValue(IdentityHom(myNewSPRValue)),
    myOldP2NewPValue(IdentityHom(myNewSPRValue)),
    myDivMaskRuleValue(DivMaskR),
    IamModuleValue(false),
    myTimeoutChecker(CheckForTimeout)
  {
    myCtorAux(theNewSPR, IsHomog, IsSatAlg);
    myTimeoutChecker.myReset(IterationVariability::high);
  }// ctor GRingInfo


  GRingInfo::GRingInfo(const SparsePolyRing& theNewSPR,
                       const SparsePolyRing& theOldSPR,
                       const FreeModule& theFM,
                       const FreeModule& theOutputFM,
                       const bool IsHomog,
                       const bool IsSatAlg,
                       const DivMaskRule& DivMaskR,
                       const CpuTimeLimit& CheckForTimeout):
    myNewSPRValue(theNewSPR),
    myOldSPRValue(theOldSPR),
    myPPMValue(NewPPMonoidEv(symbols(PPM(theNewSPR)), ordering(PPM(theNewSPR)))),
    myFreeModuleValue(theFM),
    myOutputFreeModuleValue(theOutputFM),
    myNewP2OldPValue(CreateNewP2OldPRingHom(myNewSPRValue,myOldSPRValue)),
    myOldP2NewPValue(CreateOldP2NewPRingHom(myNewSPRValue,myOldSPRValue)),
    myDivMaskRuleValue(DivMaskR),
    IamModuleValue(theNewSPR!=theOldSPR),
    myTimeoutChecker(CheckForTimeout)
  {
    //    if (!IsField(CoeffRing(theOldSPR)))
    //      CoCoA_THROW_ERROR("coefficients are not in a field", "ComputeGBasis");
    std::vector<RingElem> Y; // The grading vars
    const std::vector<RingElem>& x = indets(myNewSPRValue);
    // Fill Y
    for (long i=0; i < GradingDim(myNewSPRValue); ++i)
       Y.push_back(x[i+NumIndets(theOldSPR)]);

    const std::vector<degree> S=shifts(myFreeModuleValue);
    RingElem tmp(myNewSPRValue);
    for (long i=0; i < NumCompts(myFreeModuleValue); ++i)
    {
      tmp=power(myE(),this->myComponent(i));
      for (long j=0; j < GradingDim(myNewSPRValue); ++j)
        tmp*=power(Y[j],S[i][j]);
      myEYValue.push_back(tmp);
    }
    myCtorAux(theNewSPR, IsHomog, IsSatAlg);
    myTimeoutChecker.myReset(IterationVariability::high);
  }// ctor GRingInfo
  

  GRingInfo::GRingInfo(const SparsePolyRing& theNewSPR,
                       const SparsePolyRing& theOldSPR,
                       const FreeModule& theOutputFM,
                       const bool IsHomog,
                       const bool IsSatAlg,
                       const DivMaskRule& DivMaskR,
                       const CpuTimeLimit& CheckForTimeout):
    myNewSPRValue(theNewSPR),
    myOldSPRValue(theOldSPR),
    myPPMValue(NewPPMonoidEv(symbols(PPM(theNewSPR)), ordering(PPM(theNewSPR)))),
    myFreeModuleValue(NewFreeModule(theNewSPR,1)),
    myOutputFreeModuleValue(theOutputFM),
    myNewP2OldPValue(CreateNewP2OldPRingHom(myNewSPRValue,myOldSPRValue)),
    myOldP2NewPValue(CreateOldP2NewPRingHom(myNewSPRValue,myOldSPRValue)),
    myDivMaskRuleValue(DivMaskR),
    IamModuleValue(theNewSPR!=theOldSPR),
    myTimeoutChecker(CheckForTimeout)
  {
    //    if (!IsField(CoeffRing(theOldSPR)))
    //      CoCoA_THROW_ERROR("coefficients are not in a field", "ComputeGBasis");
    myCtorAux(theNewSPR, IsHomog, IsSatAlg);
    myTimeoutChecker.myReset(IterationVariability::high);
  }// ctor GRingInfo


  GRingInfo::GRingInfo(const SparsePolyRing& theNewSPR,
                       const SparsePolyRing& theOldSPR,
                       const bool IsHomog,
                       const bool IsSatAlg,
                       const DivMaskRule& DivMaskR,
                       const CpuTimeLimit& CheckForTimeout):
    myNewSPRValue(theNewSPR),
    myOldSPRValue(theOldSPR),
    myPPMValue(NewPPMonoidEv(symbols(PPM(theNewSPR)), ordering(PPM(theNewSPR)))),
    myFreeModuleValue(NewFreeModule(theNewSPR,1)),
    myOutputFreeModuleValue(NewFreeModule(theNewSPR,1)),
    myNewP2OldPValue(CreateNewP2OldPRingHom(myNewSPRValue,myOldSPRValue)),
    myOldP2NewPValue(CreateOldP2NewPRingHom(myNewSPRValue,myOldSPRValue)),
    myDivMaskRuleValue(DivMaskR),
    IamModuleValue(theNewSPR!=theOldSPR),
    myTimeoutChecker(CheckForTimeout)
  {
    //    if (!IsField(CoeffRing(theOldSPR)))
    //      CoCoA_THROW_ERROR("coefficients are not in a field", "ComputeGBasis");
    myCtorAux(theNewSPR, IsHomog, IsSatAlg);
    myTimeoutChecker.myReset(IterationVariability::high);
  }// ctor GRingInfo

  // GRingInfo ctors
  // ----------------------------------------------------------------------

  void GRingInfo::mySetCoeffRingType(CoeffEncoding::type CT)
  { myCoeffRingTypeValue = CT; }


  bool GRingInfo::operator==(const GRingInfo& theGRI)const
  {
    return
      (myNewSPRValue    == theGRI.myNewSPRValue
       && myOldSPRValue == theGRI.myOldSPRValue
       && myPPMValue    == theGRI.myPPMValue
       && myOutputFreeModuleValue == theGRI.myOutputFreeModuleValue
       && myEYValue     == theGRI.myEYValue
       //&& // I want to do this, the == operator is not there
       //myDivMaskRuleValue==theGRI.myDivMaskRuleValue
       );
  }//operator==



long GRingInfo::myComponent(ConstRefPPMonoidElem T)const
{
  if (!IamModule()) return 0;// True Ring
  return exponent(T,ModuleVarIndex(myNewSPRValue));
}

long GRingInfo::myPhonyComponent(ConstRefPPMonoidElem T)const
{
  if (!IamModule()) return 0;// True Ring
  return myComponent(exponent(T,ModuleVarIndex(myNewSPRValue)));
}

RingElem GRingInfo::myY(const degree& the_d)const
{
   RingElem result(one(myNewSPR()));
   for (long j=0; j < GradingDim(myNewSPR()); ++j)
      result*=power(myY(j),the_d[j]);
   return result;
}//myY


  SugarDegree GRingInfo::myNewSugar(ConstRefRingElem f) const
  {
    switch (myInputAndGrading())
    {
    case HOMOG:            // ANNA: (w)graded + homogeneous
      return NewWSugarConst(f);
    case SaturatingAlg:    // SaturatingAlg
      return NewWSugarSat(f);
    case NONHOMOG_GRADING: // ANNA: (w)graded + non-homogeneous
      {
        if (/*module && */ IsMyGradingPosPlus())
        { // ANNA: should be implemented with proper weights
          int idx = ModuleVarIndex(myNewSPR());
          return NewStdSugarNoIdx(f, idx);
        }
        return NewWDeg1CompTmp(f);
      }
    case NOGRADING:        // ANNA: GradingDim = 0 --> StandardSugarAlgorithm
      //      if (/*module && */ IsMyGradingPosPlus())
      if (IamModule())
      {
        int idx = ModuleVarIndex(myNewSPR());
        return NewStdSugarNoIdx(f, idx);
      }
      return NewStdSugar(f);
    case SaturatingAlgNoDRL: // GradingDim = 0
      if (/*module && */ IsMyGradingPosPlus())
      {
        int idx = ModuleVarIndex(myNewSPR());
        return NewStdSugarNoIdxSat(f, idx);
      }
      return NewStdSugarSat(f);
    default: CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "GRingInfo::mySugar");
    }//switch
    return NewStdSugar(f); // just to keep the compiler quiet
  }


ostream& operator<<(ostream& out, const GRingInfo& theGRI)
{
  if (!out) return out;  // short-cut for bad ostreams
  out<<"the ring is "<<theGRI.myNewSPR()<<endl
     <<" the old ring is "<<theGRI.myOldSPR()<<endl
     //<<" Input Free Module "<<theGRI.myFreeModule()<<endl
     //<<" Output Free Module "<<theGRI.myOutputFreeModule()<<endl
     <<" IamModule "<<theGRI.IamModule()<<endl
     <<" myInputAndGrading = "<<theGRI.myInputAndGrading()<<endl
     <<" myGradingPosPlusValue = "<<theGRI.IsMyGradingPosPlus()<<endl
     <<" embedding grading "
     <<" EY=\n";
  for (const RingElem& f: theGRI.myEYValue)
  { out<<f<<endl; }
  out<<endl;
  return out;
}


long ModuleVarIndex(const SparsePolyRing& P)
{
  long tmp = NumIndets(P);
  if (tmp!=0)
    return tmp-1;
  else
    return tmp;
}//ModuleVarIndex


bool AreCompatible(const GRingInfo& theGRI1,const GRingInfo& theGRI2)
{
  if (theGRI1.myNewSPRValue==theGRI2.myNewSPRValue
        &&
      theGRI1.myOldSPRValue==theGRI2.myOldSPRValue
        &&
      theGRI1.myPPMValue==theGRI2.myPPMValue
      )
       //&& // I want to do this, the == operator is not there
         //theGRI1.myDivMaskRuleValue==theGRI2.myDivMaskRuleValue
    return true;
  else
    return false;
}//AreCompatible


// A member field?
std::vector<RingElem> GRingInfo::myY()const
{
  vector<RingElem> Y;
  for (long i=0; i < GradingDim(myNewSPRValue); ++i)
    Y.push_back(indet(myNewSPRValue,i+NumIndets(myOldSPRValue)));
  return Y;
}//myY()

 // Grdim>=2, order matrix first row is 0,..,0,1
 bool GRingInfo::DetermineIfMyGradingIsPosPlus(const SparsePolyRing& theSPR)
 {
   // This checks if indeed the order is a PosPlus.
   // Another option is to SET this field at the right time.
   // Slightly more efficient, but more risky.
   return false;
   if (GradingDim(theSPR)<1)
     return false;
   ConstMatrixView OrdM = OrdMat(ordering(PPM(theSPR)));
   // JAA 2015-11-30 line above replaces the two below
   // matrix OrdM(NewDenseMat(RingQQ(),NumIndets(theSPR),NumIndets(theSPR)));
   // PPM(theSPR)->myOrdering()->myOrdMatCopy(OrdM);
   for (long i=0; i<NumIndets(theSPR)-1; ++i)
     if (OrdM(0,i)!=0)
     {
       return false;
     }
   if (OrdM(0,NumIndets(theSPR)-1)!=1)
   {
     return false;
   }
   return true;
 }

}// end namespace cocoa

// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/GBEnv.C,v 1.35 2024/06/13 12:16:19 bigatti Exp $
// $Log: GBEnv.C,v $
// Revision 1.35  2024/06/13 12:16:19  bigatti
// Summary: using the same symbols as in input poly ring
//
// Revision 1.34  2024/06/09 20:10:20  abbott
// Summary: Renamed myPrepareForNewLoop to myReset (ProgressReporter & CpuTimeLimit); remine 1705
//
// Revision 1.33  2023/03/13 19:58:30  abbott
// Summary: Removed cruft and usless include dir
//
// Revision 1.32  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.31  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.30  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.29  2020/02/18 11:28:26  abbott
// Summary: redmine 1346: new for loop syntax
//
// Revision 1.28  2019/09/27 16:48:46  abbott
// Summary: Changed MaxComponentIndex to 30000 (so it does not look like an overflow problem)
//
// Revision 1.27  2018/06/27 08:50:39  abbott
// Summary: Revised to work with new CpuTimeLimit
//
// Revision 1.26  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.25  2018/05/17 15:44:36  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.24  2017/12/01 17:15:55  bigatti
// -- just better names for arguments
//
// Revision 1.23  2017/11/29 17:53:08  bigatti
// *** empty log message ***
//
// Revision 1.22  2017/11/27 08:46:36  bigatti
// -- mainly tidying spaces (reorganized operator==)
//
// Revision 1.21  2017/09/06 11:56:28  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.20  2017/04/28 16:03:49  abbott
// Summary: Corrected enum value in case in GRingInfo::myNewSugar
//
// Revision 1.19  2017/04/28 13:54:19  bigatti
// -- minor cleaning
// -- renamed AllGraded --> HOMOG, AllAffine --> NOGRADING
//
// Revision 1.18  2016/11/11 14:15:32  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.17  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.16  2015/07/25 15:39:39  abbott
// Summary: Corrected defn of myMaxComponentIndex
//
// Revision 1.15  2015/05/20 14:44:51  bigatti
// -- renamed AmIModule --> IamModule
//
// Revision 1.14  2015/05/19 07:23:19  abbott
// Summary: Added comment about removing 32748 when possible
// Author: JAA
//
// Revision 1.13  2015/05/15 17:28:27  bigatti
// -- found (temporary) fix for myMaxComponentIndex
//
// Revision 1.12  2015/05/15 16:08:59  bigatti
// -- investigating bug about myMaxComponentIndex
//
// Revision 1.11  2015/05/13 14:26:35  abbott
// Summary: MaxComponentIndex to max representable int
// Author: JAA
//
// Revision 1.10  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.9  2014/07/31 13:10:46  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.8  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.7  2012/04/02 15:28:23  bigatti
// -- ZZ --> QQ
//
// Revision 1.6  2012/02/10 10:26:40  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.5  2012/02/08 17:10:10  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2011/03/11 16:50:04  bigatti
// -- changes  unsigned int  --> long
//
// Revision 1.3  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.2  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.1  2010/03/23 14:40:55  bigatti
// -- first import
//
