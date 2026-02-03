#ifndef CoCoA_GBEnv_H
#define CoCoA_GBEnv_H

//   Copyright (c)  2010  Anna Bigatti, Massimo Caboara

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

#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/DivMask.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SugarDegree.H"
#include "CoCoA/TmpGTypes.H"

//#include<list>
// using std::list; // for GPolyList; GPolyPtrList;
//#include<vector>
// using std::vector; // for GPolyPtrVector;

namespace CoCoA
{

  namespace CoeffEncoding
  {
    enum type {Field, FrFldOfGCDDomain};
  }


//-- class GRingInfo --------------------------------------------------*/


// Add fields for GradingIndex,ModulevarIndex.

long ModuleVarIndex(const SparsePolyRing&);

enum ComputationInputAndGradingType
{HOMOG, NOGRADING, NONHOMOG_GRADING, SaturatingAlg, SaturatingAlgNoDRL};



class GRingInfo
{

public:
  GRingInfo(const SparsePolyRing& theNewSPR,
            const bool HomogeneousInput,
	    const bool SaturatingAlgorithm,
            const DivMaskRule&,
            const CpuTimeLimit& CheckForTimeout);       // when you have polynomials only ie GB
  GRingInfo(const SparsePolyRing& theNewSPR,
            const SparsePolyRing& theOldSPR,
            const FreeModule& theFM,
            const FreeModule& theOutputFM,
            const bool HomogeneousInput,
 	    const bool SaturatingAlgorithm,
            const DivMaskRule&,
            const CpuTimeLimit& CheckForTimeout);       // when you have modules from theFM
  GRingInfo(const SparsePolyRing& theNewSPR,
            const SparsePolyRing& theOldSPR,
            const FreeModule& theOutputFM,
            const bool HomogeneousInput,
 	    const bool SaturatingAlgorithm,
            const DivMaskRule&,
            const CpuTimeLimit& CheckForTimeout);      // when you have module op on polynomials result vectorlist
  GRingInfo(const SparsePolyRing& theNewSPR,
            const SparsePolyRing& theOldSPR,
            const bool HomogeneousInput,
 	    const bool SaturatingAlgorithm,
            const DivMaskRule&,
            const CpuTimeLimit& CheckForTimeout);      // when you have module op on polynomials result polylist
  ~GRingInfo(){};
  void myCtorAux(const SparsePolyRing& theNewSPR,
                 const bool HomogeneousInput,
                 const bool SaturatingAlgorithm); ///< some assignments in common to all kind of GRingInfo

  const SparsePolyRing& myNewSPR() const {return myNewSPRValue;}
  const SparsePolyRing& myOldSPR() const {return myOldSPRValue;}
  const FreeModule& myFreeModule()  const {return myFreeModuleValue;}
  const FreeModule& myOutputFreeModule() const {return myOutputFreeModuleValue;}
  const RingHom& myNewP2OldP() const {return myNewP2OldPValue;}
  const RingHom& myOldP2NewP() const {return myOldP2NewPValue;}
  const DivMaskRule& myDivMaskRule() const {return myDivMaskRuleValue;}
  const PPMonoid& myPPM() const {return myPPMValue;}
  CoeffEncoding::type myCoeffRingType() const { return myCoeffRingTypeValue; }

  void mySetCoeffRingType(CoeffEncoding::type CT); ///< Field or FrFldOfGCDDomain
  long myComponent(const long i)const{return myMaxComponentIndex-i;}; //Complement to myMaxComponentIndex
  long myComponent(ConstRefPPMonoidElem T) const;///< the real component of a PP
  long myPhonyComponent(ConstRefPPMonoidElem T) const;///< the phony component of a PP
  bool operator==(const GRingInfo&) const;
  ///< Are equal myNewSPRValue and myDivMaskRuleValue only
  friend bool AreCompatible(const GRingInfo&,const GRingInfo&);
  friend std::ostream& operator<<(std::ostream& out, const GRingInfo& theGRI);
  bool IamModule()const{return IamModuleValue;};

  const RingElem& myE()const{return indet(myNewSPR(),ModuleVarIndex(myNewSPR()));};
  std::vector<RingElem> myEY()const{return myEYValue;};
  std::vector<RingElem> myY()const;
  SugarDegree myNewSugar(ConstRefRingElem f) const;
  RingElem myY(long i)const{return indet(myNewSPR(),i+NumIndets(myNewSPRValue)-1-GradingDim(myNewSPRValue));};
  RingElem myY(const degree& the_d)const;
  ConstRefRingElem myEY(long i)const{return myEYValue[i];};
  RingElem myE(long i)const{return power(myE(),myComponent(i));};
  //  bool amIGraded()const{return myEYValue.empty();};
  bool IsMyGradingPosPlus()const{return myGradingPosPlusValue;};
  void MySetComputationInputAndGradingType(const ComputationInputAndGradingType CIAGT)
        {myInputAndGradingValue=CIAGT;};
  //void SetMyGradingPosPlus(bool b){myGradingPosPlusValue=b;};
  ComputationInputAndGradingType myInputAndGrading()const{return myInputAndGradingValue;};
  friend ComputationInputAndGradingType  DetermineComputationType(long GrDim,
                                                                  const bool IsHomog,
                                                                  const bool IsSatAlg);
  bool DetermineIfMyGradingIsPosPlus(const SparsePolyRing&);
  void myCheckForTimeout(const char* const FnName) const { myTimeoutChecker(FnName); } 

private: // member fields
  const SparsePolyRing myNewSPRValue;///< embedding info
  const SparsePolyRing myOldSPRValue;
  const PPMonoid myPPMValue;  ///< the monoid for storing LPP's and lcm's for divisibility tests
  const FreeModule myFreeModuleValue;
  const FreeModule myOutputFreeModuleValue;
  RingHom myNewP2OldPValue;///< completely defined by myNewSPR,myOldSPRValue
  RingHom myOldP2NewPValue;
  const DivMaskRule myDivMaskRuleValue;
  const bool IamModuleValue;
  static const long myMaxComponentIndex; // max num of compts -- depends on type SmallExponent_t (see GBEnv.C)
  std::vector<RingElem> myEYValue;
  CoeffEncoding::type myCoeffRingTypeValue;
  ComputationInputAndGradingType myInputAndGradingValue;
  bool myGradingPosPlusValue;
  CpuTimeLimit myTimeoutChecker;
// FM and the two ring hom. Do we want a separate Module2RingEmbedding class?
//myReductionRule
//myLPPMonoid (?)

};// end class GRingInfo



}

#endif
