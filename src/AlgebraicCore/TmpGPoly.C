//   Copyright (c)  2005-2007  John Abbott,  Anna M. Bigatti
//   Original author: 2005-2007  Massimo Caboara

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

#include "CoCoA/TmpGPoly.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/ModuleOrdering.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/TmpGPair.H"
#include "CoCoA/VectorOps.H" // just for debugging and statistics
#include "CoCoA/assert.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"

#include <algorithm>
using std::min;
using std::max;
#include <iostream>
using std::ostream;
using std::endl;
//using std::swap;
#include <iterator>
#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;

namespace CoCoA
{  


  //---------class GPoly-------------------------------

  // WARNING: is not possible to build the zero GPoly here.
  // change the ctors to allow this possibility
  GPoly::GPoly(ConstRefRingElem the_p,
               const GRingInfo& theGRI,
               long age):  // default age=0
    myLPPForDivwMask(theGRI.myPPM(), theGRI.myDivMaskRule()),
    myLPPForOrd(LPP(the_p)),
    myLCValue(LC(the_p)),
    myPolyValue(the_p),
    myGRingInfoValue(theGRI),
    myWDeg(wdeg(LPP(the_p))),
    mySugar(uninitialized)
  {
    myLPPForDivwMask = exponents(myLPPForOrd);
    IamActive = true;
    //    IamMinimalGen = false;
    myMinimalGenLevel = -1;
    myAge = age;
    myNumTerms = NumTerms(the_p);
    myComponent = theGRI.myComponent(myLPPForDiv());
   }//ctor


 // for DYN use, the LPP and LC are given to the GPoly
 GPoly::GPoly(ConstRefRingElem the_p,
              ConstRefPPMonoidElem theLPP,
              ConstRefRingElem theLC,
              const GRingInfo& theGRI,
              long age):  // default age=0
    myLPPForDivwMask(theGRI.myPPM(), theGRI.myDivMaskRule()),
    myLPPForOrd(theLPP),
    myLCValue(theLC),
    myPolyValue(the_p),
    myGRingInfoValue(theGRI),
    myWDeg(wdeg(theLPP)),
    mySugar(uninitialized)
  {
    myLPPForDivwMask = exponents(myLPPForOrd);
    IamActive = true;
    //    IamMinimalGen = false;
    myMinimalGenLevel = -1;
    myAge = age;
    myNumTerms = NumTerms(the_p);
    myComponent = theGRI.myComponent(myLPPForDiv());
   }//ctor

  // This ctor destroys the_p
  GPoly::GPoly(RingElem& the_p,
               ConstRefPPMonoidElem theLPP,
               ConstRefRingElem theLC,
               const GRingInfo& theGRI,
               const ClearMarker,//just for operator ariety
               long age):  // default age=0
    myLPPForDivwMask(theGRI.myPPM(), theGRI.myDivMaskRule()),
    myLPPForOrd(theLPP),
    myLCValue(theLC),
    myPolyValue(theGRI.myNewSPR()),
    myGRingInfoValue(theGRI),
    myWDeg(wdeg(theLPP)),
    mySugar(uninitialized)
  {
    myLPPForDivwMask = exponents(myLPPForOrd);
    swap(the_p,myPolyValue);
    IamActive = true;
    //    IamMinimalGen = false;
    myMinimalGenLevel = -1;
    myAge = age;
    myNumTerms = NumTerms(myPolyValue);
    myComponent = theGRI.myComponent(myLPPForDiv());
  }//ctor

// This ctor destroys the_p
  GPoly::GPoly(RingElem& the_p,
               const GRingInfo& theGRI,
               const ClearMarker,// just for operator ariety
               long age):  // default age=0
    myLPPForDivwMask(theGRI.myPPM(), theGRI.myDivMaskRule()),
    myLPPForOrd(LPP(the_p)),
    myLCValue(LC(the_p)),
    myPolyValue(theGRI.myNewSPR()),
    myGRingInfoValue(theGRI),
    myWDeg(wdeg(LPP(the_p))),
    mySugar(uninitialized)
  {
    myLPPForDivwMask = exponents(myLPPForOrd);
    swap(the_p,myPolyValue);
    IamActive = true;
    //    IamMinimalGen = false;
    myMinimalGenLevel = -1;
    myAge = age;
    myNumTerms = NumTerms(myPolyValue);
    myComponent = theGRI.myComponent(myLPPForDiv());
  }//ctor


  GPoly::GPoly(const GRingInfo& theGRI):
    myLPPForDivwMask(theGRI.myPPM(), theGRI.myDivMaskRule()),
    myLPPForOrd(PPM(theGRI.myNewSPR())),
    myLCValue(CoeffRing(theGRI.myNewSPR())),
    myPolyValue(theGRI.myNewSPR()),
    myGRingInfoValue(theGRI),
    myWDeg(GradingDim(theGRI.myNewSPR())),
    mySugar(uninitialized)
  {
    IamActive = true;
    //    IamMinimalGen = false;
    myMinimalGenLevel = -1;
    myNumTerms = 0;
    myComponent = 0;
    myAge = 0;//-1?
  }//ctor


  GPoly::GPoly(const GPoly& the_gp):
      myLPPForDivwMask(the_gp.myLPPForDivwMask),
      myLPPForOrd(the_gp.myLPPForOrd),
      myLCValue(LC(the_gp)),
      myPolyValue(the_gp.myPolyValue),
      myGRingInfoValue(the_gp.myGRingInfoValue),
      myWDeg(the_gp.myWDeg),
      mySugar(the_gp.mySugar)
  {
    IamActive = the_gp.IamActive;
    myMinimalGenLevel = the_gp.myMinimalGenLevel;
    myNumTerms = the_gp.myNumTerms;
    myAge = the_gp.myAge;
    myComponent = the_gp.myComponent;
 }//ctor


  GPoly& GPoly::operator=(const GPoly& the_gp) // ANNA: should it throw if not compatible???
  {
    CoCoA_ASSERT( AreCompatible(myGRingInfoValue,the_gp.myGRingInfoValue));
    myLPPForDivwMask = the_gp.myLPPForDivwMask;
    myLPPForOrd = the_gp.myLPPForOrd;
    myLCValue = the_gp.myLCValue;
    myPolyValue = the_gp.myPolyValue;
    myWDeg = the_gp.myWDeg;
    IamActive = the_gp.IamActive;
    myNumTerms = the_gp.myNumTerms;
    myWDeg = the_gp.myWDeg;
    myAge = the_gp.myAge;
    myComponent = the_gp.myComponent;
    mySugar = the_gp.mySugar;
    return *this;
  }//operator=

  void GPoly::AssignClear(GPoly& the_gp) // ANNA: should it throw if not compatible???
  {
    CoCoA_ASSERT( AreCompatible(myGRingInfoValue,the_gp.myGRingInfoValue));
    //    swap(myLPPForDivwMask, the_gp.myLPPForDivwMask);
    myLPPForDivwMask = the_gp.myLPPForDivwMask;
    swap(myLPPForOrd, the_gp.myLPPForOrd);
    myLCValue = the_gp.myLCValue;
    swap(myPolyValue,the_gp.myPolyValue);
    myWDeg = the_gp.myWDeg;
    IamActive = the_gp.IamActive;
    myNumTerms = the_gp.myNumTerms;
    myAge = the_gp.myAge;
    myComponent = the_gp.myComponent;
    mySugar = the_gp.mySugar;
    the_gp=GPoly(myGRingInfoValue);
  }//AssignClear


  bool GPoly::operator==(const GPoly& f)const
  {
    CoCoA_ASSERT(AreCompatible(myGRingInfoValue, f.myGRingInfoValue) );
    if (myPolyValue == f.myPolyValue) return true;
    return true;
  }//operator==


  bool GPoly::operator!=(const GPoly& f)const
  {
    CoCoA_ASSERT( AreCompatible(myGRingInfoValue, f.myGRingInfoValue) );
    if (myPolyValue == f.myPolyValue) return false;
    return true;
  }//operator!=


  bool IsZero(const GPoly& f) { return CoCoA::IsZero(f.myPolyValue); }


  ostream& operator<<(ostream& out, const GPoly& f)
  {
    if (!out) return out;  // short-cut for bad ostreams

    out<<"["<<f.myPolyValue
       <<", record["
       <<"IsActive="<<f.IamActive
       <<", NumTerms="<<f.myNumTerms
       <<", Deg="<<f.myWDeg
       <<", Sugar="<<f.mySugar
       <<", Age="<<f.myAge<<" "
       <<", LPP_Comp="<<f.myComponent
       <<", LPPForDiv="<<f.myLPPForDiv()
       <<", LPPForOrd="<<f.myLPPForOrd
       <<", LC="<<f.myLCValue
       <<"]]";
    return out;
  }

// here a GPoly is updated. used in reduction and SPoly production.
void GPoly::myUpdateLenLPPLCDegComp()
{
  myNumTerms = NumTerms(myPolyValue);
  if (IsZero(*this)) // the following things are effectively undefined
  {
    myLPPForDivwMask.myAssign(one(owner(myLPPForDiv())));
    AssignOne(myLPPForOrd);
    myLCValue = 0;
/// JAA BUG???    myWDeg = 0;      // ANNA delete???
    myComponent = 0; // ANNA delete ???
  }
  else
  {
    myLPPForOrd = LPP(myPoly());//DYN here the new LPP will be computed
    myLPPForDivwMask = exponents(myLPPForOrd);
    myLCValue=LC(myPoly());//DYN here the new LC will be computed
    myWDeg = wdeg(myLPPForOrd);
    myComponent = myGRingInfoValue.myComponent(myLPPForDiv());
  }
}//myUpdateLenLPPLCDegComp


  void GPoly::myInitializeSugar(const SugarDegree& s)
  {
    CoCoA_ASSERT(!IsInitialized(mySugar));
    mySugar = s;
  }
  

  void GPoly::myAssignSPoly(const GPair& the_gp, const long the_age)
  {
    IamActive = true;
    //    IamMinimalGen = false;
    myMinimalGenLevel = -1;
    if (the_gp.IsInputPoly())
      myPolyValue = poly(the_gp.myFirstGPoly());
    else
      myPolySetSPoly(the_gp.myFirstGPoly(), the_gp.mySecondGPoly());
    myUpdateLenLPPLCDegComp();
    myAge = the_age;
    // MAX: do these things only if necessary.
    mySugar = sugar(the_gp);
   }//myAssignSPoly

 /*
//???  This does not work, I don't understand why.
 void GPoly::myAppendClear(RingElem& p)
 {
   SparsePolyRing P=owner(*this);
clog << "operator+=: myPoly " <<myPoly<< endl;
clog << "operator+=: p " <<p<< endl;
   P->myAppendClear(raw(myPoly), raw(p));
clog << "operator+=: result " <<myPoly<< endl;
   myNumTerms = NumTerms(myPoly);
 }//_myAppendClear
 */

// TEMPORARY - Dangerous, does not adjust all the fields of *this
 void GPoly::myAppendClear(RingElem& p)
 {
   SparsePolyRingPtr(owner(*this))->myAppendClear(raw(myPolyValue), raw(p));
   myNumTerms = NumTerms(myPolyValue);
 }//myAppendClear

// TEMPORARY - Dangerous, does not adjust all the fields of *this
 void GPoly::myAppendClear(GPoly& p)
 {
   SparsePolyRingPtr(owner(*this))->myAppendClear(raw(myPolyValue), raw(p.myPolyValue));
   myNumTerms = NumTerms(myPolyValue);
 }//myAppendClear


  void  GPoly::MultiplyByPP(ConstRefPPMonoidElem MultPP)
  {
    // MultPP is in PPM(owner(GPoly))
    myLPPForOrd *= MultPP; // should we make it more efficient?
    myLPPForDivwMask = exponents(myLPPForOrd);
    SparsePolyRingPtr(owner(*this))->myMulByPP(raw(myPolyValue), raw(MultPP)); // (..) because of g++ parser bug
    myWDeg = wdeg(myLPPForOrd); // should we make it more efficient?
    mySugar->myMul(MultPP);
    // The other fields stay the same.
  }//MultiplyByPP


// This procedure should rely on the procedure for polys.
// When there are orderings, it should know if
// Ord=DRL, Var=last var, in which case may just return
// exponent(LPP(*this),DH_var_index)
  long max_common_wdeg(GPoly& f,long Var)
  {
    const SparsePolyRing P = owner(f);
    RingElem tmp(f.myPolyValue);
    long result=numeric_limits<long>::max();
    for (;!IsZero(tmp);)
    {
      result=min<long>(result,exponent(LPP(tmp),Var));
      P->myDeleteLM(raw(tmp));
    }
    return result;
  }//max_common_wdeg

//   void GRingInfo::WDegLessSVar(degree& res,ConstRefPPMonoidElem T)const// wdeg(T)-wdeg(saturating to the power it has in T)
//   {
//     PPMonoid TmpPPM=owner(T);
//     long LastVarIndex=NumIndets(TmpPPM)-1;
//     PPMonoidElem T1(TmpPPM->myOne());
//     vector<long> exps(TmpPPM->myNumIndets());
//     TmpPPM->myExponents(exps,raw(T));
//     TmpPPM->myMulIndetPower(raw(T1),LastVarIndex,exps[LastVarIndex]);
//     PPMonoidElem T2(TmpPPM);
//     TmpPPM->myDiv(raw(T2),raw(T),raw(T1));
//     res=wdeg(T2);
//   }//WDegLessSVar
  
//   degree GRingInfo::WDegLessSVarA(ConstRefPPMonoidElem T)const// wdeg(T)-wdeg(saturating to the power it has in T)
//   {
//     PPMonoid TmpPPM=owner(T);
//     long LastVarIndex=NumIndets(TmpPPM)-1;
//     PPMonoidElem T1(TmpPPM->myOne());
//     vector<long> exps(TmpPPM->myNumIndets());
//     TmpPPM->myExponents(exps,raw(T));
//     TmpPPM->myMulIndetPower(raw(T1),LastVarIndex,exps[LastVarIndex]);
//     PPMonoidElem T2(TmpPPM);
//     TmpPPM->myDiv(raw(T2),raw(T),raw(T1));
//     return wdeg(T2);
//   }//WDegLessSVar

//  degree GRingInfo::WDegLessSVarB(ConstRefPPMonoidElem T01,ConstRefPPMonoidElem T02)const// wdeg(T)-wdeg(saturating to the power it has in T)
//   {
//     PPMonoid TmpPPM=owner(T01);
//     PPMonoidElem T(TmpPPM);
//     TmpPPM->myDiv(raw(T),raw(T01),raw(T02));
//     long LastVarIndex=NumIndets(TmpPPM)-1;
//     PPMonoidElem T1(TmpPPM->myOne());
//     vector<long> exps(TmpPPM->myNumIndets());
//     TmpPPM->myExponents(exps,raw(T));
//     TmpPPM->myMulIndetPower(raw(T1),LastVarIndex,exps[LastVarIndex]);
//     PPMonoidElem T2(TmpPPM);
//     TmpPPM->myDiv(raw(T2),raw(T),raw(T1));
//     return wdeg(T2);
//   }//WDegLessSVar

// void GRingInfo::WDegLessSVarSmart(degree& res,
//                                   ConstRefPPMonoidElem T1, 
//                                   ConstRefPPMonoidElem T2)const
// {
//    PPMonoid TmpPPM=owner(T1);
//    long LastVarIndex=NumIndets(TmpPPM)-1;
//    res.mySetComponent(0,exponent(T1,LastVarIndex)-exponent(T2,LastVarIndex));//WARNING: this works only for st deg
// }//WDegLessSVarSmart

// void GRingInfo::WDegLessSVarFake(degree& res,
//                                  ConstRefPPMonoidElem T1, 
//                                  ConstRefPPMonoidElem T2)const
// {
//    PPMonoid TmpPPM=owner(T1);
//    long LastVarIndex=NumIndets(TmpPPM)-1;
//    res.mySetComponent(0,exponent(T1,LastVarIndex)-exponent(T2,LastVarIndex));//WARNING: this works only for st deg
// }//WDegLessSVarFake


// WARN : this function should rely on smart_dehomog for polys.
// Using that, smart_dehomog_DRL and smart_dehomog are few lines and equal.
  void GPoly::smart_dehomog_DRL(long DH_var_index)
  {
    long mc_deg=exponent(LPPForDiv(*this),DH_var_index);
    const SparsePolyRing P = owner(*this);
    RingElem result(P);
    RingElem tmp(P);
    RingElem H2Deg = P->myMonomial(raw(one(CoeffRing(P))),
                                   raw(IndetPower(PPM(P), DH_var_index, mc_deg)));
    if (mc_deg!=0)
    {
      for (;!IsZero(myPolyValue);)
      {
        P->myDivLM(raw(tmp),raw(myPolyValue),raw(H2Deg));
        P->myDeleteLM(raw(myPolyValue));
        P->myAppendClear(raw(result),raw(tmp));
      }
      swap(myPolyValue, result);
      myLPPForOrd /= IndetPower(PPM(P), DH_var_index, mc_deg);
      myLPPForDivwMask = exponents(myLPPForOrd);
      myWDeg = wdeg(myLPPForOrd);
    }
    // myComponent and myLC... stay the same
  }//smart_dehomog_DRL


  void GPoly::smart_dehomog(long DH_var_index)
  {
    long mc_deg=max_common_wdeg(*this,DH_var_index);
    const SparsePolyRing P = owner(myPolyValue);
    RingElem result(P);
    RingElem tmp(P);
    RingElem H2Deg = P->myMonomial(raw(one(CoeffRing(P))),
                                   raw(IndetPower(PPM(P), DH_var_index, mc_deg)));

    if (mc_deg!=0)
    {
      for (;!IsZero(myPolyValue);)
      {
        P->myDivLM(raw(tmp),raw(myPolyValue),raw(H2Deg));
        P->myDeleteLM(raw(myPolyValue));
        P->myAppendClear(raw(result),raw(tmp));
      }
      swap(myPolyValue, result);
      myLPPForOrd /= IndetPower(PPM(P), DH_var_index, mc_deg);
      myLPPForDivwMask = exponents(myLPPForOrd);
      myWDeg = wdeg(myLPPForOrd);
    }
    // myComponent and myRing... stay the same
  }//smart_dehomog


//******** ReductorData ******************************************************

  ReductorData::ReductorData(GPoly* p, const long p_component,  const long count):
      myLPPForDivwMask(LPPForDivwMask(*p))
  {
    myGPolyPtr=p;
    myKey=MakeKey(*p);
    myComponent=p_component;
    myCount = count;
    IamBorelUpdated = true;
    myIamNotToBeUsedValue=false;
  }


  ReductorData::ReductorData(const ReductorData& RD):
      myLPPForDivwMask(RD.myLPPForDivwMask)
  {
    myGPolyPtr=RD.myGPolyPtr;
    myKey=RD.myKey;
    myComponent=RD.myComponent;
    myCount = RD.myCount;
    IamBorelUpdated = RD.IamBorelUpdated;
    myIamNotToBeUsedValue=RD.myIamNotToBeUsedValue;
  }


  ostream& operator<<(ostream& out, const ReductorData& RD)
  {
    if (!out) return out;  // short-cut for bad ostreams

    out<<"Record[ Key=" << RD.myKey
       <<", age=" << age(*(RD.myGPolyPtr))
       <<", LPPForDiv=" << PP(RD.myLPPForDivwMask)
       <<", LPPForOrd=" << LPPForOrd(*(RD.myGPolyPtr))
       <<", MdCmp=" << RD.myComponent
       //<<", Used = " << RD.myCount  //TMP SAT
       <<", Upd = " << RD.IamBorelUpdated
       <<", ToBeUsd = " << RD.myIamNotToBeUsedValue
       <<", Wdeg=" << wdeg(*(RD.myGPolyPtr))
       <<", Sugar=" << sugar(*(RD.myGPolyPtr))
       <<", PtrPoly="<<*(RD.myGPolyPtr)
       <<"]";
    return out;
  }

  // The Key field controls the reductors ordering
  long MakeKey(const GPoly& gp)
  {
    //if (Len(gp)<255) return Len(gp);
    // return gp.myAge+255;
    return gp.myAge;
  }




//********* Reductors *****************************************************


  Reductors::Reductors(const GRingInfo& P, const UseBorelFlag UBR):
      myGRingInfoValue(P)
  {
    myReductors.reserve(10000);
    if (UBR == DontUseBorel)
      IhaveBorelReductorsFlag = false;
    else
    {
      IhaveBorelReductorsFlag = true;
      myBorelReductors.reserve(10000);
    }
  }


  Reductors::Reductors(const GRingInfo& P):
      myGRingInfoValue(P)
  {
    myReductors.reserve(10000);
    IhaveBorelReductorsFlag = false;
  }


  const PPMonoid& PPM(const Reductors& red)
  {
    return PPM(owner(red));
  }


  const SparsePolyRing& owner(const Reductors& red)
  {
    return red.myGRingInfo().myNewSPR();
  }

  void Reductors::Insert(GPoly* p, const long count)
  {
    myReductors.push_back(ReductorData(p, component(*p), count));
    //This is useless if  myKey is  Age.
    long N = len(myReductors)-1;
    for (long i=0;i!=N;++i)
      if (myReductors[N]<myReductors[i]) std::swap(myReductors[i],myReductors[N]);
    // ANNA : Borel algorithm
    if (IhaveBorelReductorsFlag)
    {
      myBorelGPolys.push_back(*p); // ANNA
      GPoly* f = &myBorelGPolys.back();
      myBorelReductors.push_back(ReductorData(f,component(*f), count));
    }
  }


  void Reductors::myStampaReductors(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    out << "TheREDUCTORS := " << myReductors << endl;

    if (IhaveBorelReductorsFlag)
      out << "The_BOREL_REDUCTORS := " << myBorelReductors << endl;
  }


   // Find the (unique and existing) ReductorData which ptr is equal to GPolyPtr.
   // Existence is guaranteed by the fact that this procedure is used ONLY for
   // the final interreduction.
   vector<ReductorData>::iterator Reductors::find(GPoly* GPolyPtr)
   {
     for (vector<ReductorData>::iterator it=myReductors.begin();
	                                 it!=myReductors.end();
					 ++it)
       if (it->myGetGPolyPtr()==GPolyPtr)
         return it;

     CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "Reductors::find");
     return myReductors.end();// for compilator's sake.
   }//find


/*Reduces the polys in G of degree Deg(f) with the poly f */
  void Reductors::myReduceTails(const GPoly& f)
  {
    for (std::vector<ReductorData>::reverse_iterator it=myReductors.rbegin();
         it!=myReductors.rend();
         ++it )
      (it->myGPolyPtr)->myReduceTail(f);
  }


/*Reduces the polys in G of degree Deg(f) with the poly f
  WARN G is supposed to be ordered by Deg                */
  void Reductors::OrderedInterreduce(const GPoly& f)
  {
    degree d(wdeg(f));

    for (std::vector<ReductorData>::reverse_iterator it=myReductors.rbegin();
         it!=myReductors.rend()&&(wdeg(*(it->myGPolyPtr))==d);
         ++it )
      (it->myGPolyPtr)->myReduceTail(f);
  }//OrderedInterreduce


/*Reduces the polys in G with the poly f*/
  void Reductors::SuperInterreduce(const GPoly& f)
  {
    //  unsigned  int d = deg(f);
    for (std::vector<ReductorData>::reverse_iterator it=myReductors.rbegin();
         it!=myReductors.rend();
         ++it )
      (it->myGPolyPtr)->myReduceTail(f);
  }//SuperInterreduce

  // clean the reductors keeping the same GRI
  void Reductors::myClear()
  {
    myReductors.clear();
    myBorelReductors.clear();
    myBorelGPolys.clear();
  }//clean

//   void swap(Reductors& R1,Reductors& R2)
//   {
//     // TODO NEW if (R1.myPolyRing!=R2.myPolyRing) return;
//     swap(R1.myReductors,R2.myReductors);
//   }


// This function prepares the Borel Reductors for the next degree
  void Reductors::myBorelReductorsUpdateInNextDegree()
  {
    for (ReductorData& R: myBorelReductors)
//    for (vector<ReductorData>::const_iterator it = myBorelReductors.begin();
//         it != myBorelReductors.end();
//         ++it)
      R.IamBorelUpdated = false;
  }


  // forward declaration of class GPoly needed.
  void monic(std::list<GPoly>& GPL)
  {
    if (GPL.empty())
      return;
    const SparsePolyRing P(owner(GPL));
    for (GPoly& g: GPL)
    {
      g.myPolyValue= monic(g.myPolyValue);
      g.myLCValue = one(P);
    }
  }//monic


  SparsePolyRing owner(const PolyList& thePL)
  {
    CoCoA_ASSERT(!thePL.empty());
    return owner(*thePL.begin());
  }//owner

  SparsePolyRing owner(const GPolyList& theGPL)
  {
    CoCoA_ASSERT(!theGPL.empty());
    return owner(*theGPL.begin());
  }//owner

//  SparsePolyRing owner(const  std::vector<RingElem>& thePV)
//   {
//     CoCoA_ASSERT(!thePV.empty());
//     return owner(*thePV.begin());
//   }//owner

  void PolyList2PolyVectorClear(PolyList& thePL,std::vector<RingElem>& thePV)
  {
    thePV.clear();
    if (thePL.empty())
      return;
    const SparsePolyRing P(owner(thePL));
    for (RingElem& g: thePL)
    {
      thePV.push_back(one(P));
      swap(thePV.back(),g);
    }
  }//PolyList2PolyVector

  void PolyVector2PolyListClear(std::vector<RingElem>& thePV,PolyList& thePL)
  {
    thePL.clear();
    if (thePV.empty())
      return;
    const SparsePolyRing P(owner(thePV));
    for (RingElem& g: thePV)
    {
      thePL.push_back(one(P));
      swap(thePL.back(),g);
    }
  }//PolyVector2PolyList

  // power(y_1y_2..y_k,d_1d_2..d_k)=y_1^d_1y_2^d_2..y_k^d_k
  void power(RingElem& theResult,
             const std::vector<RingElem>& theV,
             const degree& the_d)
  {
    CoCoA_ASSERT(len(theV)==GradingDim(the_d));
     CoCoA_ASSERT(GradingDim(owner(theResult))==GradingDim(the_d));
     const SparsePolyRing SPR=owner(theResult);
     theResult=one(SPR);
     if (theV.empty())
       return;
     for (long j=0; j < GradingDim(SPR); ++j)
       theResult*=power(theV[j],the_d[j]);
  }//power


}// end namespace cocoa

// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpGPoly.C,v 1.50 2024/09/04 15:42:19 bigatti Exp $
// $Log: TmpGPoly.C,v $
// Revision 1.50  2024/09/04 15:42:19  bigatti
// Summary: renamed Component into component
//
// Revision 1.49  2024/05/30 15:07:20  bigatti
// Summary: renamed interreduce into myReduceTails
//
// Revision 1.48  2024/05/03 12:17:08  abbott
// Summary: Removed cruft
//
// Revision 1.47  2023/03/14 22:35:22  abbott
// Summary: New iterator for loop syntax (redmine 1436)
//
// Revision 1.46  2023/03/13 21:06:40  abbott
// Summary: New iterator "for" syntax; some cruft removed where the "for" has already been changed.  Redmine 1346
//
// Revision 1.45  2022/10/11 15:39:35  bigatti
// Summary: using new interface of exponents returning vector
//
// Revision 1.44  2022/02/18 14:12:00  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.43  2022/02/14 14:32:48  bigatti
// Summary: added date for original author
//
// Revision 1.42  2022/02/11 09:49:35  abbott
// Summary: Updated copyright notices (redmine 855)
//
// Revision 1.41  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.40  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.39  2019/11/14 17:59:31  abbott
// Summary: Removed a little cruft
//
// Revision 1.38  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.37  2018/05/18 12:24:47  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.36  2018/05/17 15:52:29  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.35  2017/11/29 17:42:08  bigatti
// -- just spaces
//
// Revision 1.34  2017/11/27 08:47:13  bigatti
// -- just spaces in a comment
//
// Revision 1.33  2017/04/26 12:53:54  bigatti
// -- some cleaning: len --> NumTerms, f.IsActive() --> IsActive(f), ..
//
// Revision 1.32  2017/04/18 09:33:32  bigatti
// -- "age" now lower-case
//
// Revision 1.31  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.30  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.29  2014/07/07 13:05:25  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.28  2013/10/28 13:22:33  bigatti
// -- myMinimalGenLevel, IsInputPoly
//
// Revision 1.27  2013/06/12 08:50:03  bigatti
// -- added IamMinimalGen
//
// Revision 1.26  2013/06/03 09:11:50  bigatti
// renamed ModuleTermOrdering into ModuleOrdering
//
// Revision 1.25  2012/10/16 09:51:53  abbott
// Replaced  RefRingElem  by  RingElem&
//
// Revision 1.24  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.23  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.22  2011/03/11 17:38:55  bigatti
// -- changed  unsigned int --> long
//
// Revision 1.21  2011/03/11 15:33:30  bigatti
// -- changed size_t into long
//
// Revision 1.20  2011/03/11 12:09:57  bigatti
// -- changed size --> len
//
// Revision 1.19  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.18  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.17  2010/07/16 09:29:58  bigatti
// -- minor cleaning and coding conventions
//
// Revision 1.16  2010/05/14 09:53:09  bigatti
// -- removed empty ctor for SugarDegree
// -- added marker for SugarDegree(uninitialized)
// -- SugarDegree for GBasis input is initialized by myPrepareGBasis
//
// Revision 1.15  2010/03/23 14:43:07  bigatti
// -- class GRingInfo estracted from TmpGPoly
//
// Revision 1.14  2009/12/03 17:41:01  abbott
// Removed useless include of FGModule.H
//
// Revision 1.13  2009/10/27 17:15:14  bigatti
// -- fixed: using sugar(g)->myWSugar() instead of wsugar(g)
//
// Revision 1.12  2009/05/14 12:27:39  abbott
// Fixed typo (changed "ESS:" into "ERR::")
//
// Revision 1.11  2009/04/27 13:38:17  bigatti
// -- changed DetermineComputationType
// -- added temporary functions (DegLess..)
//
// Revision 1.10  2009/04/27 12:31:32  bigatti
// -- added   case SaturatingAlgNoDRL
//
// Revision 1.9  2009/02/09 14:22:36  bigatti
// -- changed std::swap(RingElem, RingElem) --> swap(RingElem, RingElem)
//
// Revision 1.8  2009/01/08 06:57:57  bigatti
// -- fixed smart_dehomog (myLPPForOrd was not updated)
//
// Revision 1.7  2008/09/19 13:33:42  bigatti
// -- added: Sat algorithm (M.Caboara)
//
// Revision 1.6  2008/09/16 15:03:42  bigatti
// -- added LPPForDiv
// -- changed LPP into LPPForOrd
//
// Revision 1.5  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.4  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.3  2007/12/04 14:27:06  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.18  2007/03/08 18:50:05  bigatti
// -- disabled function  bool GRingInfo::DetermineIfMyGradingIsPosPlus(const SparsePolyRing& theSPR)
//
// Revision 1.17  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.16  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.15  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.14  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.13  2006/12/22 11:37:40  cocoa
// Removed operator<< for PolyList and GPolyList.  A few consequential changes.
//
// Revision 1.12  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.11  2006/12/04 13:55:54  cocoa
// -- added: sugar for GradingDim > 0  (called wsugar)
//
// Revision 1.10  2006/11/27 14:24:02  cocoa
// -- removed ";" in debugging printing function
//
// Revision 1.9  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.8  2006/11/09 17:46:58  cocoa
// -- version 0.9712:
// --   IdealImpl moved to SparsePolyRing from concrete rings
// -- PolyList in GTypes is now vector<RingElem> (was list)
// -- "my" coding convention applied to DistrMPoly
//
// Revision 1.7  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.6  2006/10/06 16:32:06  cocoa
// -- changed: GPoly::SPoly --> GPoly::myAssignSPoly
// -- changed: Len(const GPoly&) --> len(const GPoly&)
// -- added:   poly(const GPoly&)
// -- added:   GPoly::myUpdateLenLPPDegComp()
// -- in reduce.C added functions for computing sugar during reduction
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
// Revision 1.3  2006/08/17 09:28:56  cocoa
// -- added: sugar
// -- added: controls on homogeneity
//
// Revision 1.2  2006/07/20 14:15:42  cocoa
// -- removed unused "wdeg"
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.20  2006/05/04 14:25:16  cocoa
// -- major cleaning of FreeModule: created GradedFreeModule and moved
//    some code around
//
// Revision 1.19  2006/04/27 15:06:58  cocoa
// -- removed (wrong) check on the degree in Reductors::interreduce(const GPoly& f)
//
// Revision 1.18  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.17  2006/04/21 16:46:14  cocoa
// -- minor changes by Max
//
// Revision 1.16  2006/04/11 14:29:48  cocoa
// -- minor cleaning
//
// Revision 1.15  2006/04/11 14:16:29  cocoa
// -- reorganization of fns between reduce,SparsePolyRing,GPoly
// -- added optional "len" argument to myAssignReset in ReductionCog
//
// Revision 1.14  2006/04/10 16:24:38  cocoa
// -- changed implementation of reduce by a single GPoly
//    (I didn't really do it last time..)
//
// Revision 1.13  2006/03/22 17:58:45  cocoa
// -- changed implementation of reduce by a single GPoly
//
// Revision 1.12  2006/03/17 18:17:16  cocoa
// -- changed: use of ReductionCog for reduction (major cleanup)
//
// Revision 1.11  2006/03/16 14:34:43  cocoa
// -- moved GlobalLogput out of my way
//
// Revision 1.10  2006/03/07 17:01:21  cocoa
// -- changed: operations on coefficients are now decided by
//    "myCoeffRingType" (need a better name for it!)
//
// Revision 1.9  2006/03/02 13:43:26  cocoa
// -- changed: myPoly --> myPolyValue
//
// Revision 1.8  2006/02/14 16:23:58  cocoa
// -- defined "operator<<" for vector<RingElem>& in ring.C/H
//
// Revision 1.7  2006/02/13 14:46:45  cocoa
// -- changes by Max
//
// Revision 1.6  2006/02/13 13:50:20  cocoa
// -- changes by Max (GRingInfo)
//
// Revision 1.5  2006/01/20 15:43:30  cocoa
// -- fixed: use of RefPPMonoidElem and ConstRefPPMonoidElem
//
// Revision 1.4  2006/01/18 15:58:20  cocoa
// -- new changes by Max
//
// Revision 1.3  2006/01/17 15:44:56  cocoa
// -- chamges by Max for operations with modules
//
// Revision 1.2  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/09/28 11:50:34  cocoa
// -- new code for graded modules
//
// Revision 1.2  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
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
// Revision 1.12  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.11  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.10  2004/11/11 13:37:33  cocoa
// -- change: cout --> GlobalLogput()
//
// Revision 1.9  2004/10/29 15:55:54  cocoa
// -- changed myLPP into myLPPwMask (PPMonoidElem --> PPWithMask)
// -- added DivMask::base arguments in GPoly constructors
//
// Revision 1.8  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.7  2004/06/16 16:13:41  cocoa
// Improved I/O facilities with knock-on changes
//
// Revision 1.6  2004/05/27 16:27:21  cocoa
// -- removed ";" at the end of function bodies (g++ gives error on them)
//
// Revision 1.5  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.4  2004/03/04 11:37:18  cocoa
// -- updated code for Borel reductors:
//    ReductorData fields: myGPoly->myGPolyPtr;  new: myCount, IamBorelUpdated
//    myBorelReductors is now in Reductors (was in GReductor)
//    Reductors: field: IhaveBorelReductors;  type: enum UseBorelMarker
//
// Revision 1.3  2003/10/08 14:27:33  cocoa
// New naming convention for ring member functions.
//
// Revision 1.2  2003/10/01 10:35:32  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.14  2003/09/22 17:00:04  bigatti
// - myDeg is now defined as deg(myPoly) instead of deg(myLPP)
//
// Revision 1.13  2003/06/23 17:10:23  abbott
// Minor cleaning prior to public release.
// Improved the include directives,
//
// Revision 1.12  2003/05/29 16:31:02  bigatti
// - change: RingSpecialIndex is now int
// - change: order of SPoly arguments
//
// Revision 1.11  2003/05/28 14:14:59  bigatti
// - new code for modules
//
// Revision 1.10  2003/05/14 16:52:37  bigatti
// - new ring/PPMonoid syntax
// - new functions for "BorelReductorsPolys" and saturating algorithm
// - added myPolyRing field
// - myDeg is now of type degree
//
// Revision 1.9  2002/11/15 17:29:56  bigatti
// - ASSERT --> CoCoA_ASSERT
//
// Revision 1.8  2002/09/19 17:17:44  bigatti
// - More general structure using PolyRing
//
// Revision 1.7  2002/05/13 11:45:42  bigatti
// - new data structure for "Reductors"
//
// Revision 1.6  2002/04/17 10:46:14  bigatti
// - new function: SPoly moved here from DMP
// - new field: myLPP
//
// Revision 1.5  2002/04/15 17:26:44  bigatti
// - Max's new code
//
// Revision 1.4  2002/04/09 14:06:20  bigatti
// - SPoly now takes a GPair as argument
//
// Revision 1.3  2002/01/31 17:19:55  bigatti
// - new function: ReduceTail
//
// Revision 1.2  2001/12/12 17:56:33  bigatti
// - new structure of reduction
//
// Revision 1.1  2001/12/05 12:56:02  bigatti
// Initial revision
//
