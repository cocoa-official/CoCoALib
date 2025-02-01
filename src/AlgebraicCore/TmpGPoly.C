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

  ReductorData::ReductorData(GPoly* p, const long p_component,  const long count)//:
    //      myLPPForDivwMask(LPPForDivwMask(*p))
  {
    myGPolyPtr=p;
    myKey=MakeKey(*p);
    //    myComponent=p_component;
    myCount = count;
    //    IamBorelUpdated = true;
    myIamNotToBeUsedValue=false;
  }


  ReductorData::ReductorData(const ReductorData& RD)//:
    //      myLPPForDivwMask(RD.myLPPForDivwMask)
  {
    myGPolyPtr=RD.myGPolyPtr;
    myKey=RD.myKey;
    //    myComponent=RD.myComponent;
    myCount = RD.myCount;
    //    IamBorelUpdated = RD.IamBorelUpdated;
    myIamNotToBeUsedValue=RD.myIamNotToBeUsedValue;
  }


  ostream& operator<<(ostream& out, const ReductorData& RD)
  {
    if (!out) return out;  // short-cut for bad ostreams

    out<<"Record[ Key=" << RD.myKey
       <<", age=" << age(*(RD.myGPolyPtr))
      //       <<", LPPForDiv=" << PP(RD.myLPPForDivwMask)
       <<", LPPForDiv=" << LPPForDiv(*(RD.myGPolyPtr))
       <<", LPPForOrd=" << LPPForOrd(*(RD.myGPolyPtr))
      //       <<", MdCmp=" << RD.myComponent
       //<<", Used = " << RD.myCount  //TMP SAT
      //       <<", Upd = " << RD.IamBorelUpdated
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


  // Reductors::Reductors(const GRingInfo& P, const UseBorelFlag UBR):
  //     myGRingInfoValue(P)
  // {
  //   myReductors.reserve(10000);
  //   if (UBR == DontUseBorel)
  //     IhaveBorelReductorsFlag = false;
  //   else
  //   {
  //     IhaveBorelReductorsFlag = true;
  //     myBorelReductors.reserve(10000);
  //   }
  // }


  Reductors::Reductors(const GRingInfo& P):
      myGRingInfoValue(P)
  {
    myReductors.reserve(10000);
    //    IhaveBorelReductorsFlag = false;
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
    // if (IhaveBorelReductorsFlag)
    // {
    //   myBorelGPolys.push_back(*p); // ANNA
    //   GPoly* f = &myBorelGPolys.back();
    //   myBorelReductors.push_back(ReductorData(f, component(*f), count));
    // }
  }


  void Reductors::myStampaReductors(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    out << "TheREDUCTORS := " << myReductors << endl;

    // if (IhaveBorelReductorsFlag)
    //   out << "The_BOREL_REDUCTORS := " << myBorelReductors << endl;
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
    // myBorelReductors.clear();
    // myBorelGPolys.clear();
  }//clean

//   void swap(Reductors& R1,Reductors& R2)
//   {
//     // TODO NEW if (R1.myPolyRing!=R2.myPolyRing) return;
//     swap(R1.myReductors,R2.myReductors);
//   }


// This function prepares the Borel Reductors for the next degree
//   void Reductors::myBorelReductorsUpdateInNextDegree()
//   {
//     for (ReductorData& R: myBorelReductors)
// //    for (vector<ReductorData>::const_iterator it = myBorelReductors.begin();
// //         it != myBorelReductors.end();
// //         ++it)
//       R.IamBorelUpdated = false;
//   }


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
