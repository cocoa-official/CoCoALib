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
    myComponent = theGRI.myCompt_work(myLPPForDiv());
   }//ctor


  GPoly::GPoly(const GRingInfo& theGRI):
    myLPPForDivwMask(theGRI.myPPM(), theGRI.myDivMaskRule()),
    myLPPForOrd(PPM(theGRI.myP_work())),
    myLCValue(CoeffRing(theGRI.myP_work())),
    myPolyValue(theGRI.myP_work()),
    myGRingInfoValue(theGRI),
    myWDeg(GradingDim(theGRI.myP_work())),
    mySugar(uninitialized)
  {
    IamActive = true;
    //    IamMinimalGen = false;
    myMinimalGenLevel = -1;
    myNumTerms = 0;
    myComponent = 0;
    myAge = 0;//-1?
  }//ctor


  GPoly& GPoly::operator=(const GPoly& the_gp) // ANNA: should it throw if not compatible???
  {
// CoCoA_ASSERT( AreCompatible(myGRingInfoValue,the_gp.myGRingInfoValue));
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

  // void GPoly::AssignClear(GPoly& the_gp) // ANNA: should it throw if not compatible???
  // {
  //   CoCoA_ASSERT( AreCompatible(myGRingInfoValue,the_gp.myGRingInfoValue));
  //   //    swap(myLPPForDivwMask, the_gp.myLPPForDivwMask);
  //   myLPPForDivwMask = the_gp.myLPPForDivwMask;
  //   swap(myLPPForOrd, the_gp.myLPPForOrd);
  //   myLCValue = the_gp.myLCValue;
  //   swap(myPolyValue,the_gp.myPolyValue);
  //   myWDeg = the_gp.myWDeg;
  //   IamActive = the_gp.IamActive;
  //   myNumTerms = the_gp.myNumTerms;
  //   myAge = the_gp.myAge;
  //   myComponent = the_gp.myComponent;
  //   mySugar = the_gp.mySugar;
  //   the_gp=GPoly(myGRingInfoValue);
  // }//AssignClear


  // bool GPoly::operator==(const GPoly& f)const
  // {
  //   CoCoA_ASSERT(AreCompatible(myGRingInfoValue, f.myGRingInfoValue) );
  //   if (myPolyValue == f.myPolyValue) return true;
  //   return true;
  // }//operator==


  // bool GPoly::operator!=(const GPoly& f)const
  // {
  //   CoCoA_ASSERT( AreCompatible(myGRingInfoValue, f.myGRingInfoValue) );
  //   if (myPolyValue == f.myPolyValue) return false;
  //   return true;
  // }//operator!=


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
    myComponent = myGRingInfoValue.myCompt_work(myLPPForDiv());
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


// TEMPORARY - Dangerous, does not adjust all the fields of *this
 void GPoly::myAppendClear(RingElem& p)
 {
   SparsePolyRingPtr(myGRingInfoValue.myP_work())->myAppendClear(raw(myPolyValue), raw(p));
   myNumTerms = NumTerms(myPolyValue);
 }//myAppendClear

// TEMPORARY - Dangerous, does not adjust all the fields of *this
 void GPoly::myAppendClear(GPoly& p)
 {
   SparsePolyRingPtr(myGRingInfoValue.myP_work())->myAppendClear(raw(myPolyValue), raw(p.myPolyValue));
   myNumTerms = NumTerms(myPolyValue);
 }//myAppendClear


  namespace { // anonymous --------------------

// This procedure should rely on the procedure for polys.
// When there are orderings, it should know if
// Ord=DRL, Var=last var, in which case may just return
// exponent(LPP(*this),DH_var_index)
    // IMPLEMENT USING  SparsePolyIter <<----------------------
    // Check what it does and ADD TEST FOR "GBasisSelfSatCore;"
  long max_common_wdeg(GPoly& f,long Var)
  {
    const SparsePolyRing P = f.myGRingInfo().myP_work();
    RingElem tmp(poly(f));
    long result=numeric_limits<long>::max();
    for (;!IsZero(tmp);)
    {
      result=min<long>(result,exponent(LPP(tmp),Var));
      P->myDeleteLM(raw(tmp));
    }
    return result;
  }

  } // namespace { // anonymous --------------------


// WARN : this function should rely on smart_dehomog for polys.
// Using that, smart_dehomog_DRL and smart_dehomog are few lines and equal.
  void GPoly::smart_dehomog_DRL(long DH_var_index)
  {
    long mc_deg=exponent(LPPForDiv(*this),DH_var_index);
    const SparsePolyRing P = myGRingInfoValue.myP_work();
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
    const SparsePolyRing P = myGRingInfoValue.myP_work();
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
    myIamNotToBeUsedValue=false;
  }


  ReductorData::ReductorData(const ReductorData& RD)//:
    //      myLPPForDivwMask(RD.myLPPForDivwMask)
  {
    myGPolyPtr=RD.myGPolyPtr;
    myKey=RD.myKey;
    //    myComponent=RD.myComponent;
    myCount = RD.myCount;
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


  Reductors::Reductors(const GRingInfo& P, const CpuTimeLimit& CheckForTimeout):
    myGRingInfoValue(P),
    myTimeoutChecker(CheckForTimeout)
  {
    myReductors.reserve(10000);
    myTimeoutChecker.myReset(IterationVariability::high);
  }


  const PPMonoid& PPM(const Reductors& red)
  { return PPM(red.myGRingInfo().myP_work()); }


  void Reductors::myInsert(GPoly* p, const long count)
  {
    myReductors.push_back(ReductorData(p, component(*p), count));
    //This is useless if  myKey is  Age.
    long N = len(myReductors)-1;
    for (long i=0;i!=N;++i)
      if (myReductors[N]<myReductors[i]) std::swap(myReductors[i],myReductors[N]);
  }


  // Find the (unique and existing) ReductorData whose ptr is GPolyPtr.
  // Used ONLY for the final interreduction (guaranteed to succeed)
  vector<ReductorData>::iterator Reductors::myFind(GPoly* GPolyPtr)
  {
    for (auto it=myReductors.begin(); it!=myReductors.end(); ++it)
      if (it->myGetGPolyPtr()==GPolyPtr)  return it;
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    return myReductors.end(); // just to keep the compiler quiet
  }


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
  { myReductors.clear(); }


}// end namespace cocoa
