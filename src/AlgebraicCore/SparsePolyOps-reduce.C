//   Copyright (c)  2005-2017  John Abbott, Anna M. Bigatti
//   Author:  2005-2017  Anna M. Bigatti

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
#include "CoCoA/ReductionCog.H"
#include "CoCoA/SugarDegree.H"
#include "CoCoA/VectorOps.H" // for debugging only
#include "CoCoA/assert.H"
#include "CoCoA/interrupt.H"
//#include "CoCoA/verbosity.H"

// #include <vector>
using std::vector;

namespace CoCoA
{

// Horrible hack, very temporary - profiling only
degree HereForProfilingOnlyWDeg(ConstRefPPMonoidElem cofactor1)
{return wdeg(cofactor1);}



  //-------- headers ----------------------------------------

  // ANNA: where should these go???
  void reduce(ReductionCog& F, const Reductors& v);
  void reduce(ReductionCog& F, const GPoly& g);
  ReductionCog ChooseReductionCogPoly(const GRingInfo& GRI);
  ReductionCog ChooseReductionCogGeobucket(const GRingInfo& GRI);

  //---------------------------------------------------------

  bool ReductorData::myBorelUpdate(ConstRefPPMonoidElem pp, const Reductors& theReductors)
  {
    if (myCount < 100) return false;

    GPoly* f = myGPolyPtr;
    const PPMonoid thePPM(PPM(*f));
    const long BorelIndetIndex = NumIndets(owner(*f))-1;
    PPMonoidElem quot(thePPM);

    thePPM->myDiv(raw(quot), raw(pp), raw(PP(myLPPForDivwMask)));
    if (quot != IndetPower(thePPM, BorelIndetIndex, exponent(quot, BorelIndetIndex)))
      return false;

    f->MultiplyByPP(quot);
    myLPPForDivwMask.myAssign(pp);
    IamBorelUpdated = true;
    f->myReduceTail(theReductors);
    return true;
  }


  inline GPoly* FindReducer(const PPWithMask& pm,
                            const long comp,
                            const Reductors& theReductors)
  {
    for (ReductorData& R: theReductors.myBorelReductors)
      if ( comp == R.myComponent && IsDivisibleFast(pm, R.myLPPForDivwMask))
        if (R.IamBorelUpdated || R.myBorelUpdate(PP(pm), theReductors) )
        {
          ++(R.myCount);
          return R.myGPolyPtr;
        }
    for (const ReductorData& R: theReductors.myReductors)
    {
      if (comp == R.myComponent && IsDivisibleFast(pm, R.myLPPForDivwMask) && !R.IamNotToBeUsed())
      {
        ++(R.myCount);
        return R.myGPolyPtr;
      }
    }
    return nullptr;
  }


  GPoly* FindReducer(const GPoly& g, const Reductors& theReductors)
  {
    if ( IsZero(g) ) return nullptr;
    GRingInfo GRI(theReductors.myGRingInfo());
    return FindReducer(LPPForDivwMask(g), GRI.myComponent(LPPForDiv(g)), theReductors);
  }


  //---- ReductionCog code --------------------------------------------------

  inline GPoly* FindReducer(const ReductionCog& F, const Reductors& theReductors)
  {
    if ( IsActiveZero(F) ) return nullptr;
    GRingInfo GRI(theReductors.myGRingInfo());
    //    const PPWithMask pm(ActiveLPP(F), GRI.myDivMaskRule());
    PPWithMask pm(GRI.myPPM(), GRI.myDivMaskRule());
    pm = exponents(ActiveLPP(F));
    return FindReducer(pm, GRI.myComponent(PP(pm)), theReductors);
  }


//   inline int FindReducerIndex(const ReductionCog& F, const vector<RingElem>& v)
//   {
//     if ( IsActiveZero(F) ) return -1;
//     return FindReducerIndex(ActiveLPP(F), v);
//   }


  void ReduceActiveLM(ReductionCog& F, const Reductors& v)
  {
    GPoly* g;
    while ( (g = FindReducer(F, v)) != nullptr )
    {
      //      std::cout << "ReduceActiveLM no sugar: *g is " << *g << std::endl;
      F->myReduce(poly(*g), NumTerms(*g));
    }
  }


  void ReduceActiveLM(ReductionCog& F, SugarDegree& s, const Reductors& v)
  {
    GPoly* g;
    while ( (g = FindReducer(F, v)) != nullptr )
    {
      //      std::cout << "ReduceActiveLM s: *g is " << *g << std::endl;
      CoCoA_ASSERT( !IsZero(*g));
      CheckForInterrupt("ReduceActiveLM");
      v.myGRingInfo().myCheckForTimeout("ReduceActiveLM");
      //      std::cout << "ReduceActiveLM before s->myUpdate" << std::endl;
      s->myUpdate(F, *g);
      //      std::cout << "ReduceActiveLM after s->myUpdate" << std::endl;
      F->myReduce(poly(*g), NumTerms(*g));
    }//while
  }//ReduceActiveLM


  void reduce(ReductionCog& F, const Reductors& v)
  {
    //    std::cout << "reduce: F is " << F << std::endl;
    ReduceActiveLM(F, v);
    //    std::cout << "after ReduceActiveLM " << F << std::endl;
    while ( !IsActiveZero(F) )
    {
      //      std::cout << "reduce: F is " << F << std::endl;
      F->myMoveToNextLM();
      ReduceActiveLM(F, v);
    }
  }


  void reduce(ReductionCog& F, SugarDegree& s, const Reductors& v)
  {
    //    std::cout << "--called reduce with sugar--" << std::endl;
    //    std::cout << "reduce: F is " << F << std::endl;
    ReduceActiveLM(F, s, v);
    while ( !IsActiveZero(F) )
    {
      //      std::cout << "reduce: F is " << F << std::endl;
      F->myMoveToNextLM();
      ReduceActiveLM(F, s, v);
    }
  }

  void reduce(ReductionCog& F, SugarDegree& s, const GPoly& g)
  {
    if ( IsActiveZero(F) || ActiveLPP(F) < LPPForOrd(g) ) return;
    const GRingInfo& GRI(g.myGRingInfo());
    const PPWithMask& PMg(LPPForDivwMask(g));
    const long Componentg = GRI.myComponent(PP(PMg));
    //    PPWithMask PMF(ActiveLPP(F), GRI.myDivMaskRule());
    PPWithMask PMF(GRI.myPPM(), GRI.myDivMaskRule());
    vector<long> expv;
    PMF = exponents(expv, ActiveLPP(F));
    long ComponentF = GRI.myComponent(PP(PMF));
    while ( !IsActiveZero(F) )
    {
      if ( (ComponentF==Componentg && IsDivisibleFast(PMF, PMg) ) )
      {
        s->myUpdate(F, g);
        F->myReduce(poly(g), NumTerms(g));// MAX: here the real work is done
      }
      else
        F->myMoveToNextLM();      
      if ( IsActiveZero(F) || ActiveLPP(F) < LPPForOrd(g) ) return;
      //      PMF.myAssign(ActiveLPP(F));
      PMF = exponents(expv, ActiveLPP(F));;
      ComponentF = GRI.myComponent(PP(PMF));
    }//while
  }//reduce


  void reduce(ReductionCog& F, const GPoly& g)
  {
    // SugarDegree s(sugar(g)); // not used --> use WSugarConst
    SugarDegree s(NewWSugarConst(zero(owner(g))));
    reduce(F, s, g);
  }


  //---- ReductionCog code end ----------------------------------------------


  //-------- ChooseReductionCog... ----------------------------------------

  //ANNA where should these go???
  ReductionCog ChooseReductionCogPoly(const GRingInfo& GRI)
  {
    if ( GRI.myCoeffRingType() == CoeffEncoding::Field )
      return NewRedCogPolyField(GRI.myNewSPR());
    else if ( GRI.myCoeffRingType() == CoeffEncoding::FrFldOfGCDDomain )
      return NewRedCogPolyGCD(GRI.myNewSPR());
    else CoCoA_THROW_ERROR("Don't know what to do with these coefficients", "ChooseReductionCog");
    return ReductionCog(nullptr);  // just to keep the compiler quiet
  }


  ReductionCog ChooseReductionCogGeobucket(const GRingInfo& GRI)
  {
    if ( GRI.myCoeffRingType() == CoeffEncoding::Field )
      return NewRedCogGeobucketField(GRI.myNewSPR());
    else if ( GRI.myCoeffRingType() == CoeffEncoding::FrFldOfGCDDomain )
      return NewRedCogGeobucketGCD(GRI.myNewSPR());
    else CoCoA_THROW_ERROR("Don't know what to do with these coefficients", "ChooseReductionCog");
    return ReductionCog(nullptr);  // just to keep the compiler quiet
  }


  //-------- GPoly member functions ----------------------------------------

  void GPoly::myPolySetSPoly(const GPoly& f, const GPoly& g)
  {
    myPolyValue = poly(f);
    (owner(f))->myMulByPP(raw(myPolyValue), raw(colon(LPPForOrd(g), LPPForOrd(f))));
    ReductionCog F = ChooseReductionCogPoly(myGRingInfo());
    F->myAssignReset(myPolyValue, f.myNumTerms);
    F->myReduce(poly(g), NumTerms(g));
    F->myRelease(myPolyValue);
    // sugar will be set as the sugar of the originating pair
  }


  void GPoly::myReduceTail(const GPoly& g)
  {
    CoCoA_ASSERT( !IsZero(*this) );
    if ( LPPForOrd(*this) < LPPForOrd(g) ) return;
    // geobucket because reduce(F, g) might make many reductions by g
    ReductionCog F = ChooseReductionCogGeobucket(myGRingInfo());
    F->myAssignReset(myPolyValue, myNumTerms);
    F->myMoveToNextLM();
    reduce(F, g);
    F->myRelease(myPolyValue);
    myNumTerms = NumTerms(myPoly()); // LPP,Deg,Comp are the same
  }


  void GPoly::myReduceTail(const Reductors& theReductors)
  {
    CoCoA_ASSERT( !IsZero(*this) );
    ReductionCog F = ChooseReductionCogGeobucket(myGRingInfo());
    F->myAssignReset(myPolyValue, myNumTerms);
    F->myMoveToNextLM();
    reduce(F, mySugar, theReductors);
    F->myRelease(myPolyValue);
    myNumTerms = NumTerms(myPoly()); // LPP,Deg,Comp are the same
  }


  void GPoly::myReduce(const Reductors& theReductors)
  {
    if ( IsZero(*this) ) return;
    ReductionCog F = ChooseReductionCogGeobucket(myGRingInfo());
    //    std::cout << "myreduce1: F is " << F << std::endl;
    F->myAssignReset(myPolyValue, myNumTerms); // myPolyValue gets 0
    //    std::cout << "myreduce2: F is " << F << std::endl;
    if (IsInitialized(mySugar)) //---< AMB 2023-12 >-------
      reduce(F, mySugar, theReductors); // updates mySugar
    else
      reduce(F, theReductors); //----<>---------------------
    //    std::cout << "myreduce3: F is " << F << std::endl;
    F->myRelease(myPolyValue); // myPolyValue gets the value of F
    myUpdateLenLPPLCDegComp(); // myNumTerms, myWDeg, myComponent are updated

    if ( !IsZero(*this) && !IsOne(myLCValue) ) // makes myPolyValue monic
      if ( myGRingInfo().myCoeffRingType()==CoeffEncoding::Field )
        myGRingInfo().myNewSPR()->myDivByCoeff(raw(myPolyValue), raw(myLCValue));
    // if CoeffEncoding::Field myRelease does NOT make poly monic
    // if CoeffEncoding::FrFldOfGCDDomain myRelease makes poly content free
  }




} // end of namespace CoCoA
