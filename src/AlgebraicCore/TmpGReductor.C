//   Copyright (c)  2005-2017  John Abbott, Anna M. Bigatti
//   Authors: 2005-2007  Massimo Caboara, 2016-1017 Anna M. Bigatti

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

#include "CoCoA/TmpGReductor.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/ModuleOrdering.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H" // for dynamic_cast<RingFpImpl*>(CoeffRing.myRingPtr())
#include "CoCoA/RingQQ.H"  // for RingQQ, QQEmbeddingHom
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RealRadical.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/VectorOps.H"  // for printing lists/vectors, HasUniqueOwner
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"
#include "CoCoA/time.H"
#include "CoCoA/verbose.H"  // for VerboseLog
//#include "CoCoA/utils.H" // for LongRange

#include <algorithm>
using std::for_each;
//using std::find_if;  // for IsHomog
#include <functional>
using std::binary_function;
using std::less;
using std::mem_fun_ref; // for calling GPair::complete on GPairList
using std::not1;    // for IsHomog
using std::ptr_fun; // for IsHomog
#include <iostream>
using std::ostream;
using std::endl;
using std::flush;
#include <list>
using std::list;
#include <sstream> // for ostringstream
#include <utility>
using std::make_pair;
#include <vector>
using std::vector;

//static bool ANNA_DEBUG = true;
//static const bool MAX_DEBUG = false;


namespace CoCoA
{
  //  int GReductor::ourDefaultStatLevel = -1;

  // This function produces a new PolyRing using
  // NewPolyRing_DMPII or NewPolyRing_DMPI
  // according to the char.
  // Used internally, indets names are always x[]
//   SparsePolyRing NewSparsePolyRing(const ring& CoeffRing,
//                                    const vector<symbol>& IndetNames,
//                                    const PPOrdering& ord)
//   {
//     if (IsRingFp(CoeffRing)) // special case for SmallFpImpl
//       return NewPolyRing_DMPII(CoeffRing, IndetNames, ord);
//     else
//       return NewPolyRing_DMPI(CoeffRing, IndetNames, ord);
//   } // NewSparsePolyRing


  const GRingInfo& GetGRI(const GPolyList& theGPL)
  {
    CoCoA_ASSERT(!theGPL.empty());
    return theGPL.begin()->myGRingInfo();
  } // GRI


/*  perhaps useless
bool BoolCmpLPPPoly(ConstRefRingElem f, ConstRefRingElem g)
  {
    CoCoA_ASSERT(!IsZero(f));//BUG HUNTING  ???
    CoCoA_ASSERT(!IsZero(g));//BUG HUNTING  ???
    const SparsePolyRing P = owner(f);
    return P->myCmpLPP(raw(f), raw(g))>0;
  }
  */

  bool BoolCmpLPPGPoly(const GPoly& f, const GPoly& g)
  {
    CoCoA_ASSERT(!IsZero(f));//BUG HUNTING  ???
    CoCoA_ASSERT(!IsZero(g));//BUG HUNTING  ???
    //const PPMonoid& PPM1=PPM(owner(f));
    return PPM(owner(f))->myCmp(raw(LPPForOrd(f)),raw(LPPForOrd(g)))>0;
  }

  void GReductor::myCtorAux(const BuchbergerOpTypeFlag theBuchbergerOpType)
  //,                        const UseDynamicAlgFlag IsDynamic)
  {
    myPrepared=false;
    myAgeValue = 0;
    //    myWrongLPPFoundValue=false;

    myBuchbergerOpType=theBuchbergerOpType;

    if (!myCriteria.myCoprime) myStat.myCopLevel=1000;
    if (!myCriteria.myGM) myStat.myGMLevel=1000;
    if (!myCriteria.myBack) myStat.myBCLevel=1000;
  }


  GReductor::GReductor(const GRingInfo& theGRI,
                       const PolyList& TheInputPolys,
                       const BuchbergerOpTypeFlag theBuchbergerOpType,
                       const GBCriteria criteria):
      myGRingInfoValue(theGRI),
      myTrueReductors(theGRI),
      mySPoly(theGRI),
      myOldDeg(GradingDim(theGRI.myNewSPR())),
      myIncomingWDeg(GradingDim(theGRI.myNewSPR())),
      myStat(len(TheInputPolys)),
      myCriteria(criteria)
  {
    myCtorAux(theBuchbergerOpType);

    for (const RingElem& g: TheInputPolys)
    {
      if (!IsZero(g))
        myPolys.push_back(GPoly(g,myGRingInfoValue));
    }
    myPolys.sort(BoolCmpLPPGPoly);
    myTruncDegValue = 0;    
  } //  GReductor ctor


  // This ctor allows to avoid transforming gpolys to polys and back if
  // you want to evaluate complex expressions. The input GPolys may contain
  // zeros and/or not be sorted. Take care.
  GReductor::GReductor(const GRingInfo& theGRI,
                       const GPolyList& TheInputGPolys,
                       const BuchbergerOpTypeFlag theBuchbergerOpType,
                       const GBCriteria criteria):
      myGRingInfoValue(theGRI),
      myTrueReductors(theGRI),
      mySPoly(theGRI),
      myOldDeg(GradingDim(theGRI.myNewSPR())),
      myIncomingWDeg(GradingDim(theGRI.myNewSPR())),
      myStat(len(TheInputGPolys)),
      myCriteria(criteria)
  {
    myPolys=TheInputGPolys;
    myCtorAux(theBuchbergerOpType);
    myTruncDegValue = 0;

    // If you want to remove zeros and/or sort your input GPolys.
    // GPoly Zero(zero(myPolyRing),myGRingInfoValue,Component(LPP(zero(myPolyRing),myGRingInfoValue)),0);
    // myPolys.remove(Zero);
    // myPolys.sort(BoolCmpLPPGPoly);
  } // GReductor ctor


  ostream& operator<<(ostream& out, const GReductor& GR)
  {
    if (!out) return out;  // short-cut for bad ostreams

    out<<"The GROBNER REDUCTOR\n";
    out<<"  Reductors Len="<<GR.myReductorsLen()
        <<"  GB Len="<<GR.myGBasisLen()
        <<"  Pairs Len="<<GR.myPairsLen()
        <<"  Byte Size="<<sizeof(GR)<<"\n";
    GR.myStampaReductors(out);
    GR.myStampaGB(out);
    GR.myStampaPairs(out);
    out<<"\n";
    out<<"The Ring"<<GR.myGRingInfoValue<<"\n";
    out<<"GradingDim  is "<<GradingDim(GR.myGRingInfoValue.myNewSPR())<<endl;
    out<<"The Ring Module Index "<<ModuleVarIndex(GR.myGRingInfoValue.myNewSPR())<<"\n";
    out<<"Age "<< GR.myAgeValue <<"\n";
    out<<"Preparation done? "<<GR.myPrepared<<"\n";
    out<<"myOldDeg "<<GR.myOldDeg<<"\n";
    out<<"myIncomingWDeg " << GR.myIncomingWDeg << "\n";
    //    out<<"Is Dynamic Algorithm? "<<GR.IsDynamicAlgorithm<<"\n";
    //    out<<"Is Wrong LPP been Found? "<<GR.myWrongLPPFoundValue<<"\n";
    out<<"Cop Criterion " <<GR.myCriteria.myCoprime<<"\n";
    out<<"GM Criteria "   <<GR.myCriteria.myGM<<"\n";
    out<<"Back Criterion "<<GR.myCriteria.myBack<<"\n";
    out<<"Div Criterion " <<GR.myCriteria.myDiv<<"\n";
    out<<"Algorithm "<<GR.myBuchbergerOpType<<"\n";
    out<<"\n";
    return out;
  }

  // // Can be used if the Greductor has been updated. If not, you are missing the
  // // Spoly
  // void GReductor::GetCandidateGBasis(PolyList& thePL)const
  // {
  //   PolyList PL;
  //   for (auto& ptr: myGB)
  //     if (IsActive(*ptr))  PL.push_back((*ptr).myPoly());
  //   //if (!IsZero(mySPoly))
  //   //  PL.push_back(GetSPoly().myPoly());
  //   for (const auto& pair: myPairs)
  //     if (pair.IsInputPoly())
  //       PL.push_back(pair.myFirstGPoly().myPoly());
  //   swap(PL,thePL);
  // } // GetCandidateGBasis


  namespace { // anonymous
    void VERBOSE_NewPolyInGB(VerboseLog& VERB, long LenGB, long LenGPair, const GPoly& SPoly)
    {
      if (VerbosityLevel()>=100)
      {
        VERB(100) << "--New poly in GB:"
                  << " len(GB) = " << LenGB
                  << " len(pairs) = " << LenGPair << endl;
        VERB(101) << "--NumTerms = " << NumTerms(SPoly)
                  << " wdeg = " << wdeg(SPoly) << endl;
        VERB(150) << "--New poly is " << poly(SPoly) << endl;
      }
    }
  } // namespace anonymous


  void GReductor::myStampaGB(ostream& out)const
  {
    if (!out) return;  // short-cut for bad ostreams

    out<<"The GBASIS\n";
    for (const auto& ptr: myGB) out<<*ptr<<","<<endl;
    out<<endl;
  } // myStampaGB



  void GReductor::myStampaPPGB(ostream& out)const
  {
    if (!out) return;  // short-cut for bad ostreams

    out<<"The GBASIS\n";
    for (const auto& ptr: myGB)
      out << LPPForDiv(*ptr) << "," << endl;
    out << endl;
  } // myStampaPPGB


  void GReductor::myStampaPairs(ostream& out)const
  {
    if (!out) return;  // short-cut for bad ostreams
    out<<"The PAIRS\n"<<myPairs<<"\n";
  } // myStampaPairs


  void GReductor::myStampaReductors(ostream& out)const
  {
    if (!out) return;  // short-cut for bad ostreams

    myTrueReductors.myStampaReductors(out);
    out<<endl;
  } // myStampaReductors

  
// This procedure may be substituted by a transform_if
  void GReductor::myCopyGBasis(PolyList& GBasis)
  {
    GBasis.clear();
    for (const auto& ptr: myGB)
      if (IsActive(*ptr))
        GBasis.push_back((*ptr).myPoly());
  }


  void GReductor::myCopyMinGens(PolyList& MinGens)
  {
    MinGens.clear();
    for (const auto& ptr: myGB)
      if (IsMinimalGen(*ptr))
        MinGens.push_back((*ptr).myPoly());
  }


 // This procedure may be substituted by a transform_if
  void GReductor::myCopyGBasis(GPolyList& GBasis)
  {
    GBasis.clear();
    for (const auto& ptr: myGB)
      if (IsActive(*ptr))
        GBasis.push_back((*ptr));
  }
  

// This procedure may be substituted by a transform_if
  void GReductor::myGBasisClear(GPolyList& theGBasis)
  {
    theGBasis.clear();
    for (auto& ptr: myGB)
      if (IsActive(*ptr))
      {
        GPoly Zero(myGRingInfoValue);
        theGBasis.push_back(Zero);
        theGBasis.back().AssignClear(*ptr);
      }
    myGB.clear();
  }


  void GReductor::myCopyGBasis(VectorList& GBasis)
  {
    GBasis.clear();
    if (myGB.empty()) return;
    for (const auto& ptr: myGB)
      if (IsActive(*ptr))
        GBasis.push_back(DeEmbedPoly((*ptr).myPoly(),myGRingInfoValue));
  }


  void GReductor::myCopyMinGens(VectorList& MinGens)
  {
    MinGens.clear();
    for (const auto& ptr: myGB)
      if (IsMinimalGen(*ptr))
        MinGens.push_back(DeEmbedPoly((*ptr).myPoly(),myGRingInfoValue));
  }



//esame piu' approfondito - sia correttezza sia efficienza
//warning: il CopCriterion is false, still the coprimality is used to decide
//         which pair kill
  long GReductor::myGMInsert(GPairList& L,GPair P)
  {
    long NumPairsConsidered=0;
    bool ToBeInserted=true;
    bool erased=false;
    long P_Component=GPairComponent(P);

    GPairList::iterator it=L.begin();
    while (it!=L.end() && ToBeInserted)
    {
      ++NumPairsConsidered;
      if (P_Component==GPairComponent(*it)
          && IsDivisibleFast(LCMwMask(P), LCMwMask(*it)))
      {
        ToBeInserted=false;
        if (LCMwMask(P)==LCMwMask(*it) && (!it->IamCoprime()) && (P.IamCoprime()))
        {
          //          P.myComplete();          it->myComplete();
          *it=P;
        }
      }
      else
        if (P_Component==GPairComponent(*it)
            && IsDivisibleFast(LCMwMask(*it), LCMwMask(P)) )
        {
          it=L.erase(it);
          erased=true;
        }
      if (!erased)  ++it;
      erased=false;
    } // while
    if (ToBeInserted)
      L.push_back(P);// new pairs are sorted later
    return NumPairsConsidered;
  }


  void GReductor::myBuildNewPairs(GPairList& new_pairs)
  {
    VerboseLog VERBOSE("myBuildNewPairs");
    long standing_index = len(myGB)-1;
    long walking_index = 0;
    long inserted_pairs = 0;//STAT

    GPolyPtrList::const_iterator last=myGB.end(); --last;
    long last_component = component(**last);
    GPolyPtrList::const_iterator it;
    for (it=myGB.begin(); it!=last; ++it,++walking_index)
    {
      if (IsActive(**it)&&last_component == component(**it))
      {
        //std::clog<<"walking_component "<<Component(**it)<<endl;
        if (myCriteria.myDiv && IsDivisibleFast(LPPForDivwMask(**it), LPPForDivwMask(**last)) ) // speed is not necessary
        {
          if (VerbosityLevel() >= myStat.myPolyDeletedLevel)
            VERBOSE(myStat.myPolyDeletedLevel) <<"<"<<walking_index
                                               <<"> "<<LPPForDiv(**it)
                                               <<" DELETED BY NEW "
                <<"<"<<standing_index<<"> "<<LPPForDiv(**last)<< endl;
          (*it)->Deactivate();
          ++myStat.myPolyDeleted;
        } //second if
        //else
        ++inserted_pairs;
        ++myStat.myPInserted;
        if (myCriteria.myGM)
          myStat.myGMTouched+=myGMInsert(new_pairs, GPair(**it,**last));
        else
        {
          myStat.myGMTouched=0;
          Ordered_Insert(new_pairs,GPair(**it,**last));
        }
      };//Active if
    };//for

    myStat.myGMKilled+=inserted_pairs-len(new_pairs);
    if (VerbosityLevel() >= myStat.myGMLevel && (inserted_pairs!=len(new_pairs)))
      VERBOSE(myStat.myGMLevel) << "[GM KILLED "
                                << inserted_pairs-len(new_pairs)
                                << " OUT OF " << inserted_pairs << "]" << endl;

    long pre_cop_test=len(new_pairs);
    if (myCriteria.myCoprime)
    {
      GPairList::iterator it1=new_pairs.begin();
      while (it1!=new_pairs.end())
        if (it1->IamCoprime())
          it1=new_pairs.erase(it1);
        else
          ++it1;
    }
    if (pre_cop_test!=len(new_pairs))
      if (VerbosityLevel() >= myStat.myCopLevel)
        VERBOSE(myStat.myCopLevel) <<"[COP KILLED "<<pre_cop_test-len(new_pairs)
                                   <<" OUT OF "<<inserted_pairs<<"]" << endl;
    myStat.myCopKilled+=pre_cop_test-len(new_pairs);
    if (VerbosityLevel() >=  myStat.myNewPairsLevel)
    {
      std::ostringstream oss;
      for (const auto& pair: new_pairs)
        oss << "<"<<pair.myFirstIndex() << "," << pair.mySecondIndex() << ">, ";
      VERBOSE(myStat.myNewPairsLevel) << len(new_pairs)
                                      << " new pairs: " << oss.str() << endl;
    };
  } // myBuildNewPairs


  // void GReductor::myBuildNewPairsAll(GPairList& new_pairs)
  // {
  //   VerboseLog VERBOSE("myBuildNewPairsAll");
  //   new_pairs.clear();
  //   long standing_index=1;
  //   long walking_index=0;
  //   long inserted_pairs=0;//STAT

  //   GPolyPtrList::const_iterator BeginPlus=myGB.begin(); ++BeginPlus;
  //   for (GPolyPtrList::const_iterator last=BeginPlus;last!=myGB.end();++last,++standing_index)
  //   {
  //   long last_component=Component(**last);
  //   GPolyPtrList::const_iterator it;
  //   walking_index=0;
  //   for (it=myGB.begin(); it!=last; it++,walking_index++)
  //   {
  //     if (IsActive(**it)&&last_component==Component(**it))
  //     {
  //       if (myCriteria.myDiv && IsDivisibleFast(LPPForDivwMask(**it), LPPForDivwMask(**last)) )
  //       {
  //         if (VerbosityLevel() >= myStat.myPolyDeletedLevel)
  //           VERBOSE(1)
  //             <<"<"<<walking_index<<"> "<<LPPForDiv(**it)<<" DELETED BY NEW "
  //             <<"<"<<standing_index<<"> "<<LPPForDiv(**last)<<"\n";
  //         (*it)->Deactivate();
  //         myStat.myPolyDeleted++;
  //       } //second if
  //       //else
  //       inserted_pairs++;
  //       myStat.myPInserted++;
  //       if (myCriteria.myGM)
  //         myStat.myGMTouched+=myGMInsert(new_pairs, GPair(**it,**last));
  //       else
  //       {
  //         myStat.myGMTouched=0;
  //         Ordered_Insert(new_pairs, GPair(**it,**last));
  //       }
  //     };//Active if
  //   };//walking for
  //   } //  last for
  //   myStat.myGMKilled+=inserted_pairs-len(new_pairs);
  //   if (VerbosityLevel() >= myStat.myGMLevel && (inserted_pairs!=len(new_pairs)))
  //     VERBOSE(1) <<"[GM KILLED "<<inserted_pairs-len(new_pairs)
  //                <<" OUT OF "<<inserted_pairs<<"]\n";

  //   long pre_cop_test=len(new_pairs);
  //   if (myCriteria.myCoprime)
  //   {
  //     GPairList::iterator it1=new_pairs.begin();
  //     while (it1!=new_pairs.end())
  //       if (it1->IamCoprime())
  //         it1=new_pairs.erase(it1);
  //       else
  //         ++it1;
  //   }
  //   if (pre_cop_test!=len(new_pairs))
  //     if (VerbosityLevel() >= myStat.myCopLevel)
  //       VERBOSE(1) <<"[COP KILLED "<<pre_cop_test-len(new_pairs)
  //                  <<" OUT OF "<<inserted_pairs<<"]\n";
  //   myStat.myCopKilled+=pre_cop_test-len(new_pairs);
  //   if (VerbosityLevel() >= myStat.myNewPairsLevel)
  //   {
  //     VERBOSE(1) << "NEW PAIRS " << len(new_pairs) << "  ";
  //     for (const auto& p: new_pairs)
  //       VERBOSE(1) << "<" << p.myFirstIndex()
  //                  << "," << p.mySecondIndex() << "> " << std::endl;
  //   };
  // } // myBuildNewPairsAll


  void GReductor::myApplyBCriterion()
  {
    VerboseLog VERBOSE("myApplyBCriterion");
    const PPWithMask& newPP(LPPForDivwMask(myPolys.back()));
    const long newPP_component = component(myPolys.back());
    GPairList::iterator it=myPairs.begin();
    while (it!=myPairs.end())
    {
      ++myStat.myBTouched;
      if (GPairComponent(*it)!=newPP_component
          ||
          it->BCriterion_OK(newPP))
      {++it;}
      else
      {
        if (VerbosityLevel() >= myStat.myPolyDeletedLevel)
          VERBOSE(1)<<*it<<" KILLED BY NEW POLY "<<newPP<<" BCRIT\n";
        ++myStat.myBKilled;
        it=myPairs.erase(it);
      };// else
    };//While
  }


  void GReductor::myUpdateBasisAndPairs()
  { myUpdateBasisAndPairs(mySPoly); }


  void GReductor::myUpdateBasisAndPairs(const GPoly& inPoly)
  {
    VerboseLog VERBOSE("myUpdateBasisAndPairs");
    // --> if non-homog, inPoly must have sugar <--
    //if (!IsTrueGCDDomain(CoeffRing(mySPoly)))
    //myTrueReductors.OrderedInterreduce(mySPoly);  // ANNA
    VERBOSE(100) << wdeg(inPoly) << std::endl;
    if (myGRingInfo().myInputAndGrading()==HOMOG)
      myTrueReductors.myReduceTails(inPoly);// this must be the first op
    ++myAgeValue;
    if (IsConstant(poly(inPoly)))
    {
      VERBOSE(100) << "!!!!!!!!!! ideal(1) !!!!!!!!!!!!" << std::endl;
      myPolys.clear();
      myPolys.push_back(inPoly);
      myGB.clear();
      myGB.push_back(&myPolys.back());
      myTrueReductors.myClear();
      myTrueReductors.Insert(&myPolys.back());
      myPairs.clear();
      return;
    }
    myPolys.push_back(inPoly);// this must be the second op
    //myTrueReductors.SuperInterreduce(mySPoly);  // ANNA
    myTrueReductors.Insert(&myPolys.back());
    myGB.push_back(&myPolys.back());
    GPairList new_pairs;
    myBuildNewPairs(new_pairs);
    if (myCriteria.myBack)
      myApplyBCriterion();
    for_each(new_pairs.begin(), new_pairs.end(), [](GPair& arg) { arg.myComplete(); });
// ORIGINAL LINE   for_each(new_pairs.begin(), new_pairs.end(), mem_fun_ref(&GPair::myComplete));
    new_pairs.sort();
    myPairs.merge(new_pairs);
  }


  void GReductor::myReduceCurrentSPoly()
  {
    const char* const FnName = "myReduceCurrentSPoly";
    VerboseLog VERBOSE(FnName);
    //    CheckForInterrupt(FnName);
    //    myGRingInfo().myCheckForTimeout(FnName);
    if (!myPrepared)
    {
      VERBOSE(10)<<"GReductor: no preparation done ";
      return;
    }
    //if (myStat.myDegLevel)
    //  clog<<"GReductor:myReduceCurrentSPoly performing one pair handling";
    if (myPolys.empty() || myPairs.empty())  return;
//     if (VerbosityLevel() >= myStat.myNumPairLevel)
//       VERBOSE(myStat.myNumPairLevel) << "len(pairs) = " << len(myPairs) << endl;
    VERBOSE(myStat.myReductionLevel+1) << myPairs.front() << endl;
    if (myPairs.front().IamCoprime() && myCriteria.myCoprime)
    {
      mySPoly = GPoly(myGRingInfoValue); // initialized as 0
      VERBOSE(myStat.myCopLevel) << "coprime" << endl;
      ++myStat.myCopKilled;
      myStat.myReductionTime=0.0;
    }
    else
    {
      mySPoly.myAssignSPoly(myPairs.front(), myAgeValue);  // ??? SPoly computed only if not coprime
      if (myPairs.front().IsInputPoly())  mySPoly.mySetMinimalGen();
      VERBOSE(200) << " --before: " << poly(mySPoly) << std::endl;
      myStat.myPolyLens.push_back(make_pair(NumTerms(mySPoly),0));
      std::ostringstream VERB_s;
      if (VerbosityLevel() >= myStat.myReductionLevel)
      {
        if (IsZero(mySPoly)) VERB_s << " = 0";
        else VERB_s << "len(SPoly)=" << myStat.myPolyLens.back().first;
      }
      myStat.myReductionTime = CpuTime();
      mySPoly.myReduce(myTrueReductors); // interrupt/timeout
      myStat.myReductionTime -= CpuTime();
      VERBOSE(200) << " --after:  " << poly(mySPoly) << std::endl;
      ++myStat.myNReductions;
      myStat.myPolyLens.back().second = NumTerms(mySPoly);
      if (IsZero(mySPoly)) ++myStat.myUseless; else ++myStat.myUseful;
      // if (!IsTrueGCDDomain(CoeffRing(mySPoly)) || NumTerms(mySPoly)<=2)
      // myTrueReductors.interreduce(mySPoly);
      if (VerbosityLevel() >= myStat.myReductionLevel)
      {
        VERB_s << " --> ";
        if (IsZero(mySPoly)) VERB_s << "0";
        else
        {
          //          VERB_s << "<" << len(myGB) << ">: ";
          VERB_s << LPPForDiv(mySPoly) << "+(" << myStat.myPolyLens.back().second << ")";
          if (component(mySPoly)!=0) VERB_s <<" Comp =" << component(mySPoly) << endl;
          if (NumTerms(mySPoly)<=2) VERB_s <<" short-reducer deg =" << wdeg(mySPoly) << endl;
        }
        VERB_s << " time=" << std::floor(-myStat.myReductionTime);
        VERBOSE(myStat.myReductionLevel) << VERB_s.str() << endl;
      } // (VerbosityLevel() >= myStat.myReductionLevel)
    } // else -if not coprime
    myPairs.erase(myPairs.begin()); // erase the used gpair
    //if (!IsZero(mySPoly)) myUpdateBasisAndPairs();

//// move this check into a more appropriate position (after adding the new pairs)
    // if (!myPairs.empty())
    // {
    //   myFirstPairWDeg = wdeg(myPairs.front()); // ANNA: next pair to do
    //   if (myFirstPairWDeg != myOldDeg)
    //   {
    //     VERBOSE(myStat.myDegLevel) << "--NEW DEG: myFirstPairWDeg = "
    //                                << myFirstPairWDeg
    //                                << "\tSugar = " << sugar(myPairs.front())
    //                                << endl;
    //     myStat.myUpgradeDegStats(myOldDeg, len(myPairs));
    //     myOldDeg = myFirstPairWDeg;
    //   } //if (myFirstPairWDeg!=myOldDeg)
    // } //if (!myPairs.empty())
    // else myStat.myUpgradeDegStats(myFirstPairWDeg,0);
  } // myReduceCurrentSPoly


  void GReductor::myUpdateIncomingWDeg()
  {
    VerboseLog VERBOSE("myUpdateIncomingWDeg");
    if (!myPairs.empty())
    {
      myIncomingWDeg = wdeg(myPairs.front()); // ANNA: next pair to do
      if (myIncomingWDeg != myOldDeg)
      {
        VERBOSE(myStat.myDegLevel) << "--NEW DEG: myIncomingWDeg = "
                                   << myIncomingWDeg
                                   << "\tSugar = " << sugar(myPairs.front())
                                   << endl;
        myStat.myUpgradeDegStats(myOldDeg, len(myPairs));
        myOldDeg = myIncomingWDeg;
      } //if (myFirstPairWDeg!=myOldDeg)
    } //if (!myPairs.empty())
    else myStat.myUpgradeDegStats(myIncomingWDeg, 0);
  }
  
  
  // ANNA: called by TmpGOperations
  void GReductor::myDoGBasis()
  {
    //    CoCoA_ASSERT(myGetBuchbergerOpType() == HomogeneousAlg);
    VerboseLog VERBOSE("myDoGBasis");
    myPrepareGBasis();
    while (!myPairs.empty())
    {
      myReduceCurrentSPoly();
      if (!IsZero(mySPoly))
      {
        myUpdateBasisAndPairs();
        VERBOSE_NewPolyInGB(VERBOSE, len(myGB), len(myPairs), mySPoly);
        // if (!myPairs.empty())
        //   if (myIncomingWDeg!=myOldDeg && myTrueReductors.IhaveBorelReductors())
        //     myTrueReductors.myBorelReductorsUpdateInNextDegree();
      }
      myUpdateIncomingWDeg();
      if (myTruncDeg() != 0)
      {
        if (myPairs.empty())  mySetTruncDeg(0); // GB is complete
        else
          if (myTruncDeg() < myIncomingWDeg[0])
          {
            VERBOSE(100) << "myTruncDeg() =" << myTruncDeg()
                         << " < myIncomingWDeg[0]: "
                         << myIncomingWDeg << std::endl;
            break;
          }
      }
    } // while
    VERBOSE(100) << "--Final clean up ... " << endl;
    myFinalizeGBasis();
    if (VerbosityLevel() >= myStat.myFinalLevel)
      myStat.myStampa(VERBOSE(myStat.myFinalLevel));
  } // myDoGBasis


// ANNA: called by TmpGOperations (ComputeElimFirst)
  RingElem GReductor::myDoGBasisElimFirst(ConstRefPPMonoidElem ElimIndsProd)
  {
    VerboseLog VERBOSE("myDoGBasisElimFirst");
    //    CoCoA_ASSERT(myGetBuchbergerOpType() == HomogeneousAlg);
    myPrepareGBasis();
    //clog << "\nord = " << PPM(owner(poly(mySPoly))) << endl;
    //clog << "\ngens = " << myPairs << endl;
    while (!myPairs.empty())
    {
      myReduceCurrentSPoly();
      if (!IsZero(mySPoly))
      {
        if (IsCoprime(LPP(poly(mySPoly)), ElimIndsProd))
        {
          VERBOSE(100) << "--First Elim Poly found:"
                     << " len(GB) = " << len(myGB)
                     << " len(pairs) = " << len(myPairs)
                     << endl;
          myStat.myTotalTime-=CpuTime();
          if (VerbosityLevel() >= myStat.myFinalLevel)
            myStat.myStampa(VERBOSE(myStat.myFinalLevel));
          return poly(mySPoly);
        }
        myUpdateBasisAndPairs();
        VERBOSE_NewPolyInGB(VERBOSE, len(myGB), len(myPairs), mySPoly);
        // if (!myPairs.empty())
        //   if (myIncomingWDeg!=myOldDeg&&myTrueReductors.IhaveBorelReductors())
        //     myTrueReductors.myBorelReductorsUpdateInNextDegree();
      }
    } // while
    return zero(RingZZ()); // just to keep the compiler quiet
  } // myDoGBasisElimFirst


  // ANNA: called by TmpGOperations
  void GReductor::myDoGBasisSelfSatCore()
  {
    VerboseLog VERBOSE("myDoGBasisSelfSatCore");
    CoCoA_ASSERT(myGetBuchbergerOpType() == SaturatingAlg);
    VERBOSE(100) << "ring is " << myGRingInfoValue.myNewSPR() << std::endl;
    VERBOSE(100) << ordering(PPM(myGRingInfoValue.myNewSPR())) << endl;
    const long HIndetIndex=NumIndets(myGRingInfoValue.myNewSPR())-1;// This is OK for Ideals only
    degree SPolyPredDeg(GradingDim(myGRingInfoValue.myNewSPR()));// Used for stats in the dehmog alg
    myPrepareGBasis();
    //    double T=0.0;
    while (!myPairs.empty())
    {
      myReduceCurrentSPoly();
      if (!IsZero(mySPoly))
      {
      	SPolyPredDeg = wdeg(mySPoly);
        //mySPoly.smart_dehomog_DRL(HIndetIndex);
        mySPoly.smart_dehomog(HIndetIndex);
      	if (SPolyPredDeg!=wdeg(mySPoly))
      	{
          ++myStat.myPolyDHed;
          myStat.myDegDH += ConvertTo<long>((SPolyPredDeg-wdeg(mySPoly))[0]);
          //          if (VerbosityLevel() >= myStat.myPolyDHLevel && false)
          VERBOSE(100) << "Lower degree: "
                 <<SPolyPredDeg<< "-->"<<wdeg(mySPoly)<<" "
                 <<LPPForDiv(mySPoly)<<endl;
        } // if
        myUpdateBasisAndPairs();
        VERBOSE_NewPolyInGB(VERBOSE, len(myGB), len(myPairs), mySPoly);
        // if (!myPairs.empty())
        //   if (myIncomingWDeg!=myOldDeg&&myTrueReductors.IhaveBorelReductors())
        //     myTrueReductors.myBorelReductorsUpdateInNextDegree();
      }
    } // while
    VERBOSE(100) << "--Final clean up ... " << endl;
    myFinalizeGBasis();
  } // myDoGBasisSelfSatCore


  void GReductor::myDoGBasisRealSolve()
  {
    VerboseLog VERBOSE("myDoGBasisRealSolve");
    SparsePolyRing P = myGRingInfo().myNewSPR();
    if (!IsZZ(CoeffRing(P))) CoCoA_THROW_ERROR("must be over QQ","myDoGBasisRealSolve");
    SparsePolyRing PQQ = NewPolyRing(RingQQ(), NewSymbols(NumIndets(P)));
    RingHom P_PQQ = PolyRingHom(P, PQQ, ZZEmbeddingHom(PQQ), indets(PQQ));
    RingHom PQQ_P = PolyRingHom(PQQ, P, QQEmbeddingHom(P), indets(P));

    myPrepareGBasis();
    while (!myPairs.empty())
    {
      myReduceCurrentSPoly();
      if (!IsZero(mySPoly))
      {
//        RingElem RadSPoly = radical(P_PQQ(mySPoly.myPoly()));
        RingElem RadSPoly = RealRadical(P_PQQ(mySPoly.myPoly()));
        if (deg(RadSPoly) < deg(mySPoly.myPoly()))
          VERBOSE(70) << mySPoly.myPoly() << " [RealRadical]--> "
                       << RadSPoly << std::endl;
        GPoly GP(PQQ_P(ClearDenom(RadSPoly)), myGRingInfoValue, clear);
        GP.myInitializeSugar(sugar(mySPoly));
        myUpdateBasisAndPairs(GP);
        VERBOSE_NewPolyInGB(VERBOSE, len(myGB), len(myPairs), GP);
      }
    } // while
    VERBOSE(100) << "--Final clean up ... " << endl;
    myFinalizeGBasis();
    if (VerbosityLevel() >= myStat.myFinalLevel)
      myStat.myStampa(VERBOSE(myStat.myFinalLevel));
  } // myDoGBasisRealSolve


//   void GReductor::myDoGBasisByBatch()
//   {
//     VerboseLog VERBOSE("myDoGBasisByBatch");
//     bool NonZeroReductionPerformed=false;
//     myPrepareGBasis();

//     // process the special pairs
//     while (!myPairs.empty())
//     {
//       myReduceCurrentSPoly();
//       if (!IsZero(mySPoly))
//       {
//         myUpdateBasisOnly();
//         VERBOSE_NewPolyInGB(VERBOSE, len(myGB), len(myPairs), mySPoly);
//       }
//     } // while

//     do
//     {
//       NonZeroReductionPerformed=false;
//       myCreatePairs();
//       while (!myPairs.empty())
//       {
//         myReduceCurrentSPoly();
//         if (!IsZero(mySPoly))
//         {
//           myUpdateBasisOnly();
//           VERBOSE_NewPolyInGB(VERBOSE, len(myGB), len(myPairs), mySPoly);
//           NonZeroReductionPerformed=true;
//         }
//       } // while
//     }
//     while (NonZeroReductionPerformed);
//     myFinalizeGBasis();
//   } // myDoGBasisByBatch


//   void GReductor::myReduceUntilNonZeroRedSPoly()
//   {
//     VerboseLog VERBOSE("myReduceUntilNonZeroRedSPoly");
//     bool flag=true;
//     while ((!myPairs.empty())&&(IsZero(mySPoly)||flag))
//     {
//       flag=false;
//       myReduceCurrentSPoly();
//       if (!IsZero(mySPoly))
//       {
//         myUpdateBasisAndPairs();
//         VERBOSE_NewPolyInGB(VERBOSE, len(myGB), len(myPairs), mySPoly);
//         if (!myPairs.empty())
//           if (myFirstPairWDeg!=myOldDeg&&myTrueReductors.IhaveBorelReductors())
//             myTrueReductors.myBorelReductorsUpdateInNextDegree();
//       }
//     }
//   } // myReduceUntilNonZeroRedSPoly


//   void GReductor::myReduceUntilNonZeroRedSPoly()
//   {
//     //    VerboseLog VERBOSE("myReduceUntilNonZeroRedSPoly");
//     while (!myPairs.empty())
//     {
//       myReduceCurrentSPoly();
//       if (!IsZero(mySPoly)) return;
//     }
//   } // myReduceUntilNonZeroRedSPoly


//   void GReductor::myCreatePairs()
//   {
// clog << "myCreatePairs:begin "<<endl;
//       GPairList new_pairs;
//       myBuildNewPairsAll(new_pairs);
//       if (myCriteria.myBack)
//         myApplyBCriterion();
//       for_each(new_pairs.begin(), new_pairs.end(), mem_fun_ref(&GPair::myComplete));
//       new_pairs.sort();
//       swap(myPairs,new_pairs);
// clog << "myCreatePairs:end "<<endl;
//   } // myCreatePairs


  void GReductor::myCreateInputPolyPairs(GPolyList& theCandidateBasis)
  {
    for (const auto& p: theCandidateBasis)
      Ordered_Insert(myPairs, GPair(p));
  } // myCreateInputPolyPairs


  // Prepare the first pairs; "input poly pairs"
  void GReductor::myPrepareGBasis()
  {
    VerboseLog VERBOSE("myPrepareGBasis");
    for (auto& p: myPolys)
      p.myInitializeSugar(myGRingInfoValue.myNewSugar(poly(p)));
    myCreateInputPolyPairs(myPolys);
    myPrepared=true;
    if (myGRingInfoValue.IamModule())
      myCriteria.myCoprime = false;// CopCriterion works only for REAL ideals
    CoCoA_ASSERT(len(myPairs)!=0);
    myIncomingWDeg = wdeg(myPairs.front());
    myOldDeg = myIncomingWDeg;
    myStat.myReductionTime=0.0; // STAT
    myStat.myTotalTime=CpuTime(); // STAT
    VERBOSE(myStat.myDegLevel) <<"--myIncomingWDeg="<<myIncomingWDeg<<endl;
  } // myPrepareGBasis


  void GReductor::myFinalizeGBasis()
  {
    const char* const FnName = "myFinalizeGBasis";
    VerboseLog VERBOSE(FnName);
    myStat.myTotalTime-=CpuTime();
    // interreduction
    if (true)  // always reduced gbasis
      if (myGRingInfo().myInputAndGrading()!=HOMOG) //myUpdateBasisAndPairs
      {
        VERBOSE(105) << "interreducing..." << std::endl;
        for (auto it=myGB.begin(); it!=myGB.end(); /*++it*/)
        {
          CheckForInterrupt(FnName);
          myGRingInfo().myCheckForTimeout(FnName);
          vector<ReductorData>::iterator it1=myTrueReductors.find(*it);
          it1->mySetIamNotToBeUsed(true);
          if (FindReducer(**it, myTrueReductors) != nullptr)  // surely -->0
          //if (IsZero(**it))    //  remove it
          {
            VERBOSE(110) << "--> zero" << endl;
            it = myGB.erase(it);
            //if (it != myGB.begin()) --it;
          }
          else
          {
            (**it).myReduce(myTrueReductors); // reduces tail
            CoCoA_ASSERT_ALWAYS(!IsZero(**it));
            VERBOSE(110) << "--> non zero" << endl;
            ++it;
            it1->mySetIamNotToBeUsed(false);
          }
        } // myGB for
      } // interreduction
  } // myFinalizeGBasis


//   // ANNA: not called
//   void GReductor::myDoGBasisTEST()
//   {
//     double T=0.0;
// //clog<<"myDoGBasisTEST:begin"<<endl;
//     // Input Polynomials sorted and zero polys deleted
//     if (myPolys.empty()) return;
//     if (myGRingInfoValue.IamModule())
//       myCriteria.myCoprime = false;// CopCriterion works only for REAL ideals
//     GPoly Zero(myGRingInfoValue);
//     GPolyList::iterator it=myPolys.begin();
//     for (long i=0;it!=myPolys.end();++it)
//       Ordered_Insert(myPairs,GPair(*it));
//     myOldDeg = wdeg(myPairs.front());
//     myFirstPairWDeg = myOldDeg;
//     if (myStat.myDegLevel)
//       clog<<"\n[log] ********* Starting_Pair_Degree="<<myFirstPairWDeg<<"\n";
//     while (!myPairs.empty())
//     {
//       if (myStat.myNumPairLevel)
//         clog<<"[log] pair="<<len(myPairs)<<" "<<flush;
//       if (myStat.myReductionLevel)
//         clog<<"doing="<<myPairs.front();
//       if (myPairs.front().IamCoprime()&& myCriteria.myCoprime)
//       {
//         mySPoly=Zero;
//         if (myStat.myCopLevel)
//           clog<<"COP"<<flush;
//         ++myStat.myCopKilled;
//         T=0.0;
//        }
//        else
//        {
//          mySPoly.myAssignSPoly(myPairs.front(),myAgeValue);  // ??? SPoly computed only if not coprime
//          ++myAgeValue;
//          myStat.myPolyLens.push_back(make_pair(NumTerms(mySPoly),0));
//          if (myStat.myReductionLevel && IsZero(mySPoly))
//             clog<<"SPOLY IS ZERO BEFORE REDUCTION !"<<endl;
//          if (myStat.myReductionLevel)
//             clog<<" SPolyLen="<<myStat.myPolyLens.back().first<<" reduction="<<flush;
//          myPairs.erase(myPairs.begin());// erase the used gpair
//          T=CpuTime();
//          mySPoly.myReduce(myTrueReductors);
//          ++myStat.myNReductions;
//          T-=CpuTime();
//          myStat.myPolyLens.back().second=NumTerms(mySPoly);
//          if (IsZero(mySPoly))
//          {
//            if (myStat.myReductionLevel) clog << "0"<<flush;
//            ++myStat.myUseless;
//          }
//          else
//          {
//            if (myStat.myReductionLevel)
//              clog<<LPPForDiv(mySPoly)<<"+..<"<<len(myGB)+1<<"> Len="
//                  <<myStat.myPolyLens.back().second
//                  <<" Comp="<<Component(mySPoly)<<flush;
//            ++myStat.myUseful;
//            //if (!IsTrueGCDDomain(CoeffRing(mySPoly)) || NumTerms(mySPoly)<=2)
//            //myTrueReductors.interreduce(mySPoly);
//            if ((myStat.myReductionLevel && (NumTerms(mySPoly)<=2)))
//              clog<<" EASY REDUCTOR FOUND  LEN="<<NumTerms(mySPoly)<<" DEG="<<wdeg(mySPoly)<<endl;
//          } //  else
//        } // else -if not coprime
//        if (myStat.myReductionLevel)
//           clog <<" time="<<-T<<" \n"<<flush;
//        if (!IsZero(mySPoly)) myUpdateBasisAndPairs();
//        if (!myPairs.empty())
//        {
//          myFirstPairWDeg=wdeg(myPairs.front());//STAT
//          if (myFirstPairWDeg!=myOldDeg)
//          {
//            if (myStat.myDegLevel)
//              clog<<"\n[log] ********* Current_Pair_Degree_Now="<<myFirstPairWDeg<<"\n";
//            if (myTrueReductors.IhaveBorelReductors())
//              myTrueReductors.myBorelReductorsUpdateInNextDegree();
//            myStat.myUpgradeDegStats(myOldDeg,len(myPairs));
//            myOldDeg=myFirstPairWDeg;
//          }
//        }
//        else
//          myStat.myUpgradeDegStats(myFirstPairWDeg,0);
//     };//end main while
//     myStat.myTotalTime-=CpuTime();
//     myStat.myStampa(clog);
//   } //  end _myDoGBasis



//At the moment, no one is calling this procedure. The caller calls myDoGBasisTEST instead
//   void GReductor::myDoSATMixGBasis()
//   {
//     CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "myDoSATMixGBasis: should not be called!");
//     // MOD the homog var and the module vars conflict
//     // Input Polynomials sorted and zero polys deleted
//     if (myPolys.empty()) return;

//     if (myGRingInfoValue.IamModule())
//       myCriteria.myCoprime = false;// CopCriterion works only for REAL ideals

//     CoCoA_ASSERT(GradingDim(myGRingInfoValue.myNewSPR()) == 1);
//     const long HIndetIndex=NumIndets(myGRingInfoValue.myNewSPR())-1;// This is OK for Ideals only
// //VARIABLE SET NEVER USED    bool satjumped=false;
//     double T=0.0;//STAT
//     degree SPolyPredDeg(GradingDim(myGRingInfoValue.myNewSPR()));// Used for stats in the dehmog alg

//     GPoly Zero(myGRingInfoValue);


//     GPolyList::iterator it=myPolys.begin();
//     for (long i=0;it!=myPolys.end();++it)
//     {
//       Ordered_Insert(myPairs,GPair(*it));
//     }

//     myOldDeg = wdeg(myPairs.front());//STAT
//     myFirstPairWDeg = myOldDeg;//STAT

//     if (VerbosityLevel() >= myStat.myDegLevel)
//       clog<<"\n[log] ********* Starting_Pair_Degree="<<myOldDeg<<"\n";


//     while (!myPairs.empty())
//     {
// //VARIABLE SET NEVER USED      satjumped=false;
//       myFirstPairWDeg = wdeg(myPairs.front());//STAT

//       if (myFirstPairWDeg!=myOldDeg)
//       {
//         if (VerbosityLevel() >= myStat.myDegLevel)
//           clog<<"\n[log] New_Degree="<<myFirstPairWDeg<<"\n";
//         //      myTrueReductors.AdjustLimit(myFirstPairWDeg);
//         myStat.myUpgradeDegStats(myOldDeg,len(myPairs));
//         myOldDeg=myFirstPairWDeg;
//       }

//       if (!myPairs.empty())
//       {
//         // WARN - POSSIBLY DIFFERENT
//         if (VerbosityLevel() >= myStat.myNumPairLevel)
//           clog<<"[log] pair="<<len(myPairs)<<" "<<flush;
//         if (VerbosityLevel() >= myStat.myReductionLevel)
//           clog<<"doing="<<myPairs.front()<<" Len="
//               <<myStat.myPolyLens.back().first<<"-->"<<flush;
//         if (myPairs.front().IamCoprime()&& myCriteria.myCoprime)
//         {
//           mySPoly=Zero;
//           if (VerbosityLevel() >= myStat.myCopLevel) clog<<"COP"<<flush;
//           ++myStat.myCopKilled;
//           //--myStat.myUseless;
//           T=0.0;
//         }
//         else
//         {
//           mySPoly.myAssignSPoly(myPairs.front(), ++myAgeValue);  // ??? SPoly computed only if not coprime
//           ++myAgeValue;
//           myStat.myPolyLens.push_back(make_pair(NumTerms(mySPoly),0));
//           T=CpuTime();
//           //mySPoly.smart_dehomog(HIndetIndex);
//           mySPoly.myReduce(myTrueReductors);
//           ++myStat.myNReductions;
//           T-=CpuTime();

//           myStat.myPolyLens.back().second=NumTerms(mySPoly);
//           if (IsZero(mySPoly))
//           {
//             if (VerbosityLevel() >= myStat.myReductionLevel) clog << "0"<<flush;
//             ++myStat.myUseless;
//           }
//           else
//           {
//             if (VerbosityLevel() >= myStat.myReductionLevel)
//               clog<<LPPForDiv(mySPoly)<<"+..<"<<len(myGB)+1<<"> Len="
//                   <<myStat.myPolyLens.back().second<<flush;
//             SPolyPredDeg = wdeg(mySPoly);
//             mySPoly.smart_dehomog(HIndetIndex);

//             if (SPolyPredDeg!=wdeg(mySPoly))
//             {
// //VARIABLE SET NEVER USED              satjumped=true;
//               ++myStat.myPolyDHed;
//               myStat.myDegDH += ConvertTo<long>((SPolyPredDeg-wdeg(mySPoly))[0]);
//               if (VerbosityLevel() >= myStat.myPolyDHLevel && false)
//                 clog << "\nPOLY DEHOMOG SUCCESSFULLY    Deg "
//                      << SPolyPredDeg << "-->" << wdeg(mySPoly)
//                      << LPPForDiv(mySPoly) << endl;
//             }
//             ++myStat.myUseful;
//             //if (!IsTrueGCDDomain(CoeffRing(mySPoly)) && (NumTerms(mySPoly)<=2))
//             //if (NumTerms(mySPoly)<=2)
//             //if (!satjumped)
//             //myTrueReductors.SuperInterreduce(mySPoly);
//             //myTrueReductors.interreduce(mySPoly);
//             if ((VerbosityLevel() >= myStat.myReductionLevel && (NumTerms(mySPoly)<=2)))
//               //if (NumTerms(mySPoly)<=2)
//               clog<<"EASY REDUCTOR FOUND  LEN="<<NumTerms(mySPoly)<<" DEG="<<wdeg(mySPoly)<<endl;
//           } //  else
//         } // else -if not coprime
//         myPairs.erase(myPairs.begin());// erase the used gpair
//         if (VerbosityLevel() >= myStat.myReductionLevel)
//           clog <<" time="<<-T<<" \n"<<flush;
//         if (!IsZero(mySPoly)) myUpdateBasisAndPairs();
//       };//end if minpairs
//     };//end main while
//     //STAT
//     if (myFirstPairWDeg!=myOldDeg)
//       myStat.myDegByDeg.push_back(DegStats(myOldDeg,0,0,0,0,0,1,0));
//     else myStat.myUpgradeDegStats(myOldDeg,0);
// //     VerboseLog VERBOSE("myDoSATMixGBasis");
// //     if (myStat.myGetLevel()>-1)
// //     {
// //       VERBOSE(3) << "number of reductions = " << myStat.myNReductions << endl;
// //       VERBOSE(3) << "# GBasis = " << myStat.myUseful-myStat.myPolyDeleted<< endl;
// //     }
//      myStat.myStampa(clog);
//   } //  end myDoSATMixGBasis


//   // ANNA: not called
//   void GReductor::Rebuild(const PolyList& thePL)
//   {
//     clog << "Rebuild:begin " <<endl;
//     mySPoly=GPoly(myGRingInfoValue);
//     clog << "Rebuild:before for " <<endl;
//     myFirstPairWDeg=degree(GradingDim(myGRingInfoValue.myNewSPR()));
//     myOldDeg=degree(GradingDim(myGRingInfoValue.myNewSPR()));
//     //myStat(len(thePL),myStat.myGetLevel());
//     myPrepared=false;
//     myAgeValue = 0;
//     myStat.myNReductions=0;
//     //IsDynamicAlgorithm unchanged
//     myWrongLPPFoundValue=false;
//     //myCriteria unchanged
//     myPairs.clear();
//     myTrueReductors.myClear();
//     myGB.clear();
//     myPolys.clear();

//     for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();it++)
//     {
//       if (!IsZero(*it))
//         myPolys.push_back(GPoly(*it,myGRingInfoValue));
//     }
//     myPolys.sort(BoolCmpLPPGPoly);
//     clog << "Rebuild:end " <<endl;
//   } // Rebuild


//////////////////////////// Embedding/DeEmbedding /////////////////////


  ModOrdTypeForcing ModuleOrderType(const FreeModule& M)
  {
    if (IsOrdPosn(ordering(M))) return WDegTOPos;
    if (IsWDegPosnOrd(ordering(M))) return WDegPosTO;
    return PosWDegTO;
  } // ModOrdType


  FreeModule owner(const VectorList& theVL)
  {
    CoCoA_ASSERT(!theVL.empty());
    return owner(*theVL.begin());
  } // owner



// similar functions now moved into FreeModule.C
// FreeModule MakeNewFreeModuleForSyz2(const GPolyList& theGPL)
// {
//   if (theGPL.empty())
//     return NewFreeModule(RingQQ(),1);
//   const SparsePolyRing OldP=owner(theGPL);
//   std::vector<degree> InputShifts;
//   for (GPolyList::const_iterator it=theGPL.begin(); it!=theGPL.end(); ++it)
//      InputShifts.push_back(wdeg(*it));

//   return NewFreeModule(OldP, NewWDegPosnOrd(ordering(PPM(OldP)), InputShifts));
// } // MakeNewFreeModuleForSyz


// returns the poly ring equivalent with OldP^2, same grading
// The ordering is WDegPosnOrd if MOType==NoForcing or MOType
  SparsePolyRing MakeNewPRingForSimpleEmbeddingPosFirst(const SparsePolyRing& OldP,
                                                        bool HomogInput)
  {
    if (HomogInput)
      return MakeNewPRingForSimpleEmbedding(OldP, WDegPosTO);
    else
      return MakeNewPRingForSimpleEmbedding(OldP,PosWDegTO);
  }


  SparsePolyRing MakeNewPRingForSimpleEmbedding(const SparsePolyRing& OldP,
                                                ModOrdTypeForcing MOType)
  {
    std::vector<degree> InputShifts;
    {
      degree tmp(GradingDim(OldP));
      InputShifts.push_back(tmp);
      InputShifts.push_back(tmp);
    }
    const FreeModule FM=NewFreeModule(OldP, InputShifts, WDegPosnOrd);
    return MakeNewPRingFromModule(FM, MOType);
  } // MakeNewPRingForSimpleEmbedding


  SparsePolyRing MakeNewPRingForSimpleEmbedding(const SparsePolyRing& OldP)
  {
    return MakeNewPRingForSimpleEmbedding(OldP,NoForcing);
  } // MakeNewPRingForSimpleEmbedding


  SparsePolyRing MakeNewPRingFromModule(const FreeModule& FM)
  {
    return MakeNewPRingFromModule(FM,NoForcing);
  } // MakeNewPRingFromModule


  SparsePolyRing MakeNewPRingFromModulePosFirst(const FreeModule& FM,
                                                bool HomogInput)
  {
    if (HomogInput)
      return MakeNewPRingFromModule(FM, WDegPosTO);
    else
      return MakeNewPRingFromModule(FM, PosWDegTO);
  }


// This is OK for the non-homogeneous case
// For the homogenous case, PosTo this is inefficient, since
// the Deg rows in the To part are useless.
  SparsePolyRing MakeNewPRingFromModule(const FreeModule& FM,
                                        ModOrdTypeForcing MOType)
  {
    const ModuleOrdering MTO = ordering(FM);
    const SparsePolyRing OldP=RingOf(FM);
    const long NumOldInds=NumIndets(OldP);
    long GrDim;
    if (MOType==PosWDegTO)
      GrDim=0;// Set simple sugar on
    else
      GrDim=GradingDim(OldP);
    const long NumNewInds=NumOldInds+GrDim+1;

    ConstMatrixView OldOrdOMat = OrdMat(OldP);

    matrix NewOrdMat(NewDenseMat(RingZZ(), NumNewInds, NumNewInds));
    ////std::clog<<"NewOrdMat starts as "<<NewOrdMat<<std::endl;
    if (MOType == NoForcing)  MOType = ModuleOrderType(FM);

    switch (MOType)
    {
    case PosWDegTO:
      // Setting the module component ordering
      SetEntry(NewOrdMat, 0, NumNewInds-1, 1); 	
      // Part common to IsWDegPosnOrd and IsOrdPosn
      // Setting the Grading: the OldGrading		
      for (long i=0; i < GrDim+1; ++i)			
        for (long j=0; j < NumOldInds; ++j)		
          SetEntry(NewOrdMat, i+1, j, OldOrdOMat(i,j));
      // Setting the Grading: the NewGrading		
      for (long i=1; i < GrDim+1; ++i)			
        SetEntry(NewOrdMat, i, i+NumOldInds-1, 1); 	
      // Setting the TO ordering
      for (long i=GrDim; i < NumOldInds; ++i)
        for (long j=0; j < NumOldInds; ++j)
          SetEntry(NewOrdMat, i+1, j, OldOrdOMat(i,j));
      break;

    case WDegTOPos:
      // Part common to IsWDegPosnOrd and IsOrdPosn
      //clog<<"MakeNewPRingFromModule:case OrdPosn"<<endl;
      // Setting the Grading: the OldGrading		
      for (long i=0; i < GrDim; ++i)			
        for (long j=0; j < NumOldInds; ++j)		
          SetEntry(NewOrdMat, i, j, OldOrdOMat(i,j));
      // Setting the Grading: the NewGrading		
      for (long i=0; i < GrDim; ++i)			
        SetEntry(NewOrdMat, i, i+NumOldInds, 1);
      //clog<<"MakeNewPRingFromModule:the matrix graded "<<NewOrdMat<<endl;
      // Setting the TO	
      for (long i=GrDim; i < NumOldInds; ++i)
        for (long j=0; j < NumOldInds; ++j)
          SetEntry(NewOrdMat, i, j, OldOrdOMat(i,j));
      //clog<<"MakeNewPRingFromModule:the matrix TO "<<NewOrdMat<<endl;
      // Setting the module component ordering
      SetEntry(NewOrdMat, NumNewInds-1-GrDim, NumNewInds-1, 1); 	
      //clog<<"MakeNewPRingFromModule:the matrix TO Pos "<<NewOrdMat<<endl;
      break;

    case WDegPosTO:; // This is the default
    default:
      // Part common to IsWDegPosnOrd and IsOrdPosn
      // Setting the Grading: the OldGrading		
      for (long i=0; i < GrDim; ++i)			
        for (long j=0; j < NumOldInds; ++j)		
          SetEntry(NewOrdMat, i, j, OldOrdOMat(i,j));
    // Setting the Grading: the NewGrading		
    for (long i=0; i < GrDim; ++i)			
      SetEntry(NewOrdMat, i, i+NumOldInds, 1); 	
    // Setting the module component ordering
    SetEntry(NewOrdMat, GrDim, NumNewInds-1, 1);
    // Setting the TO ordering
    for (long i=GrDim; i < NumOldInds; ++i)
      for (long j=0; j < NumOldInds; ++j)
        SetEntry(NewOrdMat, i+1, j, OldOrdOMat(i,j));
    break;
    }
    // Filling the matrix
    for (long i=0; i < GrDim; ++i)
      for (long j=0; j < NumOldInds; ++j)
        SetEntry(NewOrdMat, NumOldInds+i+1, j, OldOrdOMat(i,j));

    //clog<<"MakeNewPRingFromModule:the matrix"<<NewOrdMat<<endl;

    const PPOrdering MatNewOrd = NewMatrixOrdering(NewOrdMat, GrDim);

    const std::vector<symbol> IndetNames = NewSymbols(NumOldInds + GrDim + 1);
//---> for DEBUGGING choose these IndetNames:
//   std::vector<symbol> IndetNames = SymbolRange("x", 0, NumOldInds-1);
//   if ( GrDim==1 ) IndetNames.push_back(symbol("s"));  // indet representing shift
//   else
//     for ( long i=0 ; i<GrDim ; ++i )
//       IndetNames.push_back(symbol("s",i));  // indet representing shift
//   IndetNames.push_back(symbol("e"));  // indet representing module component
//---> for DEBUGGING
    SparsePolyRing NewP(NewPolyRing(CoeffRing(OldP),IndetNames,MatNewOrd));
    return NewP;
  } // MakeNewPRingFromModule


  GPoly EmbedPoly(ConstRefRingElem p,
                  const GRingInfo& theGRI,
                  const long CompIndex)
  {
    //std::vector<RingElem> images;
    //for (long i=0;i!=NumIndets(theGRI.myOldSPR());++i)
    //  images.push_back(indet(theGRI.myNewSPR(), i));
    const RingHom& phi=theGRI.myOldP2NewP();
    return GPoly(phi(p)*theGRI.myE(CompIndex), theGRI);
  } // EmbedPoly


  GPoly EmbedPoly(ConstRefRingElem the_p,
                  const GRingInfo& theGRI,
                  const degree& the_d,
                  const long CompIndex)
  {
    const RingHom& phi=theGRI.myOldP2NewP();
    return GPoly(phi(the_p)*theGRI.myE(CompIndex)*theGRI.myY(the_d), theGRI);
  } // EmbedPoly


  // Embed a true polylist in a GPolyList. Component CompIndex
  GPolyList EmbedPolyList(const PolyList& PL,
                          const GRingInfo& GRI,
                          const long CompIndex)
  {
    GPolyList result;
    if (PL.empty())  return result;
    for (const RingElem& p: PL) result.push_back(EmbedPoly(p, GRI, CompIndex));
    return result;
  } // EmbedPolyList


  GPolyList EmbedPolyListNo0(const PolyList& PL,
                             const GRingInfo& GRI,
                             const long CompIndex)
  {
    GPolyList result;
    if (PL.empty())  return result;
    for (const RingElem& p: PL)
      if (!IsZero(p)) result.push_back(EmbedPoly(p, GRI, CompIndex));
    return result;
  } // EmbedPolyListNo0


  // Embed a true polylist in a GPolyList. Component CompIndex
  GPolyList EmbedPolyList(const PolyList& PL,
                          const GRingInfo& GRI,
                          const degree& d,
                          const long CompIndex)
  {
    GPolyList result;
    if (PL.empty()) return result;
    for (const RingElem& p: PL)
      result.push_back(EmbedPoly(p, GRI, d, CompIndex));
    return result;
  } // EmbedPolyList


/*
This realizes the embedding FM-->NewP of a vector v.
e_i->EY[i] and OldP2NewP gives the RingHom between BaseRing(FM) and NewP
Should be changed to avoid passing FM,NewP.
Is is here only for completeness/debug purposes.
*/
  GPoly EmbedVector(const ModuleElem& v,
                    const GRingInfo& theGRI)			
  {
    return EmbedVector(v, theGRI, 0);
  } // EmbedVector


/*
This realizes the embedding FM-->NewP of a vector v.
e_i->EY[i] and OldP2NewP gives the RingHom between BaseRing(FM) and NewP
Should be changed to avoid passing FM,NewP.
Is is here only for completeness/debug purposes.
*/
  GPoly EmbedVector(const ModuleElem& v,
                    const GRingInfo& theGRI,
                    const long StartingFromCompIndex)			
  {
    RingElem p(theGRI.myNewSPR()),eMax(power(theGRI.myE(),StartingFromCompIndex));
    const RingHom& phi=theGRI.myOldP2NewP();
    for (long i=0; i<NumCompts(owner(v)); ++i)
      p+=phi(v[i])*theGRI.myEY(i)/eMax;
    return GPoly(p, theGRI);
  } // EmbedVector


  GPolyList EmbedVectorList(const VectorList& VL, const GRingInfo& GRI)
  { return EmbedVectorList(VL, GRI, 0); }


  GPolyList EmbedVectorList(const VectorList& VL,
                            const GRingInfo& theGRI,
                            const long StartingFromCompIndex)
  {
    GPolyList outPL;
    if (VL.empty())  return outPL;
    for (const auto& v: VL)
      if (!IsZero(v))
        outPL.push_back(EmbedVector(v, theGRI, StartingFromCompIndex));
    return outPL;
  } // EmbedVectorList		


  GPolyList SyzEmbedVectorList(const VectorList& InputVectorList,
                               const GRingInfo& GRI)
  {
    GPolyList outPL;
    if (InputVectorList.empty())
      return outPL;
    const SparsePolyRing NewP=GRI.myNewSPR();
    outPL=EmbedVectorList(InputVectorList, GRI);
    long k=NumCompts(GRI.myFreeModule());
    if (GRI.myInputAndGrading()==NONHOMOG_GRADING)
      for (GPoly& p: outPL)
      { // Added by JAA 2012/10/11
        RingElem Ek = GRI.myE(k++); // previous k
        p.myAppendClear(Ek);
      }
    else
      for (GPoly& p: outPL)
      { // Added by JAA 2012/10/11
        RingElem EkY = GRI.myE(k++) * GRI.myY(wdeg(p)); // previous k
        p.myAppendClear(EkY);
      }
    return outPL;
  } // SyzEmbedVectorList


  GPolyList SyzEmbedPolyList(const PolyList& InputPolyList,
                             const GRingInfo& theGRI)
  {
    GPolyList outPL;
    if (InputPolyList.empty())  return outPL;
    const SparsePolyRing NewP=theGRI.myNewSPR();
    outPL=EmbedPolyList(InputPolyList, theGRI, 0);
    RingElem SyzPP(NewP); // Gives the right degree to p+E^i
    degree d(GradingDim(NewP));
    long k=1;
    for (GPolyList::iterator it=outPL.begin();it!=outPL.end();++it,++k) // for (auto& g: outPL) BUT MUST ALSO INCREMENT k !!
    {
      SyzPP=one(NewP);
      d=wdeg(*it);
      for (long j=0; j < GradingDim(NewP); ++j)
        SyzPP*=power(theGRI.myY(j),d[j]);
      if (theGRI.myInputAndGrading()==NONHOMOG_GRADING)
      { // Added by JAA 2012/10/11
        RingElem Ek = theGRI.myE(k);
        (*it).myAppendClear(Ek);
      }
      else
      { // Added by JAA 2012/10/11
        RingElem EkSyzPP = theGRI.myE(k)*SyzPP;
        (*it).myAppendClear(EkSyzPP);
      }
    }
    return outPL;
  } // SyzEmbedPolyList


  GPolyList IntEmbedPolyLists(const PolyList& PL1,
                              const PolyList& PL2,
                              const GRingInfo& GRI)
  {
    GPolyList Part1 = EmbedPolyListNo0(PL1, GRI, 0);
    GPolyList Part2 = EmbedPolyListNo0(PL2, GRI, 0);
    GPolyList Part3 = EmbedPolyListNo0(PL2, GRI, 1);
    GPolyList::iterator it3=Part3.begin();
    for (GPolyList::iterator it2=Part2.begin(); it2!=Part2.end(); ++it2,++it3)
      (*it2).myAppendClear(*it3);
    Part2.splice(Part2.begin(), Part1);
    return Part2;
  } // IntEmbedPolyLists


  GPolyList IntEmbedVectorLists(const VectorList& theVL1,
                                const VectorList& theVL2,
                                const GRingInfo& theGRI)
  {
    const long NC = NumCompts(owner(theVL1));
    //const SparsePolyRing NewP=theGRI.myNewSPR();
    GPolyList FirstPart=EmbedVectorList(theVL1, theGRI);
    GPolyList SecondPart=EmbedVectorList(theVL2, theGRI);
    GPolyList ThirdPart=EmbedVectorList(theVL2, theGRI, NC);
    GPolyList::iterator it1=ThirdPart.begin();
    for (GPolyList::iterator it=SecondPart.begin();it!=SecondPart.end();++it,++it1)
      (*it).myAppendClear(*it1);
    SecondPart.splice(SecondPart.begin(),FirstPart);
    return SecondPart;
  } // IntEmbedVectorLists


  // VL2 is a singleton
  GPolyList ColonEmbedVectorLists(const VectorList& theVL1,
                                  const VectorList& theVL2,
                                  const GRingInfo& theGRI)
  {
    GPolyList FirstPart;
    if (theVL1.empty())
      return FirstPart;
    const long NC = NumCompts(owner(theVL1));
    // const SparsePolyRing NewP=theGRI.myNewSPR();
    FirstPart=EmbedVectorList(theVL1, theGRI);
    GPoly GP1=EmbedVector(theVL2.front(), theGRI);
    degree d=wdeg(theVL2.front());
    GPoly GP2=EmbedPoly(one(theGRI.myOldSPR()), theGRI, d, NC);
    GP1.myAppendClear(GP2);
    FirstPart.push_back(GP1);
    return FirstPart;
  } // ColonEmbedVectorLists


  // PL2 is a singleton
  GPolyList ColonEmbedPolyLists(const PolyList& thePL1,
                                const PolyList& thePL2,
                                const GRingInfo& theGRI)
  {
    //std::clog<<"ComputeColon: start "<<endl;
     GPolyList FirstPart;
     if (thePL1.empty())
       return FirstPart;
     const SparsePolyRing NewP=theGRI.myNewSPR();
     FirstPart=EmbedPolyList(thePL1, theGRI, 0);
     //std::clog<<"ComputeColon: FirstPart "<<FirstPart<<endl;
     GPoly GP1=EmbedPoly(thePL2.front(), theGRI, 0);
     //std::clog<<"ComputeColon: GP1 "<<GP1<<endl;
     degree d=wdeg(thePL2.front());
     GPoly GP2=EmbedPoly(one(theGRI.myOldSPR()), theGRI, d, 1);
     GP1.myAppendClear(GP2);
     FirstPart.push_back(GP1);
     //std::clog<<"ComputeColon: end "<<endl;
     return FirstPart;
  } // ColonEmbedPolyLists


  ModuleElem DeEmbedPoly(ConstRefRingElem p, const GRingInfo& theGRI)
  { return DeEmbedPoly(p, theGRI, 0); }


  ModuleElem DeEmbedPoly(ConstRefRingElem p,
                         const GRingInfo& theGRI,
                         const long ComponentsLimit) // the component in p that goes to the 0 component of the output vector v. Lesser components of p go to higher component of v
  {
//std::clog<<"DeEmbedPoly:start"<<endl;
//std::clog<<"DeEmbedPoly:p  "<<p<<endl;
//std::clog<<"DeEmbedPoly:theGRI  "<<theGRI<<endl;
   const SparsePolyRing OldP=theGRI.myOldSPR();
   const SparsePolyRing NewP=theGRI.myNewSPR();
   const FreeModule FM=theGRI.myOutputFreeModule();
   ModuleElem v(FM);

//std::clog<<"DeEmbedPoly:v  "<<v<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myNewP2OldP()  "<<theGRI.myNewP2OldP()<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myOldP2NewP()  "<<theGRI.myOldP2NewP()<<endl;

   const std::vector<ModuleElem>& e = gens(FM);

   RingElem tmp(OldP);
   for (SparsePolyIter i=BeginIter(p); !IsEnded(i); ++i)
   {
//std::clog<<"DeEmbedPoly:ComponentsLimit  "<<ComponentsLimit<<endl;
//std::clog<<"DeEmbedPoly:PP(i)  "<<PP(i)<<endl;
      tmp=theGRI.myNewP2OldP()(monomial(NewP,coeff(i),PP(i)));
      CoCoA_ASSERT(theGRI.myPhonyComponent(PP(i))-ComponentsLimit < NumCompts(FM));
//std::clog<<"DeEmbedPoly:monomial(NewP,coeff(i),PP(i))  "<<monomial(NewP,coeff(i),PP(i))<<endl;
//std::clog<<"DeEmbedPoly:tmp  "<<tmp<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myPhonyComponent(PP(i))-ComponentsLimit  "<<theGRI.myPhonyComponent(PP(i))-ComponentsLimit<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myPhonyComponent(PP(i)) "<<theGRI.myPhonyComponent(PP(i))<<endl;
//std::clog<<"DeEmbedPoly:theGRI.myComponent(PP(i)) "<<theGRI.myComponent(PP(i))<<endl;
      v+=tmp*e[theGRI.myPhonyComponent(PP(i))-ComponentsLimit]; // reversed for coco4compatibility
   }
//std::clog<<"DeEmbedPoly:v "<<v<<endl;
//std::clog<<"DeEmbedPoly:end"<<endl;
   return v;
  } // DeEmbedPoly


  VectorList DeEmbedPolyList(const PolyList& PL, const GRingInfo& theGRI)
  { return DeEmbedPolyList(PL, theGRI, 0); }


  VectorList DeEmbedPolyList(const PolyList& PL,
                             const GRingInfo& theGRI,
			     const long ComponentsLimit) // Poly whose LPP has last var exponent bigger than this number disappear on DeEmbedding

  {
    VectorList VL;
    if (PL.empty())
      return VL;
    for (const RingElem& g: PL)
//    PolyList::const_iterator it;
//    for (it=PL.begin();it!=PL.end();++it)
      if (theGRI.myComponent(LPP(g)) <= theGRI.myComponent(ComponentsLimit))
        VL.push_back(DeEmbedPoly(g,theGRI,ComponentsLimit));
    return VL;
  } // DeEmbedPolyList


   RingElem DeEmbedPolyToP(ConstRefRingElem p,
                           const GRingInfo& theGRI) // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding
   {
     return theGRI.myNewP2OldP()(p);
   } // DeEmbedPolyToP


  PolyList DeEmbedPolyListToPL(const PolyList& PL,
                               const GRingInfo& theGRI,
                               const long ComponentsLimit) // Poly whose LPP has phony last var exponent lesser than this number disappear on DeEmbedding
  {
    VerboseLog VERBOSE("DeEmbedPolyListToPL");
    PolyList outPL;
    if (PL.empty())
      return outPL;
    for (const RingElem& g: PL)
      if (theGRI.myComponent(LPP(g)) <= theGRI.myComponent(ComponentsLimit)) // Copies. Copying disappears when I work with GPolys
      {
        VERBOSE(100) << "LPP(g) = " << LPP(g) << std::endl; /////////////////
//        if (IsZero(g)) continue; // cannot happen, after LPP(g)
        outPL.push_back(DeEmbedPolyToP(g,theGRI));
        // if ( IsConstant(outPL.last()) ) // redmine #1647 ////' doesn't happen
        // {
        //   outPL.clear();
        //   VERBOSE(99) << "g is constant" << std::endl; /////////////////
        //   outPL.push_back(theGRI.myOldSPR()->myOne());
        //   break;
        // }
      }
    return outPL;
  } // DeEmbedPolyList
			
      				
   void monic(PolyList& PL)
   {
     if (PL.empty())
       return;
//     const SparsePolyRing P(owner(PL));
//     PPMonoidElem tmp(PPM(P));
     for (RingElem& g: PL)
//     for (PolyList::iterator it=PL.begin();it!=PL.end();++it)
       g = monic(g);///???       g /= monomial(P,LC(g),tmp);
   }


    void PolyList2GPolyListClear(PolyList& thePL,
                                 GPolyList& theGPL,
                                 const GRingInfo& theGRI)
    {
      theGPL.clear();
      for (const RingElem& g: thePL)
//      for (PolyList::iterator it=thePL.begin();it!=thePL.end();++it)
        theGPL.push_back(GPoly(g, theGRI, clear));
      thePL.clear();
    } // PolyList2GPolyListClear


    void GPolyList2PolyListClear(GPolyList& theGPL,PolyList& thePL)
    {
      thePL.clear();
      if (theGPL.empty())
        return;
      const SparsePolyRing SPR=owner(theGPL);
      for (GPoly& g: theGPL)
//      for (GPolyList::iterator it=theGPL.begin();it!=theGPL.end();++it)
      {
        thePL.push_back(one(SPR));
        swap(thePL.back(),g.myPolyValue);
      }
      theGPL.clear();
    } // GPolyList2PolyListClear


    void PolyList2GPolyList(const PolyList& thePL,
                             GPolyList& theGPL,
                             const GRingInfo& theGRI)
    {
      theGPL.clear();
      for (const RingElem& g: thePL)
//      for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
        theGPL.push_back(GPoly(g, theGRI));
    } // PolyList2GPolyList


    void GPolyList2PolyList(const GPolyList& theGPL,PolyList& thePL)
    {
      thePL.clear();
      if (theGPL.empty())
        return;
//      const SparsePolyRing SPR=owner(theGPL);
      for (const GPoly& g: theGPL)
//      for (GPolyList::const_iterator it=theGPL.begin();it!=theGPL.end();++it)
        thePL.push_back(g.myPolyValue);
    } // GPolyList2PolyListClear


     // Just transform the PolyList in a GPolyList
    GPolyList EmbedPolyList(PolyList& thePL, const GRingInfo& theGRI)
    {
      GPolyList GPL;
      PolyList2GPolyListClear(thePL, GPL, theGRI);
      return GPL;
    }


    void DeEmbedPoly(ModuleElem& theOutputV,
                      const GPoly& the_gp)
    {
      DeEmbedPoly(theOutputV,the_gp,0l);
    } // DeEmbedPoly theOutputP

    // some copies are unavoidable when deembedding
    void DeEmbedPoly(ModuleElem& theOutputV,
                      const GPoly& the_gp,
                      const long ComponentsLimit) // the component in p that goes to the 0 component of the output vector v. Lesser components of p go to higher component of v
    {
      const FreeModule FM=owner(theOutputV);
      const GRingInfo& GRI=the_gp.myGRingInfo();
      theOutputV=zero(FM);
      const std::vector<ModuleElem>& e = gens(FM);
      const RingHom& phi=GRI.myNewP2OldP();
      for (SparsePolyIter i=BeginIter(the_gp.myPoly()); !IsEnded(i); ++i)
        theOutputV+=phi(monomial(GRI.myNewSPR(),coeff(i),PP(i)))*
                    e[GRI.myPhonyComponent(PP(i))-ComponentsLimit]; // reversed for cocoa4compatibility
    } // DeEmbedPoly  theOutputV


    void DeEmbedPolyList(VectorList& theOutputVL,
                          const GPolyList& theGPL)
    {
      DeEmbedPolyList(theOutputVL,theGPL,0);
    } // DeEmbedPolyList theOutputVL


    void DeEmbedPolyList(VectorList& theOutputVL,
                         const GPolyList& theGPL,
			 const long ComponentsLimit) // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding
    {
      theOutputVL.clear();
      if (theGPL.empty())
        return;
      const FreeModule FM=owner(theOutputVL);
      const GRingInfo& GRI(GetGRI(theGPL));
      ModuleElem tmp(FM);
      for (const GPoly& g: theGPL)
//      GPolyList::const_iterator it;
//      for (it=theGPL.begin();it!=theGPL.end();++it)
      if (GRI.myComponent(LPPForDiv(g)) <= GRI.myComponent(ComponentsLimit))			
       {
         DeEmbedPoly(tmp,g,ComponentsLimit);
         theOutputVL.push_back(tmp);
       }
    } // DeEmbedPolyList theOutputVL


    void DeEmbedPoly(RingElem& theOutputP,
                      const GPoly& the_gp) // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding
    {
      theOutputP=the_gp.myGRingInfo().myNewP2OldP()(the_gp.myPoly());
    } // DeEmbedPoly theOutputP


    void DeEmbedPolyList(PolyList& theOutputPL,
                         const GPolyList& theGPL,
			 const long ComponentsLimit) // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding
    {
       theOutputPL.clear();
       if (theGPL.empty())
         return;
       const GRingInfo& GRI(GetGRI(theGPL));
       RingElem tmp(GRI.myNewSPR());
       for (const GPoly& g: theGPL)
//       GPolyList::const_iterator it;
//       for (it=theGPL.begin();it!=theGPL.end();++it)
         if (GRI.myComponent(LPPForDiv(g)) <= GRI.myComponent(ComponentsLimit)) // Copies. Copying disappears when I work with GPolys 			
         {
           DeEmbedPoly(tmp,g);
           theOutputPL.push_back(tmp);
         }
    } // DeEmbedPolyList theOutputPL


//     void SyzEmbedGPolyList(GPolyList& theGPL)
//     {
//       const GRingInfo& GRI(GetGRI(theGPL));
//       long k=NumCompts(GRI.myFreeModule());
//       for (GPolyList::iterator it=theGPL.begin() ; it!=theGPL.end() ; ++it,++k)
//       { // Added by JAA 2012/10/11
//         RingElem EkY = GRI.myE(k)*GRI.myY(wdeg(*it));
//         (*it).myAppendClear(EkY);
//       }
//     } // SyzEmbedGPolyList


//     void IntEmbedGPolyList(GPolyList& theGPL1, GPolyList& theGPL2)
//     {
//       CoCoA_ASSERT(GetGRI(theGPL1)==GetGRI(theGPL1));
//       if (theGPL1.empty() || theGPL2.empty()) return;
//       const GRingInfo& GRI(GetGRI(theGPL1));
//       long k = NumCompts(GRI.myFreeModule());

//       for (GPolyList::iterator it=theGPL2.begin();it!=theGPL2.end();++it,++k)
//       { // Added by JAA 2012/10/11
//         RingElem EmbeddedPoly = it->myPoly()*GRI.myE(k);
//         (*it).myAppendClear(EmbeddedPoly);
//       }
//       theGPL1.splice(theGPL1.begin(),theGPL2);
//     } // IntEmbedGPolyList


    void ColonEmbedGPolyList(GPolyList& theGPL, GPoly& the_gp)
    {
      CoCoA_ASSERT(GetGRI(theGPL)==the_gp.myGRingInfo());
      const GRingInfo& GRI(the_gp.myGRingInfo());
      const long NC = NumCompts(GRI.myFreeModule());
      RingElem tmp =  GRI.myE(NC)*GRI.myY(wdeg(the_gp)); // JAA 2012-10-11
      the_gp.myAppendClear(tmp);                         // JAA 2012-10-11
      theGPL.push_back(GPoly(GRI));
      theGPL.back().AssignClear(the_gp);
    } // ColonEmbedGPolyList


    // The ordering is supposed to be Deg compatible
    RingElem homog(ConstRefRingElem the_p, const std::vector<RingElem>& theY)
    {
      const SparsePolyRing SPR=owner(the_p);
      if (GradingDim(SPR) != len(theY))
        CoCoA_THROW_ERROR("incompatible GradingDim", "homog(f, y1_yn");
      RingElem the_hp(SPR);
      RingElem tmp(SPR);
      degree MaxDeg(GradingDim(SPR));
      degree TmpDeg(GradingDim(SPR));
      // Compute the maximum degree
      for (SparsePolyIter it=BeginIter(the_p);!IsEnded(it);++it)
        MaxDeg=top(MaxDeg,wdeg(PP(it)));
      // Homogenizing
      for (SparsePolyIter it=BeginIter(the_p);!IsEnded(it);++it)
      {
        //std::clog<<"Homogenized:tmp"<<tmp<<endl;
        tmp = monomial(SPR,coeff(it),PP(it));
        TmpDeg=wdeg(PP(it));
        for (long i=0; i != len(theY); ++i)
          tmp*=power(theY[i],(MaxDeg[i]-TmpDeg[i]));
        SPR->myAddClear(raw(the_hp), raw(tmp));
      }
      return the_hp;
    } // homog (multihomog)


    void homogenized(ModuleElem& /*the_hv*/,
                     const ModuleElem& /*the_v*/,
                     const GRingInfo& /*theGRI*/)
    {
      CoCoA_THROW_ERROR(ERR::NYI, "homogenized(ModuleElem&, ModuleElem&, const GRingInfo");
    } // homogenized


/*
      const FreeModule NewFM=owner(the_hv);
      the_hv=zero(NewFM);
      const std::vector<ModuleElem>& e = gens(NewFM);
      degree maxdegree=wdeg(the_v);
      for (long i=0;i<NumCompts(owner(the_v));++i)
        for (SparsePolyIter it=BeginIter(the_v[i]); !IsEnded(it); ++it)
          the_hv=the_hv+
                   e[i]*(theGRI.myOldP2NewP())(monomial(theGRI.myNewSPR(),coeff(it),PP(it)))
                   *theGRI.myY(maxdegree-wdeg(PP(it)));
    } // Homogenize
*/

//     std::vector<long> PolyList2IndexList(const PolyList& thePL)
//     {
//       std::vector<long> outPL;
//       std::vector<long> tmp;
//       long i=0;

// //      for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
//       for (const auto& p: thePL)
//       {
//         SparsePolyIter it_p=BeginIter(p);
//         exponents(tmp,PP(it_p));
//         i=0;
//         while (i!=len(tmp) && tmp[i]==0) ++i;
//         CoCoA_ASSERT(i!=len(tmp));
//         outPL.push_back(i);
//       }
//       return outPL;
//     }


//   PPMonoidElem IndexList2PPMonoidElem(const PPMonoid& thePPM,
//                                       const std::vector<long>& theIL)
//   {
//     PPMonoidElem t(thePPM);
//     for (std::vector<long>::const_iterator it=theIL.begin();it!=theIL.end();++it)
//       t*=indet(thePPM,*it);
//     return t;
//   } // IndexList2PPMonoidElem


  std::vector<long> PPMonoidElem2IndexList(ConstRefPPMonoidElem thePP)
  {
    std::vector<long> IndexList;
    std::vector<long> tmp = exponents(thePP);
    for (long i=0; i < len(tmp); ++i)
      if (tmp[i]!=0)  IndexList.push_back(i);
    return IndexList;
  } // PPMonoidElem2IndexList


  SparsePolyRing MakeElimRingFromOld(const SparsePolyRing& theOldP,
                                     const std::vector<long>& IndetsToElim,
                                     const bool IsHomog)
  {
    ConstMatrixView weights = GradingMat(theOldP);

    matrix NewOrdMat = ElimMat(IndetsToElim, weights);
    long NewGrDim = 0;
    if (GradingDim(theOldP) != 0 && IsHomog)
    {
      NewOrdMat = ElimHomogMat(IndetsToElim, weights);
      NewGrDim = GradingDim(theOldP);
    }

    const PPOrdering& NewOrd = NewMatrixOrdering(NewOrdMat, NewGrDim);
    vector<symbol> IndetNames = NewSymbols(NumIndets(theOldP));
    SparsePolyRing NewP=NewPolyRing(CoeffRing(theOldP),IndetNames,NewOrd);
    return NewP;
  } // MakeElimRingFromOld


    bool IsHomog(const PolyList& PL)
    {
// HINT: CAN USE std::all_of here
      // morally:  return find_if(v.begin(), v.end(), not1(IsHomog)) == v.end();
      for (const RingElem& p: PL)
        if (!IsHomog(p))  return false;
      return true;
//    We *DO NOT USE* STL algorithm because ptr_fun fails when fun has arg of reference type.
//       return find_if(v.begin(), v.end(),
//                      not1(ptr_fun(static_cast<bool(*)(ConstRefRingElem)>
//                                   (CoCoA::IsHomog))))
//         == v.end();
    }



    bool IsHomog(const VectorList& VL)
    {
      for (const auto& v: VL)
        if (!IsHomog(v))  return false;
      return true;
    }


/////////////////
/*
    std::vector<degree> DegStructure(ConstRefRingElem the_p)
    {
      std::vector<degree> L;
      for (SparsePolyIter it=BeginIter(the_p);!IsEnded(it); ++it)
        L.push_back(wdeg(PP(it)));
      return L;
    } // DegStructure

   std::vector<std::vector<degree> > DegStructure(const ModuleElem& the_v)
    {
      std::vector<std::vector<degree> > L;
      std::vector<degree> tmp;
      const FreeModule FM=owner(the_v);
      const std::vector<degree> S=shifts(FM);
      for (long i=0;i<NumCompts(owner(the_v));++i)
      {
        tmp.clear();
        for (SparsePolyIter it=BeginIter(the_v[i]);!IsEnded(it); ++it)
            tmp.push_back(wdeg(PP(it))+S[i]);
        L.push_back(tmp);
      }
      return L;								
    } // DegStructure

    std::vector<std::vector<degree> > DegStructure(const PolyList& thePL)
    {
      std::vector<std::vector<degree> > L;
      for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
          L.push_back(DegStructure(*it));
      return L;
    } // DegStructure



    std::vector<std::vector<std::vector<degree> > > DegStructure(const VectorList& theVL)
    {
      std::vector<std::vector<std::vector<degree> > > L;
      for (VectorList::const_iterator it=theVL.begin();it!=theVL.end();++it)
          L.push_back(DegStructure(*it));
      return L;
    } // DegStructure

*/

    PolyList MakePolyList(ConstRefRingElem the_p)
    {
      PolyList L;
      L.push_back(the_p);
      return L;
    } // MakePolyList

    VectorList MakeVectorList(const ModuleElem& the_v)
    {
     VectorList L;
      L.push_back(the_v);
      return L;
    } // MakeVectorList


  PolyList WithoutDenominators(const PolyList& PL, SparsePolyRing Rx)
  {
    if (PL.empty()) return PL;
    CoCoA_ASSERT(HasUniqueOwner(PL));
    const SparsePolyRing Kx(owner(PL));
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(CoeffRing(Kx)));
    CoCoA_ASSERT(BaseRing(CoeffRing(Kx)) == CoeffRing(Rx));
    CoCoA_ASSERT(ordering(PPM(Kx)) == ordering(PPM(Rx)));
    PolyList outPL;
    std::vector<long> expv(NumIndets(Rx));
    RingElem f(Rx);
    for (const RingElem& p: PL)
    {
      const RingElem d(CommonDenom(p));
      f = 0;
      for (SparsePolyIter i=BeginIter(p); !IsEnded(i); ++i)
        PushBack(f, (d/den(coeff(i))*num(coeff(i))), exponents(expv,PP(i))); // same PPO
        //        f += monomial(Rx, (d/den(coeff(i))*num(coeff(i))), e);
      outPL.push_back(f);
    }
    return outPL;
  } // WithoutDenominators


  PolyList WithDenominator1Hom(const PolyList& PL, SparsePolyRing Kx)
  {
    CoCoA_ASSERT(IsFractionField(CoeffRing(Kx)));
    if (PL.empty()) return PL;
    CoCoA_ASSERT(HasUniqueOwner(PL));
    const SparsePolyRing Rx(owner(PL));
    CoCoA_ASSERT(BaseRing(CoeffRing(Kx)) == CoeffRing(Rx));
    CoCoA_ASSERT(ordering(PPM(Kx)) == ordering(PPM(Rx)));
    PolyList outPL;
    // outPL.reserve(len(PL)); // if PL is a vector
    const RingHom EmbedIntoFrF = EmbeddingHom(CoeffRing(Kx));
    const RingHom EmbedCoeff   = CoeffEmbeddingHom(Kx);
    const RingHom phi = PolyRingHom(Rx, Kx, EmbedCoeff(EmbedIntoFrF), indets(Kx));
    //    transform(PL.begin(), PL.end(), back_inserter(outPL), phi);
    return phi(PL);
  } // WithDenominator1Hom

// void ReadInt(std::istream& in, int& the_int,SkipTagType ST)
// {
//   if (ST == GetTag) SkipTag(in, "<int>");
//   in >> the_int;
//   SkipTag(in, "</int>");
// } // ReadInt

//////////////////

} // end namespace cocoa


// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpGReductor.C,v 1.132 2024/09/04 15:44:01 bigatti Exp $
// $Log: TmpGReductor.C,v $
// Revision 1.132  2024/09/04 15:44:01  bigatti
// Summary: renamed Component into component
//
// Revision 1.131  2024/05/30 14:45:22  bigatti
// Summary: fixed truncation (wdeg check moved from myReduceCurrentSPoly into myDoGBasis);
//
// Renamed myFirstPairWDeg into myIncomingWDeg
//
// Revision 1.130  2024/05/28 13:56:34  bigatti
// Summary: fixed initialization of myTruncDegValue
//
// Revision 1.129  2024/05/28 13:02:38  bigatti
// Summary: added code for truncation in myDoGBasis;
// changed myGBasis(G)/myMinGens(G)  into  myCopyGBasis(G)/myMinGens(G);
// changed myCurrentPairDeg  into myFirstPairWDeg;
// other minor changes in view of major cleaning (removing dynamic code)
//
// Revision 1.128  2024/05/03 12:17:47  abbott
// Summary: Added comment (see also redmine 1346)
//
// Revision 1.127  2024/03/22 15:58:21  bigatti
// Summary: in myFinalizeGBasis checking just the LPP instead of making reduction (to 0)
//
// Revision 1.126  2024/03/20 15:35:43  bigatti
// Summary: Simplified MakeElimRingFromOld by using GradingMat
//
// Revision 1.125  2024/03/15 15:47:15  bigatti
// Summary: added comment AllCols as a reminder for new function FirstRows
//
// Revision 1.124  2024/02/15 15:34:08  bigatti
// Summary: added verbosity in DeEmbedPolyListToPL
//
// Revision 1.123  2024/02/05 11:28:35  bigatti
// Summary: commented out using std::find_if
//
// Revision 1.122  2023/03/13 21:06:40  abbott
// Summary: New iterator "for" syntax; some cruft removed where the "for" has already been changed.  Redmine 1346
//
// Revision 1.121  2022/10/11 15:39:35  bigatti
// Summary: using new interface of exponents returning vector
//
// Revision 1.120  2022/09/14 19:18:31  abbott
// Summary: Renamed local var
//
// Revision 1.119  2022/02/18 15:39:55  abbott
// Summary: Used lambda fn to avoid "deprecated" warning from clang
//
// Revision 1.118  2022/02/18 14:12:00  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.117  2021/10/04 09:01:20  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.116  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.115  2020/02/18 11:31:28  abbott
// Summary: redmine 1346: new for loop syntax  (perhaps only partly done)
//
// Revision 1.114  2020/02/12 10:53:17  bigatti
// -- removed & after auto
//
// Revision 1.113  2020/02/11 18:13:59  bigatti
// -- converted many loops in new "foreach" syntax
//
// Revision 1.112  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.111  2020/02/11 16:12:20  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.110  2020/02/11 16:02:51  bigatti
// -- removed calls to clog
//
// Revision 1.109  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.108  2018/07/04 13:13:11  bigatti
// -- added verbosity(70) for GBasisRealSolve
//
// Revision 1.107  2018/06/27 12:15:19  abbott
// Summary: Renamed RealSolveCore to RealSolve
//
// Revision 1.106  2018/06/27 08:50:39  abbott
// Summary: Revised to work with new CpuTimeLimit
//
// Revision 1.105  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.104  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.103  2018/05/18 12:24:47  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.102  2018/05/17 15:53:54  bigatti
// -- renamed VectorOperations --> VectorOps
// -- renamed MatrixOperations --> MatrixOps
// -- sorted #includes
// -- added include SparsePolyIter
//
// Revision 1.101  2018/04/04 12:35:38  bigatti
// -- fixed GPolyList EmbedVectorList: now removing 0-generators
//
// Revision 1.100  2018/03/29 09:40:03  bigatti
// -- added check for interrupt in "ReduceActiveLM"
//    (need proper speed testing)
//
// Revision 1.99  2018/03/08 17:04:29  bigatti
// -- myUpdateBasisAndPairs: now terminating immediately if ideal is (1)
//
// Revision 1.98  2018/02/03 04:58:02  bigatti
// -- catipalized a verbosity message (NEW DEG)
//
// Revision 1.97  2017/12/18 13:12:39  abbott
// Summary: Renamed include files (underscore becomes minus)
//
// Revision 1.96  2017/12/01 21:35:08  abbott
// Summary: Now uses RealRadical
//
// Revision 1.95  2017/12/01 17:23:34  bigatti
// -- added EmbedPolyListNo0
// -- some spaces and variable names
//
// Revision 1.94  2017/11/29 17:45:05  bigatti
// -- added GBasisRealSolveCore
// -- changed myUpdateBasisAndPairs: now can take argument
// -- some spaces
//
// Revision 1.93  2017/11/24 17:46:39  bigatti
// -- renamed GBasisSelfSat --> GBasisSelfSatCore
// -- added GBasisSelfSat in cpkg5
//
// Revision 1.92  2017/11/23 12:37:02  bigatti
// -- added GBasisSelfSat
//
// Revision 1.91  2017/09/06 11:56:30  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.90  2017/05/11 09:29:16  bigatti
// -- removed commented old code for catching interrupt
//
// Revision 1.89  2017/05/02 12:16:01  bigatti
// -- removed NewSparsePolyRing (now NewPolyRing does the same)
//
// Revision 1.88  2017/04/28 13:54:19  bigatti
// -- minor cleaning
// -- renamed AllGraded --> HOMOG, AllAffine --> NOGRADING
//
// Revision 1.87  2017/04/26 12:55:08  bigatti
// -- some cleaning: len --> NumTerms, f.IsActive() --> IsActive(f), ..
// -- changed basic verbosity a bit
//
// Revision 1.86  2017/04/21 13:45:30  bigatti
// -- fixed bug in myFinalizeGBasis (problems in interreduction when
//    first poly redundant)
//
// Revision 1.85  2017/04/18 09:46:05  bigatti
// -- now using VerbosityLevel
// -- now GRStats store levels (instead of booleans)
// -- now GPairs do not store indices (coded as age(poly))
// -- changed Copyright
//
// Revision 1.84  2017/03/08 15:34:42  bigatti
// -- added verbosity (100, 105)
//
// Revision 1.83  2017/03/02 09:11:30  bigatti
// -- modified verbosity
//
// Revision 1.82  2017/02/24 08:23:43  bigatti
// -- someverbosity
//
// Revision 1.81  2017/02/14 17:34:12  bigatti
// -- added verbosity and CheckForInterrupt in myFinalizeGBasis
//
// Revision 1.80  2017/02/08 17:02:40  abbott
// Summary: Commented out variables set but never used
//
// Revision 1.79  2016/12/05 09:30:49  bigatti
// -- first verbosity entry in myDoGBasis (level 1965)
//
// Revision 1.78  2016/11/18 19:36:07  abbott
// Summary: Simplified call to CheckForInterrupt
//
// Revision 1.77  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.76  2016/09/22 15:33:37  bigatti
// -- renamed HomogElimMat into ElimHomogMat
// -- improved readability for ElimHomogMat/ElimMat (removed auxiliary functions)
//
// Revision 1.75  2016/09/08 14:14:30  bigatti
// -- commented out unused variable
//
// Revision 1.74  2016/01/26 13:52:23  bigatti
// -- just two empty lines
//
// Revision 1.73  2015/12/01 13:34:44  abbott
// Summary: Changed arg order in ElimMat and HomogElimMat; doc is out-of-date!!
//
// Revision 1.72  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.71  2015/06/30 09:55:56  bigatti
// -- added interrupt mechanism
//
// Revision 1.70  2015/05/20 14:44:51  bigatti
// -- renamed AmIModule --> IamModule
//
// Revision 1.69  2015/03/04 10:33:56  bigatti
// -- added myDoGBasisElimFirst
// -- minor fixes
//
// Revision 1.68  2014/10/24 12:52:27  bigatti
// -- only white spaces
//
// Revision 1.67  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.66  2014/07/31 13:10:45  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.65  2014/07/30 14:10:32  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.64  2014/07/22 10:47:33  bigatti
// -- using NewSymbols for elim
//
// Revision 1.63  2014/07/14 15:08:16  abbott
// Summary: Changed include of tmp.H into UtilsTemplate.H
// Author: JAA
//
// Revision 1.62  2014/07/09 14:27:53  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.61  2014/07/08 08:39:48  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.60  2014/07/07 13:06:32  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.59  2014/04/30 16:16:04  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.58  2014/04/17 13:38:57  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.57  2014/04/11 15:44:27  abbott
// Summary: Renamed MatrixArith to MatrixOperations (in includes)
// Author: JAA
//
// Revision 1.56  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.55  2013/12/03 12:09:11  bigatti
// -- now incrementing myAgeValue only if added to GB
//
// Revision 1.54  2013/10/28 13:18:05  bigatti
// -- myCreateSpecialPairs --> myCreateInputPolyPairs
//
// Revision 1.53  2013/06/12 08:52:19  bigatti
// -- added mySetMinimalGen
// -- commented out unused code
//
// Revision 1.52  2013/06/03 09:11:50  bigatti
// renamed ModuleTermOrdering into ModuleOrdering
//
// Revision 1.51  2013/05/27 17:00:21  bigatti
// -- new names for ModuleOrdering
//
// Revision 1.50  2013/03/26 14:56:06  abbott
// Updated the conversion fns (in ptic removed procedure "convert");
// numerous consequential changes.
//
// Revision 1.49  2013/02/21 17:17:13  bigatti
// -- moved NewFreeModuleForSyz to FreeModule
//
// Revision 1.48  2013/02/14 17:34:35  bigatti
// -- cleaned up code for elimination matrices
//
// Revision 1.47  2012/10/24 12:23:27  abbott
// Removed use of std::ptr_fun.
//
// Revision 1.46  2012/10/16 09:54:27  abbott
// Replaced  RefRingElem  by  RingElem&
// Modified many calls to  myAppendClear  (because a temporary value will
// not match a non-const reference in a fn call).
//
// Revision 1.45  2012/10/03 12:26:08  bigatti
// -- changed homogenized --> homog: now returns a RingElem (instead of void)
//
// Revision 1.44  2012/05/28 16:47:23  bigatti
// -- added useful comment for debugging on anonymous symbols
//
// Revision 1.43  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.42  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.41  2012/05/10 14:44:53  abbott
// Employed anonymous symbols in MakeNewPRingFromModule (fixing issue #155).
//
// Revision 1.40  2012/03/30 17:31:04  bigatti
// -- making ordering matrices over QQ
//
// Revision 1.39  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.38  2012/02/08 16:14:49  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.37  2012/01/26 16:49:51  bigatti
// -- using apply instead of transform
//
// Revision 1.36  2011/12/07 16:00:07  bigatti
// -- removed double "#include <functional>"
//
// Revision 1.35  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.34  2011/05/27 12:13:14  bigatti
// -- added CoCoA_ASSERT for empty list of generators
//
// Revision 1.33  2011/05/24 08:16:00  bigatti
// -- only spaces
//
// Revision 1.32  2011/03/22 11:26:53  bigatti
// -- unsigned long --> long
// -- myAge --> myAgeValue
//
// Revision 1.31  2011/03/11 17:39:10  bigatti
// -- changed  unsigned int --> long
//
// Revision 1.30  2011/03/11 15:39:05  bigatti
// -- changes size_t --> long
// -- changes size --> len
//
// Revision 1.29  2011/03/11 14:50:39  abbott
// Changed some size_t into long.
//
// Revision 1.28  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.27  2011/03/08 18:06:33  bigatti
// -- removed useless vector<long> copy in MakeElimRingFromOld
//
// Revision 1.26  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.25  2010/08/06 17:11:02  bigatti
// -- restored owner(const VectorList&)  // should be moved into TmpGTypes
//
// Revision 1.24  2010/08/06 10:28:15  bigatti
// -- commented out myDoAFFGBasis and owner(VL): unused
//
// Revision 1.23  2010/08/05 16:05:36  bigatti
// -- simpler code after new class criteria
//
// Revision 1.22  2010/07/16 09:29:58  bigatti
// -- minor cleaning and coding conventions
//
// Revision 1.21  2010/05/14 09:53:09  bigatti
// -- removed empty ctor for SugarDegree
// -- added marker for SugarDegree(uninitialized)
// -- SugarDegree for GBasis input is initialized by myPrepareGBasis
//
// Revision 1.20  2010/03/23 14:43:07  bigatti
// -- class GRingInfo estracted from TmpGPoly
//
// Revision 1.19  2009/11/26 16:59:19  bigatti
// -- minor cleaning of WithoutDenominators and WithDenominator1Hom
// -- sorted includes
//
// Revision 1.18  2009/10/26 16:29:57  bigatti
// -- minor change in creting a zero GPoly
// -- added printing of sugar when starting a new degree
//
// Revision 1.17  2009/09/22 13:35:55  bigatti
// -- following coding conventions in function names Matrix --> Mat
// -- forced all matrices to be over RingZ
//
// Revision 1.16  2009/05/14 12:30:33  abbott
// Added missing colon.
//
// Revision 1.15  2009/04/27 14:11:22  bigatti
// -- added BuchbergerOpTypeMarker flag
// -- added cocoa_assert
//
// Revision 1.14  2009/03/20 15:33:26  bigatti
// -- just comments and spaces
//
// Revision 1.13  2009/01/30 13:41:50  bigatti
// -- enum instead of bool arguments
//
// Revision 1.12  2008/09/22 15:58:58  bigatti
// -- commented out unused "double T=0.0"
//
// Revision 1.11  2008/09/19 14:08:16  bigatti
// -- modified GRStats (M.Caboara)
//
// Revision 1.10  2008/09/16 15:03:42  bigatti
// -- added LPPForDiv
// -- changed LPP into LPPForOrd
//
// Revision 1.9  2008/07/09 16:10:11  abbott
// Added missing & (for C++ const reference).
//
// Revision 1.8  2008/04/21 11:23:11  abbott
// Separated functions dealing with matrices and PPOrderings into a new file.
// Added matrix norms, and completed adjoint.
//
// Revision 1.7  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.6  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.5  2007/11/09 10:45:52  bigatti
// -- [caboara] preparation for self-saturating algorithm
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.2  2007/03/12 16:31:10  bigatti
// -- fixed iterators issues (M.Abshoff)
//
// Revision 1.1  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.30  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.29  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.28  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.27  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.26  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.25  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.24  2007/02/12 16:44:05  bigatti
// -- restore "GlobalLogput()" in logging info -- to be discussed with Max
//
// Revision 1.23  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.22  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.21  2006/12/04 13:55:54  cocoa
// -- added: sugar for GradingDim > 0  (called wsugar)
//
// Revision 1.20  2006/11/27 14:24:49  cocoa
// -- reorganised #include files
//
// Revision 1.19  2006/11/24 17:39:14  cocoa
// -- reorganized includes of header files
//
// Revision 1.18  2006/11/22 15:33:39  cocoa
// -- changed handling of timings and number of reductions (M.Caboara)
//
// Revision 1.17  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.16  2006/11/14 17:20:23  cocoa
// -- uncommented myDoAFFGBasis for RingWeyl (should be updated soon)
//
// Revision 1.15  2006/11/09 17:36:06  cocoa
// -- bug fix for WDegPos ordering (by M.Caboara)
//
// Revision 1.14  2006/10/06 17:32:26  cocoa
// -- wip: evolution of Groebner Framework (Max)
//
// Revision 1.13  2006/10/06 16:32:06  cocoa
// -- changed: GPoly::SPoly --> GPoly::myAssignSPoly
// -- changed: Len(const GPoly&) --> len(const GPoly&)
// -- added:   poly(const GPoly&)
// -- added:   GPoly::myUpdateLenLPPDegComp()
// -- in reduce.C added functions for computing sugar during reduction
//
// Revision 1.12  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.11  2006/10/06 10:48:34  cocoa
// Removed #include references to GradedFreeModule.H
//
// Revision 1.10  2006/10/06 10:16:26  cocoa
// Replaced two dynamic_casts by calls to the new function IsRingFp.
//
// Revision 1.9  2006/09/27 14:33:50  cocoa
// -- fixed check for RingFp in NewSparsePolyRing
//
// Revision 1.8  2006/08/17 09:32:26  cocoa
// -- added: flags for homogeneous input
//
// Revision 1.7  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.6  2006/07/20 14:43:56  cocoa
// -- just spaces and a CoCoA_ASSERT
//
// Revision 1.5  2006/07/18 11:01:00  cocoa
// -- changed symbol for the names of the "module" indeterminates
// -- added variables for calling ring homomorphisms
// -- moved IsHomogeneous to SparsePolyRing/FreeModule
//
// Revision 1.4  2006/06/20 17:17:45  cocoa
// -- reinserted TotalTime printing: this is to compare timings for
//    copying the result
//
// Revision 1.3  2006/06/13 16:40:55  cocoa
// -- simpler code in CoCoAServer to print length of the result and time
//    (easier for benchmarcks)
//
// Revision 1.2  2006/06/09 13:24:53  cocoa
// -- just a "\t" before printing "TotalTime"
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.24  2006/05/16 09:01:26  cocoa
// -- added "const" to arg of homogenized
// -- changed ordering matrix for modules (compatibility with GradedFreeModule)
//
// Revision 1.23  2006/05/12 17:03:16  cocoa
// -- swapped arguments in homogenized
//
// Revision 1.22  2006/05/12 10:46:53  cocoa
// -- commented out printouts
//
// Revision 1.21  2006/05/11 16:00:22  cocoa
// -- fixed spelling of "homogenize"
//
// Revision 1.20  2006/05/04 14:25:16  cocoa
// -- major cleaning of FreeModule: created GradedFreeModule and moved
//    some code around
//
// Revision 1.19  2006/05/02 14:37:12  cocoa
// -- only changes to logging info (by M.Abshoff)
//
// Revision 1.18  2006/04/26 16:59:27  cocoa
// -- restored printing of "time" (for benchmarks)
//
// Revision 1.17  2006/04/26 09:50:40  cocoa
// -- added GReductor::ourDefaultStatLevel variable to allow CoCoAServer
//    to set statistics level
//
// Revision 1.16  2006/04/21 16:45:17  cocoa
// -- new functions by Max
//
// Revision 1.15  2006/04/11 16:59:52  cocoa
// -- commented out "sort" in myUpdateBasisAndPairs: probably unnecessary
//
// Revision 1.14  2006/04/11 16:44:19  cocoa
// -- added functions for Elim
// -- minor cleanup
//
// Revision 1.13  2006/04/10 17:02:47  cocoa
// -- changed: BCriterion_OK now uses PPWithMask instead of ConstRefPPMonoidElem
//
// Revision 1.12  2006/03/17 18:17:16  cocoa
// -- changed: use of ReductionCog for reduction (major cleanup)
//
// Revision 1.11  2006/03/02 13:45:57  cocoa
// -- changes by Max
//
// Revision 1.10  2006/02/13 14:46:45  cocoa
// -- changes by Max
//
// Revision 1.9  2006/02/13 13:53:04  cocoa
// -- changes by Max (GRingInfo)
//
// Revision 1.8
// date: 2006/01/20 15:43:30;  author: cocoa;  state: Exp;  lines: +6 -3
// -- fixed: use of RefPPMonoidElem and ConstRefPPMonoidElem

// Revision 1.7
// date: 2006/01/18 15:58:20;  author: cocoa;  state: Exp;  lines: +133 -131
// -- new changes by Max

// Revision 1.6  2006/01/17 15:44:56  cocoa
// -- chamges by Max for operations with modules
//
// Revision 1.5  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.4  2005/12/31 12:22:17  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.3  2005/11/17 13:21:54  cocoa
// -- used "zero(myPolyRing)"
//
// Revision 1.2  2005/11/15 13:05:42  cocoa
// Forgot to check these before checking in the new timing functions.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
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
// Revision 1.11  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.10  2004/11/02 15:58:33  cocoa
// -- added some (commented) lines for new DivMask::SingleBitWrapImpl
//
// Revision 1.9  2004/10/29 16:02:52  cocoa
// -- added field myDivMaskImplPtr for creating GPolys with the same
//    DivMask::base
// -- calls to LCMwMask (instead of LCM) for fast divisibility tests
// -- function IsDivisible had wrong semantics --> swapped arguments everywhere
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
// Revision 1.6  2004/05/27 16:52:12  cocoa
// -- removed ";" at the end of function bodies (g++ gives error on them)
//
// Revision 1.5  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.4  2004/03/04 11:37:18  cocoa
// -- updated code for Borel reductors:
//    ReductorData fields: myGPoly->myGPolyPtr;  new: myCount, IamBorelUpdated
//    myBorelReductors is now in Reductors (was in GReductor)
//    Reductors: field: IhaveBorelReductors;  type: enum UseBorelMarker
//
// Revision 1.3  2003/10/09 12:48:17  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/10/01 10:35:31  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.15  2003/09/22 17:30:10  bigatti
// - new field myOrdPoly to order the GPairs
// - ==> not ordered list of pairs which are completed and sorted later
//
// Revision 1.14  2003/06/23 17:10:00  abbott
// Minor cleaning prior to public release.
// Improved the include directives,
//
// Revision 1.13  2003/05/29 16:56:56  bigatti
// - change: myRingSpecialIndex is now an int
// - added:  flag for criterion (to be disables for non-comm Groebner)
// - fix:    all zero polynomials are now removed from input list
// - change: each GReductor may now set its own statistics level
//
// Revision 1.12  2003/05/28 14:29:43  bigatti
// - new code for modules
//
// Revision 1.11  2003/05/14 17:09:28  bigatti
// - new functions for "BorelReductorsPolys" and saturating algorithm
// - optimization of PPMonoidElem function calls
//
// Revision 1.10  2002/11/15 16:46:40  bigatti
// - clear code after new PP.H (with DivMask)
//
// Revision 1.9  2002/09/19 17:23:32  bigatti
// - Cleaner code based on PolyRing
//
// Revision 1.8  2002/05/13 11:48:36  bigatti
// - new data structure for "Reductors"
//
// Revision 1.7  2002/04/29 17:06:29  bigatti
// - removed "true" from "BuildNewPairs" [if (true||myStat.print_poly_deleted)..]
//
// Revision 1.6  2002/04/15 17:18:59  bigatti
// - Max's new code
//
// Revision 1.5  2002/04/09 14:20:27  bigatti
// - SPoly now takes a GPair as argument
// - CKR criterion
//
// Revision 1.4  2002/01/31 17:23:39  bigatti
// - interreduce
//
// Revision 1.3  2002/01/11 10:34:34  bigatti
// - removed all "--myGB.end()"
//
// Revision 1.2  2001/12/12 18:16:39  bigatti
// - new structure of reduction
//
// Revision 1.1  2001/12/05 13:25:45  bigatti
// Initial revision
//

// Ordinamento delle coppie - cambiare ?
// inline e controllo efficienza
// Sugar
// liste-vettori   poly con rc
// risistemazioni delle stampa
// GPairs con i cofattori
// poly const * f - per theGPoly e per Pairs, se non c'e' il rc
// current paur? SPoly? or the actual structure using the last of .. .. vedere se c'e' rc
// ordinamento dei riduttori
// le statistiche deg by deg sono appiccicate con lo sputo, in quanto non
//   matchano con la struttura.
