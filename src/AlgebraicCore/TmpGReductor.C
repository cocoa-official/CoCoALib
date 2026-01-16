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
    myTruncDegValue = ourNoTruncValue;
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
    myTruncDegValue = ourNoTruncValue;

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
  void GReductor::myCopyGBasis(PolyList& outG)
  {
    outG.clear();
    for (const auto& ptr: myGB)
      if (IsActive(*ptr))  outG.push_back((*ptr).myPoly());
  }


  void GReductor::myCopyMinGens(PolyList& outG)
  {
    outG.clear();
    for (const auto& ptr: myGB)
      if (IsMinimalGen(*ptr))  outG.push_back((*ptr).myPoly());
  }


 // This procedure may be substituted by a transform_if
  void GReductor::myCopyGBasis_GPoly(GPolyList& outG)
  {
    outG.clear();
    for (const auto& ptr: myGB)
      if (IsActive(*ptr))  outG.push_back((*ptr));
  }
  

// This procedure may be substituted by a transform_if
  void GReductor::myGBasisClear(GPolyList& outG)
  {
    outG.clear();
    for (auto& ptr: myGB)
      if (IsActive(*ptr))
      {
        GPoly Zero(myGRingInfoValue);
        outG.push_back(Zero);
        outG.back().AssignClear(*ptr);
      }
    myGB.clear();
  }


  void GReductor::myCopyGBasis_module(VectorList& outG)
  {
    outG.clear();
    if (myGB.empty()) return;
    for (const auto& ptr: myGB)
      if (IsActive(*ptr))
        outG.push_back(DeEmbedPoly((*ptr).myPoly(), myGRingInfoValue));
  }


  void GReductor::myCopyMinGens_module(VectorList& outG)
  {
    outG.clear();
    for (const auto& ptr: myGB)
      if (IsMinimalGen(*ptr))
        outG.push_back(DeEmbedPoly((*ptr).myPoly(), myGRingInfoValue));
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
      if (myTruncDeg() != ourNoTruncValue)
      {
        if (myPairs.empty())  mySetTruncDeg(ourNoTruncValue); // GB complete
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
    if (!IsZZ(CoeffRing(P)))  CoCoA_THROW_ERROR2(ERR::BadRing, "must be over QQ");
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
//     CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
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
        CoCoA_THROW_ERROR2(ERR::BadArg, "incompatible GradingDim");
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
      CoCoA_THROW_ERROR1(ERR::NYI);
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


  SparsePolyRing MakeElimRing(const SparsePolyRing& P,
                              const std::vector<long>& IndetsToElim,
                              const bool IsHomog)
  {
    ConstMatrixView weights = GradingMat(P);

    matrix NewOrdMat = ElimMat(IndetsToElim, weights);
    long NewGrDim = 0;
    if (GradingDim(P) != 0 && IsHomog)
    {
      NewOrdMat = ElimHomogMat(IndetsToElim, weights);
      NewGrDim = GradingDim(P);
    }

    const PPOrdering& NewOrd = NewMatrixOrdering(NewOrdMat, NewGrDim);
    vector<symbol> IndetNames = NewSymbols(NumIndets(P));
    SparsePolyRing NewP=NewPolyRing(CoeffRing(P),IndetNames,NewOrd);
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
