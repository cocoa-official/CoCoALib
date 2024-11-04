//   Copyright (c)  2014-2017  John Abbott and Anna M. Bigatti
//   Authors:  2014-2017  John Abbott, Anna M. Bigatti

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


#include "CoCoA/SparsePolyOps-implicit.H"

#include "CoCoA/DenseMatrix.H"  // for elim (weights) and ImplicitDirectOrd2
#include "CoCoA/MatrixForOrdering.H"  // for elim (weights) and ImplicitDirectOrd2
#include "CoCoA/MatrixOps.H"  // for ImplicitDirectOrd2
#include "CoCoA/MatrixView.H"  // for ImplicitDirectOrd2
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPMonoidHom.H" // for PPMonoidHom
#include "CoCoA/PPOrdering.H"
#include "CoCoA/PPOrdering.H"  // for elim (weights) and ImplicitDirectOrd2
#include "CoCoA/QBGenerator.H"
#include "CoCoA/QuotientRing.H" // for NewZZmod
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"  // for elim (weights)
#include "CoCoA/RingZZ.H"  // for ord matrix (ImplicitDirectOrd2)
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpGOperations.H"  // for ComputeElimFirst
#include "CoCoA/VectorOps.H" // (tmp) for printing
#include "CoCoA/ideal.H" // for NF: checking result
#include "CoCoA/interrupt.H" // for CheckForInterrupt
#include "CoCoA/random.H"
#include "CoCoA/time.H" // (tmp) for timings
#include "CoCoA/verbose.H" // for VerboseLog

#include <map>
using std::map;
#include <string>
using std::string;
//#include <vector>
using std::vector;

namespace CoCoA
{

// ImplicitDirect

  namespace // anonymous for file local fns
  {

    RingElem ComputeImage(const PPMonoidElem& t, const vector<RingElem>& L)
    {
      RingElem ans = one(owner(L[0]));
      const vector<long> exp = exponents(t);
      for (int i=0; i < len(exp); ++i)
        ans *= power(L[i],exp[i]);
      return ans;
    }

    void reduce(RingElem& NewRed, RingElem& NewPolyx, map<PPMonoidElem, int>& FindReducerIndex, const vector<RingElem>& reducer, const vector<RingElem>& polyx)
    {
      //      RingHom embed1 = CoeffEmbeddingHom(owner(NewRed));
      //      RingHom embed2 = CoeffEmbeddingHom(owner(NewPolyx));
      while (true)
      {
        if (IsZero(NewRed)) return;
        const PPMonoidElem PP = LPP(NewRed);
        const int i = -1 + FindReducerIndex[PP]; // remove shift of 1
        if (i == -1) return;
        const RingElem c = LC(NewRed)/LC(reducer[i]);
        //        NewRed -= embed1(c)*reducer[i];
        //        NewPolyx -= embed2(c)*polyx[i];
        RingElem red = reducer[i];
        RingElem pol = polyx[i];
        SparsePolyRingPtr(owner(red))->myMulByCoeff(raw(red), raw(c));
        SparsePolyRingPtr(owner(pol))->myMulByCoeff(raw(pol), raw(c));
        NewRed -= red;
        NewPolyx -= pol;
      }
    }

    // this allows the computation in QQ[...]
    // in fact it is probably better to compute in Fp[...] and lift
    SparsePolyRing NewPolyRingForImplicit(const ring& K, const string& name, long n)
    {
      if (IsQQ(K))
        return NewPolyRing_DMPI(K, SymbolRange(name,1,n));
      return NewPolyRing_DMPII(K, SymbolRange(name,1,n));
      //this should be done properly: check if K is SmallFpImpl
    }

    SparsePolyRing NewPolyRingForImplicit(const ring& K, long n, PPOrdering ord)
    {
      if (IsQQ(K))
        return NewPolyRing_DMPI(K, NewSymbols(n), ord);
      return NewPolyRing_DMPII(K, NewSymbols(n), ord);
      //this should be done properly: check if K is SmallFpImpl
    }

    SparsePolyRing NewPolyRingForImplicit(const ring& K, const string& name, PPOrdering ord)
    {
      if (IsQQ(K))
        return NewPolyRing_DMPI(K, SymbolRange(name,1,NumIndets(ord)), ord);
      return NewPolyRing_DMPII(K, SymbolRange(name,1,NumIndets(ord)), ord);
      //this should be done properly: check if K is SmallFpImpl
    }

  } // end of anonymous namespace ///////////////////////////////////


  RingElem ImplicitDirect(const std::vector<RingElem>& ParamDescrOrig)
  {
    if (ParamDescrOrig.empty() || len(ParamDescrOrig) != 1+NumIndets(owner(ParamDescrOrig[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirectLPP");
    const ring& Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

    ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
    RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
    vector<RingElem> ParamDescr = phi(ParamDescrOrig);
    // Originally the code was like this (ie. compute in the originally given ring)
    // if (ParamDescr.empty() || len(ParamDescr) != 1+NumIndets(owner(ParamDescr[0])))
    //   CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirect");
    // ring P = owner(ParamDescr[0]);
    // const int n = len(ParamDescr);
    ring Kx = NewPolyRingForImplicit(CoeffRing(Porig), "x", n);
    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;

    reducer.push_back(one(Kt));
    polyx.push_back(one(Kx));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    QBGenerator QBG(PPM(Kx));
    PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    while (true)
    {
      const PPMonoidElem PP = QBG.myCorners().front();
      QBG.myCornerPPIntoQB(PP);
      RingElem NewReducer = ComputeImage(PP, ParamDescr);
      RingElem NewPolyx = monomial(Kx, PP);
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
    }
  }


  PPMonoidElem ComputeLPP(const PPMonoidElem& t, const vector<PPMonoidElem>& VecLPP)
  {
    const PPMonoid& PPM = owner(VecLPP[0]);
    PPMonoidElem ans = one(PPM);
    const int n = len(VecLPP);
    CoCoA_ASSERT(n == NumIndets(owner(t)));
    const vector<long> exp = exponents(t);
    for (int i=0;i<n;++i)
      ans *= power(VecLPP[i],exp[i]);
    return ans;
  }


  RingElem ImplicitDirectLPP(const std::vector<RingElem>& ParamDescrOrig)
  {
    if (ParamDescrOrig.empty() || len(ParamDescrOrig) != 1+NumIndets(owner(ParamDescrOrig[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirectLPP");
    const ring& Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

    ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
    RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
    vector<RingElem> ParamDescr = phi(ParamDescrOrig);
    ring Kx = NewPolyRingForImplicit(CoeffRing(Porig), "x",n);
    vector<PPMonoidElem> VecLPP; VecLPP.reserve(len(ParamDescr)); for(int i=0;i<len(ParamDescr);++i)VecLPP.push_back(LPP(ParamDescr[i]));
    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;

    polyx.push_back(one(Kx));
    reducer.push_back(one(Kt));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    QBGenerator QBG(PPM(Kx));
    PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    while (true)
    {
      using std::list;
      // pick corner elem giving smallest LPP in image
//      list<PPMonoidElem>::const_iterator BestIter = QBG.myCorners().begin();
      PPMonoidElem BestCorner = QBG.myCorners().front();
//      PPMonoidElem BestLPP = ComputeLPP(*BestIter, VecLPP);
      PPMonoidElem BestLPP = ComputeLPP(BestCorner, VecLPP);

//      for(list<PPMonoidElem>::const_iterator it=QBG.myCorners().begin(); it != QBG.myCorners().end(); ++it)
        for (const PPMonoidElem& t: QBG.myCorners())
      {
        const PPMonoidElem tmp = ComputeLPP(t, VecLPP);
          if (tmp < BestLPP) { BestLPP = tmp; BestCorner = t; }
      }
//      const PPMonoidElem PP = *BestIter;
      QBG.myCornerPPIntoQB(BestCorner);
      RingElem NewReducer = ComputeImage(BestCorner,ParamDescr);
      RingElem NewPolyx = monomial(Kx, BestCorner);
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
    }
  }



  int ChooseIndet(const PPMonoidElem& t, const vector<RingElem>& ParamDescr, std::map<PPMonoidElem, int>& FindPP, const std::vector<RingElem>& polyx)
    {
      const vector<long> expv = exponents(t);
      const int n = len(expv); // NumIndets(PPM)
      PPMonoid PPMt = owner(t);
      int BestI = -1;
      long BestCost = 0; // should never be used
      for (int i=0; i<n; ++i)
      {
        if (expv[i] == 0) continue;
        //        const PPMonoidElem& x_i = indet(PPMt,i);
        //        const PPMonoidElem tfactor = t/x_i;
        const PPMonoidElem tfactor = t/indet(PPMt,i);
        const int j = FindPP[tfactor];
        const long EstCost = NumTerms(ParamDescr[i])*NumTerms(polyx[j]);
        if (BestCost == 0 || BestCost > EstCost)
        {
          BestCost = EstCost;
          BestI = i;
        }
      }
      return BestI;
    }

//   void ComputeImage(RingElem& NewReducer, RingElem& NewPolyx, const PPMonoidElem& t, const vector<RingElem>& L, const std::map<PPMonoidElem, int>& FindFactorIndex, const std::vector<RingElem>& polyx)
//     {
// ///      RingElem ans = one(owner(L[0]));
//       const vector<long> exp = exponents(t);
//       const int n = len(exp); // NumIndets(PPM)
//       PPMonoid PPM = owner(t);
//     int BestI = -1;
//     long BestCost = 0; // should never be used
//     for (int i=0; i<n; ++i)
//     {
//       if (exp[i] == 0) continue;
//       const int j = FindFactorIndex(t/indet(PPM,i));
//       const long EstCost = NumTerms(L[i])*NumTerms(polyx[j]);
//       if (BestI >=0 && EstCost > BestCost) continue;
//       BestCost = EstCost;
//       BestI = i;
//       BestFactor = t/indet(PPM,i);
//       BestPoly = j;
//     }
//     const ring& P = owner(L[0]);
//     NewPolyx = L[BestI]*polyx[BestPoly];
//     NewReducer = reducer[FindReducerIndex[BestFactor]]*indet(PP,BestI);
//     }

  // PPMonoidElem ComputeLPP(const PPMonoidElem& t, const vector<PPMonoidElem>& VecLPP)
  // {
  //   const PPMonoid& PPM = owner(VecLPP[0]);
  //   PPMonoidElem ans = one(PPM);
  //   const int n = len(VecLPP);
  //   CoCoA_ASSERT(n == NumIndets(owner(t)));
  //   const vector<long> exp = exponents(t);
  //   for (int i=0;i<n;++i)
  //     ans *= power(VecLPP[i],exp[i]);
  //   return ans;
  // }

  RingElem ImplicitDirectLPP2Homog(const std::vector<RingElem>& ParamDescrOrig);
  RingElem ImplicitDirectLPP2HomogBis(const std::vector<RingElem>& ParamDescrOrig);
  
  RingElem ImplicitDirectLPP2(const std::vector<RingElem>& ParamDescrOrig)
  {
    if (ParamDescrOrig.empty() || len(ParamDescrOrig) != 1+NumIndets(owner(ParamDescrOrig[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirectLPP");
    if (IsHomog(ParamDescrOrig))
    {
      long i;
      long d = deg(ParamDescrOrig[0]);
      for (i=1; i<len(ParamDescrOrig); ++i)
        if (d != deg(ParamDescrOrig[i])) break;
      if (i==len(ParamDescrOrig)) return ImplicitDirectLPP2HomogBis(ParamDescrOrig);
    }
    const ring& Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

    ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
    RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
    vector<RingElem> ParamDescr = phi(ParamDescrOrig);
//    ring Kx = NewPolyRing(CoeffRing(P), SymbolRange("x",1,n));
    ring Kx = NewPolyRingForImplicit(CoeffRing(Porig), "x",n);
    vector<PPMonoidElem> VecLPP; VecLPP.reserve(len(ParamDescr)); for(int i=0;i<len(ParamDescr);++i)VecLPP.push_back(LPP(ParamDescr[i]));
    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;
    map<PPMonoidElem,int> FindFactorIndex;

    polyx.push_back(one(Kx));
    FindFactorIndex[LPP(polyx[0])] = 0;
    reducer.push_back(one(Kt));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    QBGenerator QBG(PPM(Kx));
    PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    while (true)
    {
      using std::list;
      // pick corner elem giving smallest LPP in image
      list<PPMonoidElem>::const_iterator BestIter = QBG.myCorners().begin();
      PPMonoidElem BestLPP = ComputeLPP(*BestIter, VecLPP);

      for(list<PPMonoidElem>::const_iterator it=QBG.myCorners().begin(); it != QBG.myCorners().end(); ++it)
//      for(int i=1;i<len(QBG.myCorners());++i)
      {
        const PPMonoidElem tmp = ComputeLPP(*it, VecLPP);
          if (tmp < BestLPP) { BestLPP = tmp; BestIter = it; }
      }
      const PPMonoidElem PP = *BestIter;
      QBG.myCornerPPIntoQB(PP);
      const int var = ChooseIndet(PP, ParamDescr, FindFactorIndex, polyx);
      const int ReducerIndex = FindFactorIndex[PP/indet(PPM(Kx),var)];
      RingElem NewPolyx = indet(Kx,var)*polyx[ReducerIndex];
      RingElem NewReducer = reducer[ReducerIndex]*ParamDescr[var];
//////      ComputeImage(NewReducer, NewPolyx, PP,ParamDescr, FindFactorIndex, polyx);
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
      FindFactorIndex[PP] = len(reducer)-1;
    }
  }

  RingElem ImplicitDirectOrd2HomogBis(const std::vector<RingElem>& ParamDescrOrig);

  RingElem ImplicitDirectOrd2(const std::vector<RingElem>& ParamDescrOrig)
  {
    if (ParamDescrOrig.empty() || len(ParamDescrOrig) != 1+NumIndets(owner(ParamDescrOrig[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirectOrd2");
    if (IsHomog(ParamDescrOrig))
    {
      long i;
      long d = deg(ParamDescrOrig[0]);
      for (i=1; i<len(ParamDescrOrig); ++i)
        if (d != deg(ParamDescrOrig[i])) break;
      if (i==len(ParamDescrOrig)) return ImplicitDirectOrd2HomogBis(ParamDescrOrig);
    }
    const ring& Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

    ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
    RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
    vector<RingElem> ParamDescr = phi(ParamDescrOrig);

    //  ordering compatible with LPP in Kt
    vector< vector<long> >  VV(n, vector<long>(n));
    for (long i=0; i<n; ++i)  exponents(VV[i], LPP(ParamDescr[i]));
    
    matrix M = MakeTermOrdMat(StdDegRevLexMat(n-1) * transpose(NewDenseMat(RingZZ(),VV)));
    //    std::cout << M << std::endl;
 
    /// we should always check that there are no constants in ParamDescr
    ring Kx = NewPolyRingForImplicit(CoeffRing(Porig), "x", NewMatrixOrdering(M, 1));

    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;
    map<PPMonoidElem,int> FindFactorIndex;

    polyx.push_back(one(Kx));
    FindFactorIndex[LPP(polyx[0])] = 0;
    reducer.push_back(one(Kt));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    QBGenerator QBG(PPM(Kx));
    PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    while (true)
    {
      using std::list;
      // pick corner elem giving smallest LPP in image
      const PPMonoidElem PP = QBG.myCorners().front(); // special ord
      QBG.myCornerPPIntoQB(PP);
      const int var = ChooseIndet(PP, ParamDescr, FindFactorIndex, polyx);
      const int ReducerIndex = FindFactorIndex[PP/indet(PPM(Kx),var)];
      RingElem NewPolyx = indet(Kx,var)*polyx[ReducerIndex];
      RingElem NewReducer = reducer[ReducerIndex]*ParamDescr[var];
////  ComputeImage(NewReducer, NewPolyx, PP,ParamDescr, FindFactorIndex, polyx);
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
      FindFactorIndex[PP] = len(reducer)-1;
    }
  }


  RingElem ImplicitDirectLPP2Homog(const std::vector<RingElem>& ParamDescrOrig)
  {
    using std::list;
    VerboseLog VERBOSE("ImplicitDirectLPP2Homog");
    VERBOSE(90) << "start" << std::endl;
    const ring& Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

    ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
    RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
    vector<RingElem> ParamDescr = phi(ParamDescrOrig);
    ring Kx = NewPolyRingForImplicit(CoeffRing(Porig), "x",n);
    vector<PPMonoidElem> VecLPP; VecLPP.reserve(len(ParamDescr)); for(int i=0;i<len(ParamDescr);++i)VecLPP.push_back(LPP(ParamDescr[i]));
    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;
    map<PPMonoidElem,int> FindFactorIndex;

    reducer.push_back(one(Kt));
    polyx.push_back(one(Kx));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    FindFactorIndex[LPP(polyx[0])] = 0;
    QBGenerator QBG(PPM(Kx));
    PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    while (true)
    {
      list<PPMonoidElem>::const_iterator BestIter = QBG.myCorners().begin();
      PPMonoidElem BestLPP = *BestIter;
      for(list<PPMonoidElem>::const_iterator it=QBG.myCorners().begin(); it != QBG.myCorners().end(); ++it)
      {
        const PPMonoidElem tmp = *it;
          if (tmp < BestLPP) { BestLPP = tmp; BestIter = it; }
      }
      const PPMonoidElem PPx = *BestIter;
//       if (deg(PPx) != CurrentDeg)
//       {
//         CurrentDeg = deg(PPx);
//         //        std::cout << "deg(PPx) = " << deg(PPx) << std::endl;
//       }
      QBG.myCornerPPIntoQB(PPx);
      const int var = ChooseIndet(PPx, ParamDescr, FindFactorIndex, polyx);
      const int ReducerIndex = FindFactorIndex[PPx/indet(PPM(Kx),var)];
      RingElem NewPolyx   = polyx[ReducerIndex]  *indet(Kx,var);
      RingElem NewReducer = reducer[ReducerIndex]*ParamDescr[var];
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
      FindFactorIndex[PPx] = len(reducer)-1;
    }
  }


  RingElem ImplicitDirectLPP2HomogBis(const std::vector<RingElem>& ParamDescrOrig)
  {
    using std::list;
    VerboseLog VERBOSE("ImplicitDirectLPP2HomogBis");
    VERBOSE(90) << "start" << std::endl;
    const ring& Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

    ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
    RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
    vector<RingElem> ParamDescr = phi( ParamDescrOrig);
    ring Kx = NewPolyRingForImplicit(CoeffRing(Porig), "x",n);
    vector<PPMonoidElem> VecLPP; VecLPP.reserve(len(ParamDescr)); for(int i=0;i<len(ParamDescr);++i)VecLPP.push_back(LPP(ParamDescr[i]));
    vector<RingElem> reducerPrev;
    vector<RingElem> polyxPrev;
    map<PPMonoidElem,int> FindReducerIndexPrev;
    map<PPMonoidElem,int> FindFactorIndexPrev;
    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;
    map<PPMonoidElem,int> FindFactorIndex;

    reducer.push_back(one(Kt));
    polyx.push_back(one(Kx));
    reducerPrev.push_back(one(Kt));
    polyxPrev.push_back(one(Kx));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    FindFactorIndex[LPP(polyx[0])] = 0;
    FindReducerIndexPrev[LPP(reducer[0])] = 1; // deliberately add 1
    FindFactorIndexPrev[LPP(polyx[0])] = 0;
    QBGenerator QBG(PPM(Kx));
    PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    long CurrentDeg = 0;
    while (true)
    {
      list<PPMonoidElem>::const_iterator BestIter = QBG.myCorners().begin();
      PPMonoidElem BestLPP = *BestIter;
      for(list<PPMonoidElem>::const_iterator it=QBG.myCorners().begin(); it != QBG.myCorners().end(); ++it)
      {
        const PPMonoidElem tmp = *it;
          if (tmp < BestLPP) { BestLPP = tmp; BestIter = it; }
      }
      const PPMonoidElem PPx = *BestIter;
      if (deg(PPx) != CurrentDeg)
      {
        CurrentDeg = deg(PPx);
        reducerPrev.clear();
        polyxPrev.clear();
        FindReducerIndexPrev.clear();
        FindFactorIndexPrev.clear();
        reducerPrev.swap(reducer);
        polyxPrev.swap(polyx);
        FindReducerIndexPrev.swap(FindReducerIndex);
        FindFactorIndexPrev.swap(FindFactorIndex);
        reducer.push_back(one(Kt));
        polyx.push_back(one(Kx));
        FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
        FindFactorIndex[LPP(polyx[0])] = 0;
      }
      QBG.myCornerPPIntoQB(PPx);
      const int var = ChooseIndet(PPx, ParamDescr, FindFactorIndexPrev, polyxPrev);
      const int ReducerIndex = FindFactorIndexPrev[PPx/indet(PPM(Kx),var)];
      RingElem NewPolyx   = polyxPrev[ReducerIndex]  *indet(Kx,var);
      RingElem NewReducer = reducerPrev[ReducerIndex]*ParamDescr[var];
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
      FindFactorIndex[PPx] = len(reducer)-1;
    }
  }


  RingElem ImplicitDirectOrd2HomogBis(const std::vector<RingElem>& ParamDescrOrig)
  {
    using std::list;
    VerboseLog VERBOSE("ImplicitDirectOrd2HomogBis");
    VERBOSE(90) << "start" << std::endl;

    const ring& Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

    ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
    RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
    vector<RingElem> ParamDescr = phi(ParamDescrOrig);
    //    ring Kx = NewPolyRingForImplicit(CoeffRing(Porig), "x",n);
    vector<PPMonoidElem> VecLPP; VecLPP.reserve(len(ParamDescr)); for(int i=0;i<len(ParamDescr);++i)VecLPP.push_back(LPP(ParamDescr[i]));
    //  ordering compatible with LPP in Kt
    vector< vector<long> > VV(n, vector<long>(n));
    for (long i=0; i<n; ++i)  exponents(VV[i], LPP(ParamDescr[i]));    
    matrix M = MakeTermOrdMat(StdDegRevLexMat(n-1) * transpose(NewDenseMat(RingZZ(),VV)));
    //    std::cout << M << std::endl;
 
    /// we should always check that there are no constants in ParamDescr
    ring Kx = NewPolyRingForImplicit(CoeffRing(Porig), "x", NewMatrixOrdering(M, 1));

    vector<RingElem> reducerPrev;
    vector<RingElem> polyxPrev;
    map<PPMonoidElem,int> FindReducerIndexPrev;
    map<PPMonoidElem,int> FindFactorIndexPrev;
    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;
    map<PPMonoidElem,int> FindFactorIndex;

    reducer.push_back(one(Kt));
    polyx.push_back(one(Kx));
    reducerPrev.push_back(one(Kt));
    polyxPrev.push_back(one(Kx));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    FindFactorIndex[LPP(polyx[0])] = 0;
    FindReducerIndexPrev[LPP(reducer[0])] = 1; // deliberately add 1
    FindFactorIndexPrev[LPP(polyx[0])] = 0;
    QBGenerator QBG(PPM(Kx));
    PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    degree CurrentDeg = wdeg(PP1);
    while (true)
    {
//       list<PPMonoidElem>::const_iterator BestIter = QBG.myCorners().begin();
//       PPMonoidElem BestLPP = *BestIter;
//       for(list<PPMonoidElem>::const_iterator it=QBG.myCorners().begin(); it != QBG.myCorners().end(); ++it)
//       {
//         const PPMonoidElem tmp = *it;
//           if (tmp < BestLPP) { BestLPP = tmp; BestIter = it; }
//       }
// //       const PPMonoidElem PPx = *BestIter;
//       if (*BestIter != QBG.myCorners().front())
//       {
//         std::cout  << "ppx is " << QBG.myCorners().front() <<
//                    << "  *BestIter is " << *BestIter << std::endl;
//       }
      const PPMonoidElem PPx = QBG.myCorners().front(); // special ord
      if (wdeg(PPx) != CurrentDeg)
      {
        CurrentDeg = wdeg(PPx);
        reducerPrev.clear();
        polyxPrev.clear();
        FindReducerIndexPrev.clear();
        FindFactorIndexPrev.clear();
        reducerPrev.swap(reducer);
        polyxPrev.swap(polyx);
        FindReducerIndexPrev.swap(FindReducerIndex);
        FindFactorIndexPrev.swap(FindFactorIndex);
        reducer.push_back(one(Kt));
        polyx.push_back(one(Kx));
        FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
        FindFactorIndex[LPP(polyx[0])] = 0;
      }
      QBG.myCornerPPIntoQB(PPx);
      const int var = ChooseIndet(PPx, ParamDescr, FindFactorIndexPrev, polyxPrev);
      const int ReducerIndex = FindFactorIndexPrev[PPx/indet(PPM(Kx),var)];
      RingElem NewPolyx   = polyxPrev[ReducerIndex]  *indet(Kx,var);
      RingElem NewReducer = reducerPrev[ReducerIndex]*ParamDescr[var];
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
      FindFactorIndex[PPx] = len(reducer)-1;
    }
  }




  //  Procedure for curves.  but there is no stopping criterion :-(
//   std::vector<RingElem> BM_param(const std::vector<RingElem>& ParamDescrOrig)
//   {
//     if (ParamDescrOrig.empty())
//       CoCoA_THROW_ERROR(ERR::BadArg, "BM_param");
//     const ring& Porig = owner(ParamDescrOrig[0]);
//     const int n = len(ParamDescrOrig);

//     ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
//     RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
//     vector<RingElem> ParamDescr = phi(ParamDescrOrig);
//     ring Kx = NewPolyRingForImplicit(CoeffRing(Porig), "x",n);
//     vector<PPMonoidElem> VecLPP; VecLPP.reserve(len(ParamDescr)); 
//     for (int i=0;i<n;++i) VecLPP.push_back(LPP(ParamDescr[i]));
//     vector<RingElem> reducer;
//     vector<RingElem> polyx;
//     map<PPMonoidElem,int> FindReducerIndex;
//     map<PPMonoidElem,int> FindFactorIndex;

//     reducer.push_back(one(Kt));
//     polyx.push_back(one(Kx));
//     FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
//     FindFactorIndex[LPP(polyx[0])] = 0;
//     QBGenerator QBG(PPM(Kx));
//     PPMonoidElem PP1 = QBG.myCorners().front();
//     QBG.myCornerPPIntoQB(PP1);
//     vector<RingElem> GB(0);
//     while (!QBG.myCorners().empty())
//     {
//       std::cout << "QBG = " << QBG << std::endl;
      
//       using std::list;
//       const PPMonoidElem t = QBG.myCorners().front();
//       const int var = ChooseIndet(t, ParamDescr, FindFactorIndex, polyx);
//       const int ReducerIndex = FindFactorIndex[t/indet(PPM(Kx),var)];
//       RingElem NewPolyx = indet(Kx,var)*polyx[ReducerIndex];
//       RingElem NewReducer = reducer[ReducerIndex]*ParamDescr[var];
//       reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
//       //      if (IsZero(NewReducer)) return NewPolyx;
//       if (IsZero(NewReducer))
//       {
//         QBG.myCornerPPIntoAvoidSet(t);
//         GB.push_back(NewPolyx);
//         std::cout << "--GB = " << GB << std::endl;
//       }
//       else
//       {
//         QBG.myCornerPPIntoQB(t);
//         reducer.push_back(NewReducer);
//         polyx.push_back(NewPolyx);
//         FindReducerIndex[LPP(NewReducer)] = len(reducer);
//         FindFactorIndex[t] = len(reducer)-1;
//       }
//     }
//     return GB;
//   }


  RingElem ImplicitDirectWithCond(const std::vector<RingElem>& ParamDescrOrig, const std::vector<RingElem>& cond)
  {
    if (ParamDescrOrig.empty() || len(ParamDescrOrig)+len(cond) != 1+NumIndets(owner(ParamDescrOrig[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirectWithCond");
    if (cond.empty()) return ImplicitDirect(ParamDescrOrig);
    const ring& Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

//     ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t",NumIndets(Porig));
//     RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
//     vector<RingElem> ParamDescr = phi(ParamDescrOrig);
//     vector<RingElem> relations = phi(cond);
    ring Kt = Porig;
    vector<RingElem> ParamDescr = ParamDescrOrig;
    vector<RingElem> relations = cond;
    ideal RelationIdeal(relations);
    // Originally the code was like this (ie. compute in the originally given ring)
    // if (ParamDescr.empty() || len(ParamDescr) != 1+NumIndets(owner(ParamDescr[0])))
    //   CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirect");
    // const ring& P = owner(ParamDescr[0]);
    // const int n = len(ParamDescr);
    ring Kx = NewPolyRingForImplicit(CoeffRing(Kt), "x",n);
    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;

    reducer.push_back(one(Kt));
    polyx.push_back(one(Kx));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    QBGenerator QBG(PPM(Kx));
//     PPMonoidElem PP1 = QBG.myCorners().front();
//     QBG.myCornerPPIntoQB(PP1);
    QBG.myCornerPPIntoQB(QBG.myCorners().front());
    while (true)
    {
      const PPMonoidElem t = QBG.myCorners().front();
      QBG.myCornerPPIntoQB(t);
      RingElem NewReducer = NF(ComputeImage(t,ParamDescr), RelationIdeal);
      RingElem NewPolyx = monomial(Kx, t);
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
    }
  }


  RingElem ImplicitDirectWithCondLPP(const std::vector<RingElem>& ParamDescrOrig, const std::vector<RingElem>& cond)
  {
    if (ParamDescrOrig.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (len(ParamDescrOrig)+len(cond) != 1+NumIndets(owner(ParamDescrOrig[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirectWithCondLPP");
    if (cond.empty()) return ImplicitDirectLPP(ParamDescrOrig);
    ring Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

    ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
    RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
    vector<RingElem> ParamDescr = phi(ParamDescrOrig);
    vector<RingElem> relations = phi(cond);
    ideal RelationIdeal(relations);
    // Originally the code was like this (ie. compute in the originally given ring)
    // if (ParamDescr.empty() || len(ParamDescr) != 1+NumIndets(owner(ParamDescr[0])))
    //   CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirect");
    // const ring& P = owner(ParamDescr[0]);
    // const int n = len(ParamDescr);
    ring Kx = NewPolyRingForImplicit(CoeffRing(Kt), "x",n);
    vector<PPMonoidElem> VecLPP; VecLPP.reserve(len(ParamDescr)); for(int i=0;i<len(ParamDescr);++i)VecLPP.push_back(LPP(ParamDescr[i]));
    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;

    reducer.push_back(one(Kt));
    polyx.push_back(one(Kx));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    QBGenerator QBG(PPM(Kx));
    PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    while (true)
    {
      using std::list;
      // pick corner elem giving smallest LPP in image
      list<PPMonoidElem>::const_iterator BestIter = QBG.myCorners().begin();
      PPMonoidElem BestLPP = ComputeLPP(*BestIter, VecLPP);

      for(list<PPMonoidElem>::const_iterator it=QBG.myCorners().begin(); it != QBG.myCorners().end(); ++it)
//      for(int i=1;i<len(QBG.myCorners());++i)
      {
        const PPMonoidElem tmp = ComputeLPP(*it, VecLPP);
          if (tmp < BestLPP) { BestLPP = tmp; BestIter = it; }
      }
      const PPMonoidElem t = *BestIter;
      QBG.myCornerPPIntoQB(t);
      RingElem NewReducer = NF(ComputeImage(t,ParamDescr), RelationIdeal);
      RingElem NewPolyx = monomial(Kx, t);
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
    }
  }


  RingElem ImplicitDirectWithCondOrd2(const std::vector<RingElem>& ParamDescrOrig, const std::vector<RingElem>& cond)
  {
    if (ParamDescrOrig.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (len(ParamDescrOrig)+len(cond) != 1+NumIndets(owner(ParamDescrOrig[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirectWithCondOrd2");
    if (cond.empty()) return ImplicitDirectOrd2(ParamDescrOrig);
    CoCoA_THROW_ERROR(ERR::NYI, "ImplicitDirectOrd2");
    return zero(RingZZ()); // just to keep the compiler quiet
  }
  

  RingElem ImplicitDirectWithCondLPP2(const std::vector<RingElem>& ParamDescrOrig, const std::vector<RingElem>& cond)
  {
    if (ParamDescrOrig.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (len(ParamDescrOrig)+len(cond) != 1+NumIndets(owner(ParamDescrOrig[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitDirectWithCondLPP2");
    if (cond.empty()) return ImplicitDirectLPP2(ParamDescrOrig);
    const ring& Porig = owner(ParamDescrOrig[0]);
    const int n = len(ParamDescrOrig);

    // faster if creating a new ring (for locality?!?)
    ring Kt = NewPolyRingForImplicit(CoeffRing(Porig), "t", NumIndets(Porig));
    RingHom phi = PolyAlgebraHom(Porig, Kt, indets(Kt));
    vector<RingElem> ParamDescr = phi(ParamDescrOrig);
    vector<RingElem> relations = phi(cond);
    // ring Kt = Porig;
    // const vector<RingElem>& ParamDescr = ParamDescrOrig;
    // const vector<RingElem>& relations = cond;
    ideal RelationIdeal(relations);
    ring Kx = NewPolyRingForImplicit(CoeffRing(Kt), "x", n);
    vector<PPMonoidElem> VecLPP;
    VecLPP.reserve(len(ParamDescr));
    for(int i=0;i<len(ParamDescr);++i) VecLPP.push_back(LPP(ParamDescr[i]));
    vector<RingElem> reducer;
    vector<RingElem> polyx;
    map<PPMonoidElem,int> FindReducerIndex;
    map<PPMonoidElem,int> FindFactorIndex; // -- LPP2

    polyx.push_back(one(Kx));
    reducer.push_back(one(Kt));
    FindReducerIndex[LPP(reducer[0])] = 1; // deliberately add 1
    QBGenerator QBG(PPM(Kx));
    PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    while (true)
    {
      using std::list;
      // pick corner elem giving smallest LPP in image
      list<PPMonoidElem>::const_iterator BestIter = QBG.myCorners().begin();
      PPMonoidElem BestLPP = ComputeLPP(*BestIter, VecLPP);

      for(list<PPMonoidElem>::const_iterator it=QBG.myCorners().begin(); it != QBG.myCorners().end(); ++it)
//      for(int i=1;i<len(QBG.myCorners());++i)
      {
        const PPMonoidElem tmp = ComputeLPP(*it, VecLPP);
          if (tmp < BestLPP) { BestLPP = tmp; BestIter = it; }
      }
      const PPMonoidElem PP = *BestIter;
      QBG.myCornerPPIntoQB(PP);
      const int var = ChooseIndet(PP, ParamDescr, FindFactorIndex, reducer);
      const int ReducerIndex = FindFactorIndex[PP/indet(PPM(Kx),var)];
      RingElem NewPolyx = polyx[ReducerIndex]*indet(Kx,var);
      RingElem NewReducer = NF(reducer[ReducerIndex]*ParamDescr[var], RelationIdeal);
      reduce(NewReducer, NewPolyx, FindReducerIndex, reducer, polyx);
      if (IsZero(NewReducer)) return NewPolyx;
      reducer.push_back(NewReducer);
      polyx.push_back(NewPolyx);
      FindReducerIndex[LPP(NewReducer)] = len(reducer);
      FindFactorIndex[PP] = len(reducer)-1;
    }
  }



  //-------------------------------------------------------
  // ImplicitByPoints

  namespace // anonymous for file local fns
  {

    RingElem eval(const RingElem& f, const vector<RingElem>& pt)
    {
      const ring& P = owner(f);
      RingElem ans = zero(CoeffRing(P));
      RingHom phi = CoeffEmbeddingHom(P);
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
      {
        RingElem ThisTerm = coeff(it);
        for (int i=0; i < NumIndets(P); ++i)
          ThisTerm *= power(pt[i], exponent(PP(it),i));
        ans += ThisTerm;
      }
      return ans;
    }

    // currently just use random points
    vector<RingElem> GetNewSamplePt(const vector<RingElem>& ParamDescr, const vector< vector<RingElem> >& SamplePts)
    {
      const ring& R = owner(ParamDescr[0]);
      const int NumParams = NumIndets(R);
      vector<RingElem> params(NumParams);
      vector<RingElem> NewSamplePt(len(ParamDescr));
      TryAgain:
      for (int i=0; i < NumParams; ++i)
      {
        params[i] = RingElem(CoeffRing(R),RandomLong(-999,999));
      }
      for (int i=0; i < len(ParamDescr); ++i)
        NewSamplePt[i] = eval(ParamDescr[i], params);
      if (find(SamplePts.begin(), SamplePts.end(), NewSamplePt) != SamplePts.end()) goto TryAgain;
      return NewSamplePt;
    }

  } // end of anonymous namespace

  RingElem ImplicitByPoints(const std::vector<RingElem>& ParamDescr)
  {
    if (ParamDescr.empty() || len(ParamDescr) != 1+NumIndets(owner(ParamDescr[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitByPoints");
    const ring& P = owner(ParamDescr[0]);
    const long n = len(ParamDescr);
    ring Kx = NewPolyRing(CoeffRing(P), SymbolRange("x", 1, n));
    RingHom phi = CoeffEmbeddingHom(Kx);

    QBGenerator QBG(PPM(Kx));
    const PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);

    const RingElem one = monomial(Kx, PP1);
    vector<RingElem> polyx; polyx.push_back(one);
    vector< vector<RingElem> > RowReducers(1); // essentially the BM matrix

    vector< vector<RingElem> > SamplePts;  // should be a std::set!!!
    int NumPts = 0;
    while (true)
    {
      ++NumPts;
      const vector<RingElem> NewSamplePt = GetNewSamplePt(ParamDescr, SamplePts);
      SamplePts.push_back(NewSamplePt);

      for (long i=0; i < NumPts; ++i)
      {
        RowReducers[i].push_back(eval(polyx[i], NewSamplePt)); // clear but inefficient
      }
      if (IsZero(RowReducers[NumPts-1][NumPts-1])) return polyx[NumPts-1];
      const PPMonoidElem t = QBG.myCorners().front();
      QBG.myCornerPPIntoQB(t);
      RingElem NewPolyx = monomial(Kx, t);
      vector<RingElem> NewRow;
      for(long i=0; i < NumPts ; ++i)
        NewRow.push_back(eval(NewPolyx, SamplePts[i]));

      // just gaussian reduction -- we know that matrix is upper triangular
      for (long i=0; i < NumPts; ++i)
      {
        if (IsZero(NewRow[i])) continue;
        const RingElem c = NewRow[i]/RowReducers[i][i];
        NewPolyx -= phi(c)*polyx[i];
        for (long j=i; j < NumPts; ++j)
          NewRow[j] -= c*RowReducers[i][j];
      }
      polyx.push_back(NewPolyx);
      RowReducers.push_back(NewRow);
    }
  }



  RingElem ImplicitByPoints2(const std::vector<RingElem>& ParamDescr)
  {
    using std::list;
    if (ParamDescr.empty() || len(ParamDescr) != 1+NumIndets(owner(ParamDescr[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitByPoints");
    const ring& P = owner(ParamDescr[0]);
    const long p = ConvertTo<long>(characteristic(P));
    CoCoA_ASSERT(IsPrime(p));
    SmallFpImpl ModP(p);
    const long n = len(ParamDescr);
    ring Kx = NewPolyRing(CoeffRing(P), SymbolRange("x", 1, n));
    RingHom phi = CoeffEmbeddingHom(Kx);
    vector<PPMonoidElem> VecLPP; VecLPP.reserve(len(ParamDescr)); for(int i=0;i<len(ParamDescr);++i)VecLPP.push_back(LPP(ParamDescr[i]));

    QBGenerator QBG(PPM(Kx));
    const PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);

    const RingElem one = monomial(Kx, PP1);
    vector<RingElem> polyx; polyx.push_back(one);
///JAA    vector< vector<SmallFpImpl::value_t> > RowReducers(1); // essentially the BM matrix
    vector< vector<SmallFpImpl::value> > RowReducers2(1); // essentially the BM matrix

    vector< vector<RingElem> > SamplePts;  // should be a std::set!!!
    int NumPts = 0;

    while (true)
    {
      ++NumPts;
      const vector<RingElem> NewSamplePt = GetNewSamplePt(ParamDescr, SamplePts);
      SamplePts.push_back(NewSamplePt);

      for (long i=0; i < NumPts; ++i)
      {
///JAA        RowReducers[i].push_back(ModP.myReduce(ConvertTo<long>(eval(polyx[i], NewSamplePt)))); // clear but inefficient
        RowReducers2[i].push_back(ModP.myReduce(ConvertTo<long>(eval(polyx[i], NewSamplePt)))); // clear but inefficient
      }
///JAA      CoCoA_ASSERT(ModP.myExport(RowReducers[NumPts-1][NumPts-1]) == ModP.myExport(RowReducers2[NumPts-1][NumPts-1]));
      if (IsZero(RowReducers2[NumPts-1][NumPts-1])) return polyx[NumPts-1];

      // pick corner elem giving smallest LPP in image
      list<PPMonoidElem>::const_iterator BestIter = QBG.myCorners().begin();
      PPMonoidElem BestLPP = ComputeLPP(*BestIter, VecLPP);

      for(list<PPMonoidElem>::const_iterator it=QBG.myCorners().begin(); it != QBG.myCorners().end(); ++it)
//      for(int i=1;i<len(QBG.myCorners());++i)
      {
        const PPMonoidElem tmp = ComputeLPP(*it, VecLPP);
          if (tmp < BestLPP) { BestLPP = tmp; BestIter = it; }
      }
      const PPMonoidElem t = *BestIter;
      QBG.myCornerPPIntoQB(t);
      RingElem NewPolyx = monomial(Kx, t);
///JAA      vector<SmallFpImpl::value_t> NewRow;
      vector<SmallFpImpl::NonRedValue> NewRow2;
      for(long i=0; i < NumPts ; ++i)
      {
///JAA        NewRow.push_back(ModP.myReduce(ConvertTo<long>(eval(NewPolyx, SamplePts[i]))));
        NewRow2.push_back(ModP.myReduce(ConvertTo<long>(eval(NewPolyx, SamplePts[i]))));
///JAA        CoCoA_ASSERT(ModP.myExport(NewRow[i]) == ModP.myExport(ModP.myNormalize(NewRow2[i])));
      }

      // just gaussian reduction -- we know that matrix is upper triangular
      long IterCount = 0;
      for (long i=0; i < NumPts; ++i)
      {
///JAA        NewRow[i] = ModP.myNormalize(NewRow[i]);
        const SmallFpImpl::value coordi = ModP.myNormalize(NewRow2[i]);
///JAA        NewRow2[i] = coordi;  /// USELESS???
        if (IsZero(coordi)) continue;
///JAA        const SmallFpImpl::value_t c = ModP.myDiv(NewRow[i], RowReducers[i][i]);
        const SmallFpImpl::value c2 = ModP.myNegate(ModP.myDiv(coordi, RowReducers2[i][i]));
///JAA        CoCoA_ASSERT(ModP.myExport(ModP.myNegate(c)) == ModP.myExport(c2));
///JAA        NewRow[i] = 0;
        NewRow2[i] = zero(SmallFp); // USELESS????
///JAA        NewPolyx -= RingElem(Kx,c)*polyx[i];
        NewPolyx += ModP.myExportNonNeg(c2)*polyx[i];
        for (long j=i+1; j < NumPts; ++j)
        {
///JAA          NewRow[j] += (p-c)*RowReducers[i][j];
          NewRow2[j] += c2 * RowReducers2[i][j];
//          NewRow2[j] = add(NewRow2[j], mul(c2,RowReducers2[i][j]));
        }
        ++IterCount;
        if (IterCount == ModP.myMaxIters())
        {
          IterCount = 0;
          for (long j=i+1; j < NumPts; ++j)
          {
///JAA            NewRow[j] = ModP.myHalfNormalize(NewRow[j]);
            NewRow2[j] = ModP.myHalfNormalize(NewRow2[j]);
///JAA            CoCoA_ASSERT(ModP.myExport(ModP.myNormalize(NewRow[j])) == ModP.myExport(ModP.myNormalize(NewRow2[j])));
          }
        }
      }
      polyx.push_back(NewPolyx);
///JAA      NewRow[NumPts-1] = ModP.myNormalize(NewRow[NumPts-1]);
      NewRow2[NumPts-1] = ModP.myNormalize(NewRow2[NumPts-1]);
///JAA      CoCoA_ASSERT(ModP.myExport(NewRow[NumPts-1]) == ModP.myExport(ModP.myNormalize(NewRow2[NumPts-1])));
      // CoCoA_ASSERT(NewRow[NumPts-1] == 0);
///JAA      RowReducers.push_back(NewRow);
      RowReducers2.push_back(vector<SmallFpImpl::value>(NumPts));
    }
  }




    // currently just use random points
  vector<SmallFpImpl::value> GetNewSamplePt3(const SmallFpImpl& ModP, const vector<RingElem>& ParamDescr, const vector< vector<SmallFpImpl::value> >& SamplePts)
    {
      const ring& R = owner(ParamDescr[0]);
      const int NumParams = NumIndets(R);
      vector<RingElem> params(NumParams);
      vector<RingElem> ImagePt(len(ParamDescr));
      vector<SmallFpImpl::value> CandidatePt(len(ParamDescr));
      TryAgain:
      for (int i=0; i < NumParams; ++i)
      {
        params[i] = RingElem(CoeffRing(R),RandomLong(0,9999));
      }
      for (int i=0; i < len(ParamDescr); ++i)
        ImagePt[i] = eval(ParamDescr[i], params);
      for (int i=0; i < len(ParamDescr); ++i)
      {
        CandidatePt[i] = ModP.myReduce(ConvertTo<long>(ImagePt[i]));
      }
      if (find(SamplePts.begin(), SamplePts.end(), CandidatePt) != SamplePts.end()) goto TryAgain;
      return CandidatePt;
    }

  RingElem ImplicitByPoints3(const std::vector<RingElem>& ParamDescr)
  {
    using std::list;
    if (ParamDescr.empty() || len(ParamDescr) != 1+NumIndets(owner(ParamDescr[0])))
      CoCoA_THROW_ERROR(ERR::BadArg, "ImplicitByPoints");
    const ring& P = owner(ParamDescr[0]);
    const long p = ConvertTo<long>(characteristic(P));
    CoCoA_ASSERT(IsPrime(p));
    SmallFpImpl ModP(p);
    const long n = len(ParamDescr);
    ring Kx = NewPolyRing(CoeffRing(P), SymbolRange("x", 1, n));
    RingHom phi = CoeffEmbeddingHom(Kx);
    vector<PPMonoidElem> VecLPP; VecLPP.reserve(len(ParamDescr)); for(int i=0;i<len(ParamDescr);++i)VecLPP.push_back(LPP(ParamDescr[i]));
    map<PPMonoidElem,int> FindFactorIndex;

    QBGenerator QBG(PPM(Kx));
    const PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);
    vector<long> QBFactorPP;  QBFactorPP.push_back(-1);
    vector<long> QBFactorVar; QBFactorVar.push_back(-1);

//POLY    const RingElem one = monomial(Kx, PP1);
//POLY    vector<RingElem> polyx; polyx.push_back(one);
    FindFactorIndex[one(PPM(Kx))] = 0;
    vector< vector<SmallFpImpl::value> > PolyCoeffs; PolyCoeffs.push_back(vector<SmallFpImpl::value>(1,ModP.myReduce(1)));

    vector< vector<SmallFpImpl::value> > RowReducers(1); // essentially the BM matrix

    vector< vector<SmallFpImpl::value> > SamplePts;  // should be a std::set!!!
    int NumPts = 0;

    while (true)
    {
      ++NumPts;
//      if (NumPts%10 == 0) clog<<"NumPts="<<NumPts<<endl;
      const vector<SmallFpImpl::value> NewSamplePt = GetNewSamplePt3(ModP, ParamDescr, SamplePts);
      SamplePts.push_back(NewSamplePt);

//POLY      vector<RingElem> NewPt(n); for(int i=0;i<n;++i)NewPt[i]=RingElem(CoeffRing(P),ModP.myExport(NewSamplePt[i]));

      vector<SmallFpImpl::value> QBvalue(len(QBG.myQB()));
      QBvalue[0] = ModP.myReduce(1);
      for (long i=1; i < len(QBG.myQB()); ++i)
        QBvalue[i] = ModP.myMul(QBvalue[QBFactorPP[i]], NewSamplePt[QBFactorVar[i]]);

      for (long i=0; i < NumPts; ++i)
      {
//POLY        RowReducers[i].push_back(ModP.myReduce(ConvertTo<long>(eval(polyx[i], NewSamplePt)))); // clear but inefficient
        SmallFpImpl::NonRedValue EvalAtNewPt = PolyCoeffs[i][0];
        long IterCount=0;
        const long ReduceNow = ModP.myMaxIters();
        for (int j=1;j<len(PolyCoeffs[i]);++j)
        {
///          EvalAtNewPt = add(EvalAtNewPt, mul(QBvalue[j], PolyCoeffs[i][j]));
          EvalAtNewPt += QBvalue[j] * PolyCoeffs[i][j];
          ++IterCount;
          if (IterCount == ReduceNow) { EvalAtNewPt = ModP.myHalfNormalize(EvalAtNewPt); IterCount = 0; }
//          EvalAtNewPt = ModP.myAdd(EvalAtNewPt, ModP.myMul(QBvalue[j],PolyCoeffs[i][j]));
        }

        RowReducers[i].push_back(ModP.myNormalize(EvalAtNewPt));
////        RowReducers[i].push_back(ModP.myReduce(ConvertTo<long>(eval(polyx[i], NewPt)))); // clear but inefficient
////        if (ModP.myExport(EvalAtNewPt) != ModP.myExport(RowReducers[i].back())) CoCoA_THROW_ERROR("BUNGLED", "eval check");
      }
///JAA      CoCoA_ASSERT(ModP.myExport(RowReducers[NumPts-1][NumPts-1]) == ModP.myExport(RowReducers2[NumPts-1][NumPts-1]));
//      if (IsZero(RowReducers[NumPts-1][NumPts-1])) return polyx[NumPts-1];
      if (IsZero(RowReducers[NumPts-1][NumPts-1]))
      {
        // convert answer to a polynomial
        RingElem ans(Kx);
        for (long i=0; i < len(PolyCoeffs[NumPts-1]); ++i)
          ans += monomial(Kx, ModP.myExportNonNeg(PolyCoeffs[NumPts-1][i]), QBG.myQB()[i]);
        return ans;
      }

      // pick corner elem giving smallest LPP in image
      list<PPMonoidElem>::const_iterator BestIter = QBG.myCorners().begin();
      PPMonoidElem BestLPP = ComputeLPP(*BestIter, VecLPP);

      for(list<PPMonoidElem>::const_iterator it=QBG.myCorners().begin(); it != QBG.myCorners().end(); ++it)
//      for(int i=1;i<len(QBG.myCorners());++i)
      {
        const PPMonoidElem tmp = ComputeLPP(*it, VecLPP);
        if (tmp < BestLPP) { BestLPP = tmp; BestIter = it; }
      }
      const PPMonoidElem PP = *BestIter;
      QBG.myCornerPPIntoQB(PP);
      const vector<PPMonoidElem>& QB = QBG.myQB();
      FindFactorIndex[PP] = len(QB)-1;

///JAA      vector<SmallFpImpl::value_t> NewRow;
      const vector<long> expv = exponents(PP);
      int BestVar = -1; // just to keep compiler quiet
      long BestFactorPP = -1;
      for (int var=0; var < n; ++var)
      {
        if (expv[var] > 0)
        {
          const int ReducerIndex = FindFactorIndex[PP/indet(PPM(Kx),var)];
          if (ReducerIndex > BestFactorPP) { BestFactorPP = ReducerIndex; BestVar = var; }
        }
      }
      QBFactorPP.push_back(BestFactorPP);
      QBFactorVar.push_back(BestVar);
      const PPMonoidElem x = indet(PPM(Kx), BestVar);
      vector<SmallFpImpl::NonRedValue> NewRow;
      for(long j=0; j < NumPts ; ++j)
      {
//        NewRow.push_back(ModP.myReduce(ConvertTo<long>(eval(NewPolyx, SamplePts[i]))));
        NewRow.push_back(ModP.myMul(SamplePts[j][BestVar], RowReducers[BestFactorPP][j]));
      }
//POLY      RingElem NewPolyx = polyx[BestFactorPP]*indet(Kx,BestVar);
      vector<SmallFpImpl::NonRedValue> NewPolyCoeffs2(1+NumPts);
////NOT  NEEDED      NewPolyCoeffs[NumPts] = ModP.myReduce(1);
      const vector<SmallFpImpl::value>& OldPolyCoeffs = PolyCoeffs[BestFactorPP];
      for(int it=0;it<len(OldPolyCoeffs);++it)
      {
        if (IsZero(OldPolyCoeffs[it])) continue;
        const int j = FindFactorIndex[x*QB[it]];
        NewPolyCoeffs2[j] = OldPolyCoeffs[it];
      }

      // just gaussian reduction -- we know that matrix is upper triangular
      long IterCount = 0;
      for (long i=0; i < NumPts; ++i)
      {
        const SmallFpImpl::value coordi = ModP.myNormalize(NewRow[i]);
        if (IsZero(coordi)) continue;
        const SmallFpImpl::value q = ModP.myNegate(ModP.myDiv(coordi, RowReducers[i][i]));
        NewRow[i] = zero(SmallFp); // USELESS, but reassuring????
//POLY        NewPolyx += ModP.myExport(q)*polyx[i];
        for (long j=0; j<len(PolyCoeffs[i]);++j)
        {
///          NewPolyCoeffs2[j] = add(NewPolyCoeffs2[j], mul(q,PolyCoeffs[i][j]));  /// SLLOOOOOOWWWWWW
          NewPolyCoeffs2[j] += q * PolyCoeffs[i][j];
        }
        for (long j=i+1; j < NumPts; ++j)
        {
///          NewRow[j] = add(NewRow[j], mul(q,RowReducers[i][j]));
          NewRow[j] += q * RowReducers[i][j];
        }
        ++IterCount;
        if (IterCount == ModP.myMaxIters())
        {
          IterCount = 0;
          for (long j=i+1; j < NumPts; ++j)
          {
            NewRow[j] = ModP.myHalfNormalize(NewRow[j]);
          }
        for (long j=0; j<=NumPts;++j)
        {
          NewPolyCoeffs2[j] = ModP.myHalfNormalize(NewPolyCoeffs2[j]);

        }
        }
      }
//POLY      polyx.push_back(NewPolyx);
      vector<SmallFpImpl::value> NewPolyCoeffs(NumPts+1); for(long j=0;j<=NumPts;++j) NewPolyCoeffs[j]=ModP.myNormalize(NewPolyCoeffs2[j]);
      PolyCoeffs.push_back(NewPolyCoeffs);
      NewRow[NumPts-1] = ModP.myNormalize(NewRow[NumPts-1]);
      // CoCoA_ASSERT(IsZero(ModP.myNormalize(NewRow[NumPts-1]) == 0);
      RowReducers.push_back(vector<SmallFpImpl::value>(NumPts));
    }
  }


  //----------------------------------------------------------------------
  //-- SLICING ALGORITHM
  //----------------------------------------------------------------------

  namespace // anonymous for file local fns
  {

  template <typename T>
  std::vector<T> first(const std::vector<T>& v, long n)
  {
    std::vector<T> w;
    if (n>len(v)) CoCoA_THROW_ERROR("vector too short", "first(v, n)");
    for (long i=0; i<n; ++i) w.push_back(v[i]);
    return w;
  }

  template <typename T>
  std::vector<T> last(const std::vector<T>& v, long n)
  {
    std::vector<T> w;
    if (n>len(v)) CoCoA_THROW_ERROR("vector too short", "first(v, n)");
    for (long i=len(v)-n; i<len(v); ++i) w.push_back(v[i]);
    return w;
  }

  template <typename T>
  std::vector<T> concat(const std::vector<T>& v1, const std::vector<T>& v2)
  {
    std::vector<T> w;
    for (long i=0; i<len(v1); ++i) w.push_back(v1[i]);
    for (long i=0; i<len(v2); ++i) w.push_back(v2[i]);
    return w;
  }

  }  // end of anonymous namespace

  //------------------------------
  class ImplicitMill
  {
  public:
    enum FinalCall_t { elim, elim1, elimth, IDWC, IDWCLPP, IDWCLPP2, IDWCOrd2 };
    
      
  public:
    ImplicitMill(const vector<RingElem>& ParamDescr,
                 long NumXEndRec,
                 const string& FinalCall);
    
    //  private:
    vector<RingElem> myParamDescr;
    long myNumXEndRec;
    FinalCall_t myFinalCall;
    ring myKt;
    ring myKx;
    ring myKxt;  // for elim(1)(t) with/without weights, or homogenizing indet
    RingHom myHomT_XT;
    RingHom myHomXT_X;
    vector<RingElem> myXEndRec; // myKx: x[0],..,x[myNumXEndRec]
    mutable vector<long> myMinNumSlices;
    mutable vector<long> myXValues;
  };


  ImplicitMill::ImplicitMill(const vector<RingElem>& ParamDescr,
                             long NumXEndRec,
                             const string& FinalCall):
    myHomT_XT(IdentityHom(owner(ParamDescr[0]))),
    myHomXT_X(IdentityHom(owner(ParamDescr[0])))
  {
    myNumXEndRec = NumXEndRec;
    const ring& Ktorig = owner(ParamDescr[0]);
    //    myKt = owner(ParamDescr[0]);
    long NumT = NumIndets(Ktorig);
    long NumX = len(ParamDescr);
    if (NumT != NumX-1) CoCoA_THROW_ERROR("NumT != NumX-1", "ImplicitMill");
    myKx = NewPolyRingForImplicit(CoeffRing(Ktorig), "x",NumX);
    myKt = NewPolyRingForImplicit(CoeffRing(Ktorig), "t",NumT);
    myParamDescr = PolyAlgebraHom(Ktorig,myKt,indets(myKt))( ParamDescr );

    myMinNumSlices = vector<long>(NumX, 0);
    myXValues      = vector<long>(NumX, 0);  // x[i] = n

    if (FinalCall == "elim") myFinalCall = elim;
    else if (FinalCall == "elimth") myFinalCall = elimth;
    else if (FinalCall == "elim1") myFinalCall = elim1;
    else if (FinalCall == "IDWC") myFinalCall = IDWC;
    else if (FinalCall == "IDWCLPP") myFinalCall = IDWCLPP;
    else if (FinalCall == "IDWCLPP2") myFinalCall = IDWCLPP2;
    else if (FinalCall == "IDWCOrd2") myFinalCall = IDWCOrd2;
    else CoCoA_THROW_ERROR("either elim(1)(th), IDWC, IDWCLPP2, or IDWCLPP ","ImplicitMill");

    myXEndRec = first(indets(myKx), myNumXEndRec);

    switch (myFinalCall)
    {
    case ImplicitMill::IDWC:
    case ImplicitMill::IDWCLPP:
    case ImplicitMill::IDWCLPP2:
    case ImplicitMill::IDWCOrd2:
      break;
    case ImplicitMill::elim:   //   ring for elim -->
    case ImplicitMill::elim1:   //   ring for elim -->
      {
      matrix W=NewDenseMat(RingQQ(), 1, NumT+NumXEndRec);
      if (myFinalCall==elim1)
        for (long i=0; i<NumXEndRec; ++i) SetEntry(W,0,i,1);
      else 
        for (long i=0; i<NumXEndRec; ++i) SetEntry(W,0,i,deg(ParamDescr[i]));
      for (long i=0; i<NumT; ++i) SetEntry(W, 0, NumXEndRec+i, 1);
      myKxt = NewPolyRingForImplicit(CoeffRing(myKt),
                                     NumXEndRec+NumT,
                                     NewMatrixOrdering(MakeTermOrdMat(W),1));
      
      myHomT_XT = PolyAlgebraHom(myKt, myKxt, last(indets(myKxt),NumT));
      myHomXT_X = PolyAlgebraHom(myKxt, myKx,
                                 concat(myXEndRec,
                                        vector<RingElem>(NumT,zero(myKx))));
      }
      break;
    case ImplicitMill::elimth:  //   ring for elim -->
      {
      matrix W=NewDenseMat(RingQQ(), 1, NumT+NumXEndRec+1);
      for (long i=0; i<NumXEndRec; ++i) SetEntry(W,0,i,deg(ParamDescr[i]));
      for (long i=0; i<=NumT; ++i) SetEntry(W, 0, NumXEndRec+i, 1);
      myKxt = NewPolyRingForImplicit(CoeffRing(myKt),
                                     NumXEndRec+NumT+1,
                                NewMatrixOrdering(MakeTermOrdMat(W),1));
      
      myHomT_XT = PolyAlgebraHom(myKt, myKxt, first(last(indets(myKxt),NumT+1),NumT));
      myHomXT_X = PolyAlgebraHom(myKxt, myKx,
                                 concat(myXEndRec,
                                        vector<RingElem>(NumT+1,one(myKx))));
      }
      break;
    default:
      CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "ImplicitMill ctor");
    }
  }


  //------------------------------

  RingElem reconstruction(const vector<RingElem>& F, const vector<RingElem>& L)
  {
    long d = len(L);
    if (d==1) return one(owner(F[0]));
    RingElem s(owner(F[0]));
    RingElem PrLWithout_i(owner(F[0]));
    for (long i=0; i<d; ++i)
    {
      PrLWithout_i = one(owner(F[0]));
      for (long j=0; j<d; ++j)  if (i!=j) PrLWithout_i *= L[j];
      s += (PrLWithout_i*F[i]) / NR(PrLWithout_i, vector<RingElem>(1,L[i]));
    }
    //    std::cout << " s = " << s << std::endl;
    
    return s;
  }


  namespace{ // anonymous

    RingElem CallElim(const ImplicitMill& IM, long /*NumX*/)
    {
      vector<RingElem> v;
      long j=0;
      for (; j< IM.myNumXEndRec; ++j)
        v.push_back(indet(IM.myKxt,j) - IM.myHomT_XT(IM.myParamDescr[j]));
      for (; j<NumIndets(IM.myKx); ++j)
        v.push_back(IM.myXValues[j] - IM.myHomT_XT(IM.myParamDescr[j]));
      ideal J = elim(ideal(v), IM.myHomT_XT(indets(IM.myKt) ));
//       //----------
//       std::cout << "\n---------\n" << v << "\n---------\n" << std::endl;
//         RingElem ff = ComputeHypersurface(v,
//                                         LPP(product(IM.myHomT_XT( indets(IM.myKt) ))));
//         if (len(gens(J)) != 1)
//           std::cout << "\n len(gens(J)) != 1" << std::endl;
//       if (IM.myHomXT_X(monic(gens(J)[0]))
//           != IM.myHomXT_X(monic(ff)))
//         std::cout << IM.myHomXT_X(monic(gens(J)[0])) << "\n != \n"
//                   << IM.myHomXT_X(monic(ff)) << std::endl;
//       RingElem ff = CallImplicitDirectWithCondLPP2(IM, NumX);
//       if (monic(IM.myHomXT_X(gens(J)[0])) != monic(ff))
//         std::cout << monic(IM.myHomXT_X(gens(J)[0])) << "\n != \n"
//                   << monic(ff) << std::endl;
      //----------
      //      return monic(IM.myHomXT_X(gens(J)[0]));
      return IM.myHomXT_X(monic(gens(J)[0]));
    }
    

    RingElem CallElimTH(const ImplicitMill& IM, long /*NumX*/)
    {
      vector<RingElem> v;
      long j=0;
      RingElem h = indets(IM.myKxt)[NumIndets(IM.myKxt)-1];
      for (; j< IM.myNumXEndRec; ++j)
        v.push_back(homog(indet(IM.myKxt,j) - IM.myHomT_XT(IM.myParamDescr[j]), h));
      for (; j<NumIndets(IM.myKx); ++j)
        v.push_back(homog(IM.myXValues[j] - IM.myHomT_XT(IM.myParamDescr[j]), h));
      //      std::cout <<" v =" << v << std::endl;
      RingElem ff = ComputeElimFirst(v,
                                     LPP(product(IM.myHomT_XT( indets(IM.myKt) ))));
      return IM.myHomXT_X(monic(ff));
    }
    

    RingElem CallImplicitDirectWithCond(const ImplicitMill& IM, long NumX)
    {
      vector<RingElem> cond;
      for (long i=NumX; i<NumIndets(IM.myKx); ++i)
        cond.push_back(IM.myParamDescr[i] - IM.myXValues[i]);
      RingElem f = monic(ImplicitDirectWithCond(first(IM.myParamDescr,NumX),cond));
      return PolyAlgebraHom(owner(f), IM.myKx, IM.myXEndRec)(f);
    }
    

    RingElem CallImplicitDirectWithCondLPP(const ImplicitMill& IM, long NumX)
    {
      vector<RingElem> cond;
      for (long i=NumX; i<NumIndets(IM.myKx); ++i)
        cond.push_back(IM.myParamDescr[i] - IM.myXValues[i]);
      RingElem f = monic(ImplicitDirectWithCondLPP(first(IM.myParamDescr,NumX), cond));
      return PolyAlgebraHom(owner(f), IM.myKx, IM.myXEndRec)(f);
    }
    

    RingElem CallImplicitDirectWithCondLPP2(const ImplicitMill& IM, long NumX)
    {
      vector<RingElem> cond;
      for (long i=NumX; i<NumIndets(IM.myKx); ++i)
        cond.push_back(IM.myParamDescr[i] - IM.myXValues[i]);
      RingElem f = monic(ImplicitDirectWithCondLPP2(first(IM.myParamDescr,NumX), cond));
      return PolyAlgebraHom(owner(f), IM.myKx, IM.myXEndRec)(f);
    }
    

    RingElem CallImplicitDirectWithCondOrd2(const ImplicitMill& IM, long NumX)
    {
      vector<RingElem> cond;
      for (long i=NumX; i<NumIndets(IM.myKx); ++i)
        cond.push_back(IM.myParamDescr[i] - IM.myXValues[i]);
      RingElem f = monic(ImplicitDirectWithCondOrd2(first(IM.myParamDescr,NumX), cond));
      return PolyAlgebraHom(owner(f), IM.myKx, IM.myXEndRec)(f);
    }
    

    RingElem SliceCoreRec(const ImplicitMill& IM, long NumX)
  {  // X = [x[0], .., x[NumX-1],   a[NumX], .., a[N-1]]
    VerboseLog VERBOSE("SliceCoreRec");
    ring KX = IM.myKx;
    ring KT = IM.myKt;
    if (NumX <= IM.myNumXEndRec)
      switch (IM.myFinalCall)
      {
      case ImplicitMill::elim:  return CallElim(IM, NumX);
      case ImplicitMill::elim1: return CallElim(IM, NumX);
      case ImplicitMill::elimth: return CallElimTH(IM, NumX);
      case ImplicitMill::IDWC: return CallImplicitDirectWithCond(IM, NumX);
      case ImplicitMill::IDWCLPP: return CallImplicitDirectWithCondLPP(IM,NumX);
      case ImplicitMill::IDWCLPP2: return CallImplicitDirectWithCondLPP2(IM, NumX);
      case ImplicitMill::IDWCOrd2: return CallImplicitDirectWithCondOrd2(IM, NumX);
      default: CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "SliceCoreRec");
      }
    RingElem x_n = indet(KX, NumX-1);
    vector<RingElem> Slices;
    vector<RingElem> ResultSlices;
    for (long i=1; i<=50; ++i)
    {
      //      std::cout << " " << i << "(" << NumX << ")" << std::flush;
      Slices.push_back(x_n - i);
      IM.myXValues[NumX-1] = i; // [0,.., 0, *i*, a[NumX], .., a[N-1]]
      ResultSlices.push_back(SliceCoreRec(IM, NumX-1));
      if (i < IM.myMinNumSlices[NumX-1]) continue;
      RingElem candidate = monic(reconstruction(ResultSlices, Slices));
      if (IM.myMinNumSlices[NumX-1]==0 &&
          deg(candidate,NumX-1) == len(Slices)-1)
      {
        //        std::cout << ":LowD" << std::flush;
        continue;
      }
      vector<RingElem> v, img(NumIndets(KX), zero(KT)); 
      for (long j=0; j<NumX; ++j) img[j] = IM.myParamDescr[j];
      for (long j=NumX; j<len(IM.myXValues); ++j) img[j] = IM.myXValues[j];
      for (long j=NumX; j<len(IM.myXValues); ++j) v.push_back(IM.myParamDescr[j]-IM.myXValues[j]);
      RingElem sbst = NF(PolyAlgebraHom(KX,KT,img)(candidate), ideal(KT,v));
      if (IsZero(sbst))
      {
        if (IM.myMinNumSlices[NumX-1]==0) IM.myMinNumSlices[NumX-1] = i-1;
        //        std::cout << "\n--> MinNumSlices = " << MinNumSlices;
        //        std::cout << "  time: " << CpuTime()-T << std::endl;
        return candidate;
      }
      VERBOSE(50) << "    !!!wrong candidate!!! one more slice...";
    }
    CoCoA_THROW_ERROR("More than 50 slices: probably bad slicer!","SliceCoreRec");
    return zero(KX);
  }
  }
  

  RingElem SliceCore(const vector<RingElem>& ParamDescr,
                     long RecDepth,
                     const string& FinalCall)
  {
    ImplicitMill IM(ParamDescr, len(ParamDescr)-RecDepth, FinalCall);
    return SliceCoreRec(IM, len(ParamDescr));
  }
  

namespace // anonymous for file local defns
{ 
  RingElem SliceCoreModp(SparsePolyRing QQx,
                         const vector<RingElem>& ParamDescr,
                         long RecDepth,
                         const string& FinalCall, long p)
  {
    const SparsePolyRing& QQt(owner(ParamDescr[0]));
    const long s = NumIndets(QQt);
    const SparsePolyRing Fpt(NewPolyRing_DMPII(NewZZmod(p), SymbolRange("t",1,s)));
    const RingHom phi = PolyRingHom(QQt,Fpt, QQEmbeddingHom(Fpt), indets(Fpt));
    const RingElem f = SliceCore(phi(ParamDescr), RecDepth, FinalCall);
    const SparsePolyRing& Fpx(owner(f));
    const PPMonoidHom psi(GeneralHom(PPM(Fpx), first(indets(PPM(QQx)),NumIndets(Fpx))));
    RingElem ans(QQx);
    for (SparsePolyIter it=BeginIter(f) ; !IsEnded(it) ; ++it )
      ans += monomial(QQx, ConvertTo<BigInt>(coeff(it)), psi(PP(it)));
    return ans;
  }


  bool IsZeroEvalHorner(ConstRefRingElem f, const vector<RingElem>& ParamDescr)
  {
    VerboseLog VERBOSE("IsZeroEvalHorner");
    const SparsePolyRing& P1 = owner(f);
    const SparsePolyRing& P2 = owner(ParamDescr[0]);
    SparsePolyRing P = NewPolyRing(CoeffRing(P1),
                                   NewSymbols(NumIndets(P1)+NumIndets(P2)));
    RingHom phi1 = PolyAlgebraHom(P1, P, first(indets(P), NumIndets(P1)));
    RingHom phi2 = PolyAlgebraHom(P2, P, last(indets(P), NumIndets(P2)));
    const vector<RingElem>& ParamDescr_P = phi2(ParamDescr);
    RingElem eval_f = phi1(f);
    for (long i=0; i<NumIndets(P1); ++i)
    {
      VERBOSE(50) << "horner-"<<i << std::endl;
      std::vector<RingElem> c = CoeffVecWRT(eval_f, indet(P,i));
      eval_f = zero(P);
      for (long d=len(c)-1; d>=0; --d)
        eval_f = eval_f*ParamDescr_P[i] + c[d];
    }
    return IsZero(eval_f);
  }
  


}  // namespace // anonymous


  

  RingElem SliceCoreQQ(const vector<RingElem>& ParamDescr,
                       long RecDepth,
                       const string& FinalCall)
  {
    VerboseLog VERBOSE("SliceCoreQQ");
    SparsePolyRing QQx(RingQQt(len(ParamDescr)));
    const SparsePolyRing& QQt = owner(ParamDescr[0]);
    RingHom eval = PolyAlgebraHom(QQx, QQt, ParamDescr);
    long PrimeCount = 0;
    long p = 46349;
    RingElem fCRT(QQx);
    BigInt fModulus(p);
//     RingElem ParDDenom = CommonDenom(ParamDescr[0]);
//     for (long i=1; i<len(ParamDescr); ++i)
//       ParDDenom = lcm(ParDDenom, CommonDenom(ParamDescr[i]));
    RingElem ParDDenom = CommonDenom(ParamDescr);
    while (IsZero(fCRT))
    {
      CheckForInterrupt("SliceCoreQQ");
      p = PrevPrime(p);
      VERBOSE(20) << ++PrimeCount << ": prime is " << p << std::endl;
      if (IsDivisible(ParDDenom, p))
      {
        VERBOSE(80) << "   UGLY PRIME: going to another prime" << std::endl;
        continue;
      }
      fCRT = SliceCoreModp(QQx, ParamDescr, RecDepth, FinalCall, p);
      fModulus = p;
      try
      {
        VERBOSE(80) << "----  RatReconstruct" << std::endl;
        const RingElem f = RatReconstructPoly(fCRT, fModulus);
        VERBOSE(80) << "----  IsZero(eval(f))" << std::endl;
        if (IsZero(eval(f))) return f;
      }
      catch (const CoCoA::ErrorInfo& err)
      {
        if (err != ERR::CannotReconstruct) throw;
      }
    }
    while (true)
    {
      CheckForInterrupt("SliceCoreQQ");
      p = PrevPrime(p);
      VERBOSE(20) << ++PrimeCount << ": prime is " << p << std::endl;
      if (IsDivisible(ParDDenom, p))
      {
        VERBOSE(80) << "   UGLY PRIME: going to another prime" << std::endl;
        continue;
      }
      const RingElem fp = SliceCoreModp(QQx,ParamDescr, RecDepth, FinalCall, p);
      CRTPoly(fCRT, fModulus,  fCRT, fModulus,  fp, BigInt(p));
      try
      {
        VERBOSE(80) << "----  RatReconstruct" << std::endl;
        const RingElem f = RatReconstructPoly(fCRT, fModulus);
        double t=CpuTime();
        bool b = IsZero(eval(f));
        VERBOSE(90) << "---- IsZero time: " << CpuTime()-t << std::endl;
        // t = CpuTime();
        // b = IsZeroEvalHorner(f, ParamDescr);
        // std::cout << std::boolalpha; // so that bools print out as true/false
        // std::cout << "---- IsZeroEvalHorner " << b
        //          << " time: " << CpuTime()-t << std::endl;
        if (b) return f;
      }
      catch (const CoCoA::ErrorInfo& err)
      {
        if (err != ERR::CannotReconstruct) throw;
      }
    }
  }
  
  


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-implicit.C,v 1.13 2024/07/15 17:11:24 bigatti Exp $
// $Log: SparsePolyOps-implicit.C,v $
// Revision 1.13  2024/07/15 17:11:24  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.12  2024/03/15 15:56:15  bigatti
// Summary: using function elim
//
// Revision 1.11  2023/03/13 21:06:40  abbott
// Summary: New iterator "for" syntax; some cruft removed where the "for" has already been changed.  Redmine 1346
//
// Revision 1.10  2023/03/09 22:41:07  abbott
// Summary: Tidied use of exponents (redmine 1694)
//
// Revision 1.9  2022/10/11 15:44:00  bigatti
// Summary: using new interface of exponents returning vector
//
// Revision 1.8  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.7  2022/02/04 21:41:06  abbott
// Summary: Changed name MakeTermOrd to MakeTermOrdMat (redine 854)
//
// Revision 1.6  2021/10/04 08:59:41  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.5  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.4  2019/12/21 16:43:43  abbott
// Summary: Cleaned (most ring& instead of ring); removed some cruft
//
// Revision 1.3  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.2  2018/05/17 15:50:13  bigatti
// -- renamed MatrixOperations --> MatrixOps
// -- sorted includes
// -- renamed VectorOperations --> VectorOps
// -- added include SparsePolyIter
//
// Revision 1.1  2018/04/06 15:16:11  bigatti
// -- renamed TmpImplicit.C
//
// Revision 1.50  2018/03/29 09:36:40  bigatti
// -- added member functions myTestIsRadical, myTestIsPrimary and flags
//
// Revision 1.49  2018/02/27 17:30:23  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.48  2018/02/27 10:55:37  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.47  2017/09/06 11:56:30  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.46  2017/07/07 09:52:59  bigatti
// -- changed NextPrime --> PrevPrime
//
// Revision 1.45  2017/06/28 12:48:18  bigatti
// -- all cout now are VERBOSE(..)
//
// Revision 1.44  2017/05/15 10:24:37  bigatti
// -- detecting UGLY PRIME using divisibility (instead of try/catch)
//
// Revision 1.43  2017/05/11 10:58:10  bigatti
// -- now catching ugly primes (as for MinPolyModular)
//
// Revision 1.42  2017/05/11 08:46:46  bigatti
// -- added error ERR::CannotReconstruct
// -- cleaned up code accordingly
//
// Revision 1.41  2017/04/18 09:52:51  bigatti
// -- added verbosity
//
// Revision 1.40  2017/04/10 15:51:19  bigatti
// -- debug printing set to false (to be substituted by verbose)
//
// Revision 1.39  2017/03/08 15:35:09  bigatti
// -- added verbosity (90)
//
// Revision 1.38  2016/11/03 12:25:25  abbott
// Summary: Changed IsRadical (for PPMonoidElem) into IsSqFree
//
// Revision 1.37  2016/07/19 11:35:57  bigatti
// -- static bool for printing algorithm name
//
// Revision 1.36  2016/06/24 14:27:41  bigatti
// -- renamed CRT_poly --> CRTPoly
//
// Revision 1.35  2016/04/26 14:33:00  abbott
// Summary: commented out some debug print commands
//
// Revision 1.34  2016/04/22 16:06:59  bigatti
// -- added first implementation of IsZeroEvalHorner
//
// Revision 1.33  2016/03/09 10:17:10  abbott
// Summary: Commented out some logging print commands
//
// Revision 1.32  2016/02/17 10:30:37  bigatti
// -- added a return to keep compiler quiet (NYI)
//
// Revision 1.31  2016/01/26 13:51:00  bigatti
// -- added ImplicitDirectOrd2
//
// Revision 1.30  2015/12/08 13:56:09  abbott
// Summary: Updated Mario's code!  Very many changes!
//
// Revision 1.29  2015/11/04 12:14:58  abbott
// Summary: Several consequential changes (after revising SmallFpImpl)
//
// Revision 1.28  2015/10/01 10:15:05  bigatti
// -- moved CRT_poly, RatReconstructPoly --> SparsePolyRing.C
//
// Revision 1.27  2015/06/26 14:59:51  abbott
// Summary: Revised after changing interface to ErrorInfo (inherited from CoCoA::exception)
// Author: JAA
//
// Revision 1.26  2015/06/11 16:57:24  bigatti
// -- using new functions monomial(ring, pp) and monomial(ring, expv)
//
// Revision 1.25  2015/05/04 13:17:16  bigatti
// -- just new names: FindReducerIndex, NewPolyx, polyx,...
//
// Revision 1.24  2015/04/24 15:34:02  bigatti
// -- added ImplicitDirectLPP2Homog, ImplicitDirectLPP2HomogBis
// -- changed ring names into Kx and Kt
//
// Revision 1.23  2015/03/09 13:13:27  bigatti
// -- added SliceCoreQQ (and all needed functions in anonymous namespace)
//
// Revision 1.22  2015/03/05 16:49:32  bigatti
// -- changed back to creating a new ring (faster) in ImplicitDirectWithCondLPP2
//
// Revision 1.21  2015/03/04 11:18:19  bigatti
// -- added elim1 (no weights), elimth (truncated + homogeneous)
//
// Revision 1.20  2015/02/12 15:59:03  bigatti
// -- in ImplicitDirectWithCondLPP2: ChooseIndet(..poly) --> ChooseIndet(..reducer)
//
// Revision 1.19  2015/02/12 09:38:11  bigatti
// -- changed all ring names into Kt and Kx
//
// Revision 1.18  2015/02/06 13:48:03  bigatti
// -- allowing direct computation in QQ[...] using DMPI (instead of just DMPII)
//
// Revision 1.17  2014/12/18 15:57:45  bigatti
// -- renamed myRingX and T into myKx myKt in ImplicitMill
// -- myKt is now DMPII
// -- fixed usage of myKt in ImplicitDirectWithCondLPP2
//
// Revision 1.16  2014/12/18 14:20:40  bigatti
// -- renamed rings into Kx, Kt in ImplicitDirectWithCondLPP2
// -- no new ring Kt (using the input one)
// -- indexing x from 1 in SliceCore
//
// Revision 1.15  2014/12/11 15:36:34  bigatti
// -- SliceCore output in "x" (instead of anonymous symbols)
//
// Revision 1.14  2014/12/10 11:55:21  bigatti
// -- new base ring for SliceCore
//
// Revision 1.13  2014/12/04 08:55:43  bigatti
// -- minor optimization to reduce (use myAddMul instead?)
//
// Revision 1.12  2014/12/03 15:55:31  bigatti
// -- trivial change in ChooseIndet
//
// Revision 1.11  2014/12/03 13:55:02  bigatti
// -- added ImplicitDirectWithCondLPP2
//
// Revision 1.10  2014/11/28 15:22:55  bigatti
// -- added arg to SliceCore for algorithm of final calls
// -- now SliceCore uses ImplicitMill
// -- minor fix in ImplicitDirectWithCondLPP
//    (now calling ImplicitDirectLPP is cond is empty)
//
// Revision 1.9  2014/11/27 11:32:25  abbott
// Summary: Added ImplicitDirectWithCondLPP
// Author: JAA
//
// Revision 1.8  2014/11/27 09:45:37  bigatti
// -- SliceCoreRec: choose CallElim or CallImplicitDirect (at compiletime)
//
// Revision 1.7  2014/11/17 10:34:53  abbott
// Summary: Added ImplicitDirectLPP2 and ImplicitDirectWithCond
// Author: JAA
//
// Revision 1.6  2014/11/13 16:39:52  bigatti
// -- SliceCore: creating NewPolyRing_DMPII, removed printouts
//
// Revision 1.5  2014/11/13 16:30:40  abbott
// Summary: Now ImplicitDirect uses DMPII
// Author: JAA
//
// Revision 1.4  2014/11/13 15:45:05  abbott
// Summary: update to ImplicitDirectLPP
// Author: JAA
//
// Revision 1.3  2014/11/11 09:19:59  bigatti
// -- added SliceCore
//
// Revision 1.2  2014/11/04 18:20:29  abbott
// Summary: Added anon namespaces
// Author: JAA
//
// Revision 1.1  2014/11/04 18:10:13  abbott
// Summary: new tmp code for implicitization
// Author: JAA
//
//
