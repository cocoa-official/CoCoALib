//   Copyright (c)  2020  John Abbott,  Anna Bigatti
//   Original author: 2020 Julian Danner (interreduced transcoded from CoCoA-5)

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


//#include "CoCoA/SparsePolyOps-vector.H"
#include "CoCoA/SparsePolyOps-RingElem.H"

#include "CoCoA/PPMonoidHom.H"
#include "CoCoA/ReductionCog.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/geobucket.H" // for myMul
#include "CoCoA/interrupt.H"
#include "CoCoA/verbose.H"


#include <algorithm>
using std::max;     // for CoeffVecWRT
#include <map>
using std::map;
#include <utility>
using std::make_pair;
using std::pair;
//#include <vector>
using std::vector;


namespace CoCoA
{

  // Deliberately not defined!
  // void interreduce(std::vector<RingElem>& v)
  // {
  //   std::swap(v, interreduced(v));
  // }


  // naive impl (orig transcoded from CoCoA-5 by Julian Danner)
  std::vector<RingElem> interreduced(std::vector<RingElem> v)
  {
    if (v.empty())  { return v; } // ??? or error???
    if (!HasUniqueOwner(v))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    VerboseLog VERBOSE("interreduced");
    //delete possible zeros in v
    const ring& P = owner(v[0]);
    const ring& k = CoeffRing(P);
    const bool CaseFp = IsFiniteField(k);
    const bool CaseQQ = IsQQ(k);
    v.erase(std::remove(v.begin(), v.end(), zero(P)), v.end());

    // this local fn is used in call to sort
    const auto CompareLPPs = [](const RingElem& f, const RingElem& g) { return LPP(f)<LPP(g); };
    long count = 0;
    while (true)
    {
      VERBOSE(90) << "round " << ++count << std::endl; // NB *always* incrs count!
      sort(v.begin(), v.end(), CompareLPPs);
      if (VerbosityLevel()>=99) std::cout << "L: " << v << std::endl;
      vector<RingElem> ans;
      RingElem rem;
      bool NewLPPfound = false;
      for (const auto& f: v)
      {
        CheckForInterrupt("interreduced");
        rem = NR(f, ans);
        if (IsZero(rem))  { continue; }
        const bool LPP_changed = (LPP(rem) != LPP(f));
        if (LPP_changed)
        {
          // Tidy up the new poly (redmine 1642)  ??? Worth it ???
          if (CaseFp)  { rem = monic(rem); }  // Anna: should this be done only before return?
          else if (CaseQQ)   { rem = prim(rem);}
        }
        ans.push_back(rem);
        NewLPPfound = (NewLPPfound || LPP_changed);
      }
      if (!NewLPPfound)  { return ans; }
      swap(v, ans); // quicker than: v = ans;
    }
  }


  // (by Anna M Bigatti inspired by proof of termination)
  // it does not seem to work well: sometimes a little faster, sometimes a lot slower
  std::vector<RingElem> interreduced_LT(std::vector<RingElem> v)
  {
    if (v.empty())  { return v; } // ??? or error???
    if (!HasUniqueOwner(v))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    VerboseLog VERBOSE("interreduced_LT");
    //delete possible zeros in v
    const ring& P = owner(v[0]);
    const ring& k = CoeffRing(P);
    const bool CaseFp = IsFiniteField(k);
    const bool CaseQQ = IsQQ(k);
    v.erase(std::remove(v.begin(), v.end(), zero(P)), v.end());

    // this local fn is used in call to sort
    const auto CompareLPPs = [](const RingElem& f, const RingElem& g) { return LPP(f)<LPP(g); };
    long count = 0;
    while (true)
    {
      ++count
      VERBOSE(90) << "round " << count << std::endl;
      if (VerbosityLevel()>=99) std::cout << "L: " << v << std::endl;
      sort(v.begin(), v.end(), CompareLPPs);
      vector<RingElem> ans;
      RingElem rem;
      bool NewLPPfound = false;
      for (const auto& f: v)
      {
        CheckForInterrupt("interreduced_LT");
        rem = NR_LT(f, ans);
        if (IsZero(rem))  { continue; }
        const bool LPP_changed = (LPP(rem) != LPP(f));
        if (LPP_changed)
        {
          // Tidy up the new poly (redmine 1642)  ??? Worth it ???
          if (CaseFp)  { rem = monic(rem); }  // Anna: should this be done only before return?
          else if (CaseQQ)   { rem = prim(rem);}
        }
        ans.push_back(rem);
        NewLPPfound = (NewLPPfound || LPP_changed);
      }
      if (!NewLPPfound)  { return ans; }
      swap(v, ans); // quicker than: v = ans;
    }
  }


  RingElem IndetsProd(const std::vector<RingElem>& L)
  {
    if (L.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    const ring& P = owner(L[0]);
    if (!IsSparsePolyRing(P))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const vector<long> VarIndices = IndetsIn(L);
    RingElem ans = one(P);
    for (auto k: VarIndices)
      ans *= indet(P,k);
    return ans;
  }


  // indices of indets appearing in (non-empty) L
  std::vector<long> IndetsIn(const std::vector<RingElem>& L)
  {
    if (L.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!IsSparsePolyRing(owner(L[0])))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRing& P = owner(L[0]);
    const int nvars = NumIndets(P);
    const long n = len(L);
    for (long i=1; i < n; ++i)
      if (owner(L[i]) != P)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    vector<bool> IndetSeen(nvars);
    int NumSeen = 0;
    for (long i=0; i < n; ++i)
    {
      const vector<long> IndetsInThisPoly = IndetsIn(L[i]);
      for (long xj: IndetsInThisPoly)
      {
        if (IndetSeen[xj]) continue;
        IndetSeen[xj] = true;
        if (++NumSeen == nvars) goto double_break;
      }
    }
    double_break:
    vector<long> ans; ans.reserve(nvars); // potentially wasteful reserve
    for (int i=0; i < nvars; ++i)
      if (IndetSeen[i]) ans.push_back(i);
    return ans;
  }


  //----------------------------------------------------------------------
  //??? the following functions to compute NR will be replaced by GBMill

  int FindReducerIndex(ConstRefPPMonoidElem pp, const std::vector<RingElem>& v)
  {
    const long nelems = len(v);
    for (long i=0; i < nelems; ++i)
      if (IsDivisible(pp, LPP(v[i])))
        return i;
    return -1;
  }


  inline int FindReducerIndex(const ReductionCog& F, const std::vector<RingElem>& v)
  {
    if ( IsActiveZero(F) ) return -1;
    return FindReducerIndex(ActiveLPP(F), v);
  }


  void ReduceActiveLM(ReductionCog& F, const std::vector<RingElem>& v)
  {
    int i;
    while ( (i = FindReducerIndex(F, v) ) != -1)
    {
      CheckForInterrupt("ReduceActiveLM");
      F->myReduce(v[i]);
    }
  }


  void reduce(ReductionCog& F, const std::vector<RingElem>& v)
  {
    ReduceActiveLM(F, v);
    while ( !IsActiveZero(F) )
    {
      F->myMoveToNextLM();
      ReduceActiveLM(F, v);
    }
  }
  //--------------------------------------------

  RingElem NR(ConstRefRingElem f, const std::vector<RingElem>& v)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P) )  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (!HasUniqueOwner(v) || (!v.empty() && owner(v.front()) != P))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if ( IsZero(f) ) return f;
    if ( v.empty() ) return f;
    RingElem ans(f);
    ReductionCog F = NewRedCogGeobucketField(owner(ans));
    F->myAssignReset(ans);
    reduce(F, v);
    F->myRelease(ans);
    return ans;
  }

  RingElem NR_LT(ConstRefRingElem f, const std::vector<RingElem>& v)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P) )  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (!HasUniqueOwner(v) || (!v.empty() && owner(v.front()) != P))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if ( IsZero(f) ) return f;
    if ( v.empty() ) return f;
    RingElem ans(f);
    ReductionCog F = NewRedCogGeobucketField(owner(ans));
    F->myAssignReset(ans);
    ReduceActiveLM(F, v);
    F->myRelease(ans);
    return ans;
  }

  // naive implementation
  namespace { // anonymous

    RingElem DivAlgLM(vector<RingElem>& QuotRem, ConstRefRingElem f, const vector<RingElem>& v)
    {
      const SparsePolyRing P = owner(f);
      RingElem m(P);
      RingElem r =f;
      
      long j = FindReducerIndex(LPP(r), v);
      while (j != -1)
      {
        P->myDivLM(raw(m), raw(r), raw(v[j])); // m =LM(r)/LM(v[i]); no checks
        QuotRem[j] += m;
        r -= m * v[j];
        if (IsZero(r)) break;
        j = FindReducerIndex(LPP(r), v);
      }
      return r;
    }

  } // anonymous
  
  // RingElem NormalRemainder(ConstRefRingElem f, const vector<RingElem>& v)
  // {
  //   if (IsZero(f)) return f;
  //   const SparsePolyRing P = owner(f);
  //   RingElem ansNR(P);
  //   RingElem tmpNR =f;
  
  //   tmpNR = NRLM(f, v);
  //   while (!IsZero(tmpNR))
  //   {
  //     P->myMoveLMToBack(raw(ansNR), raw(tmpNR));
  //     tmpNR = NRLM(tmpNR, v);
  //   }
  //   return ansNR;
  // }


  std::vector<RingElem> TmpDivAlg(ConstRefRingElem f, const std::vector<RingElem>& v)
  {
    const SparsePolyRing& P = owner(f);
    if (!IsPolyRing(P) )  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (!HasUniqueOwner(v) || (!v.empty() && owner(v.front()) != P))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    //    if ( IsZero(f) ) return f;
    //    if ( v.empty() ) return f;
    RingElem ansNR(P);
    RingElem tmpNR =f;
    vector<RingElem> QuotRem(len(v)+1, zero(P));
    while (!IsZero(tmpNR))
    {
      tmpNR = DivAlgLM(QuotRem, tmpNR, v);
      if (IsZero(tmpNR)) break;
      P->myMoveLMToBack(raw(ansNR), raw(tmpNR)); // ansNR+=LM(tmpNR)
    }
    swap(QuotRem[len(v)], ansNR);
    return QuotRem;
  }

  ///----------------  coefficients vec -------------------------------------------

  
  namespace // anonymous for file local defns
  {

    class ByDecreasingPP
    {
    public:
      bool operator()(const PPMonoidElem& A, const PPMonoidElem& B) const
        {
          return A > B;
        }
    };

    // This fn is needed in a call to std::transform
    CoeffPP CoeffPPCtor(const pair<PPMonoidElem, RingElem>& arg)
    {
      return CoeffPP(arg.second, arg.first);
    }

  } // end of namespace anonymous


  std::vector<CoeffPP> CoefficientsWRT(ConstRefRingElem f, const std::vector<long>& indets)
  {
    if (!IsSparsePolyRing(owner(f)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRing P = owner(f);
    for (long i=0; i < len(indets); ++i)
      if (indets[i] < 0 || indets[i] >= NumIndets(P))  CoCoA_THROW_ERROR1(ERR::BadIndex);

    // Force the sorting criterion in the map
    typedef map<PPMonoidElem, RingElem, ByDecreasingPP> CoeffTable_t;
    CoeffTable_t CoeffTable;
    PPMonoidHom projection = RestrictionHom(PPM(P), indets);
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      const PPMonoidElem t = projection(PP(it));
      CoeffTable_t::iterator pos = CoeffTable.find(t);
      if (pos == CoeffTable.end())
      {
        CoeffTable.insert(make_pair(t, monomial(P, coeff(it), PP(it)/t)));
        continue;
      }
      RingElem m = monomial(P, coeff(it), PP(it)/t);
      P->myMoveLMToBack(raw(pos->second), raw(m));
//      pos->second += monomial(P, coeff(it), PP(it)/t);
    }
    vector<CoeffPP> ans; ans.reserve(len(CoeffTable));
    transform(CoeffTable.begin(), CoeffTable.end(), back_inserter(ans), CoeffPPCtor);
    // NOTE: CoeffTable will automatically be in decreasing order by PP!
    // (see 23.2.4/10 in C++11)
    return ans;
  }


  std::vector<CoeffPP> CoefficientsWRT(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (owner(f) != owner(x))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsSparsePolyRing(owner(f)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    vector<long> indices(1);
    if (!IsIndet(indices[0], x))  CoCoA_THROW_ERROR1(ERR::ReqIndet);
    return CoefficientsWRT(f, indices);
  }


  std::vector<RingElem> CoeffVecWRT(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (owner(f) != owner(x))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsSparsePolyRing(owner(f)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (IsZero(f)) return vector<RingElem>();
    vector<long> indices(1);
    if (!IsIndet(indices[0], x))  CoCoA_THROW_ERROR1(ERR::ReqIndet);
    vector<CoeffPP> CoeffPPList = CoefficientsWRT(f, indices);
    long Degf = 0;
    for (const CoeffPP& term: CoeffPPList)
      Degf = max(Degf, StdDeg(term.myPP));
//    for (int i=0; i < len(CoeffList); ++i)
//      Degf = max(Degf, StdDeg(CoeffList[i].myPP));
    vector<RingElem> ans(1+Degf, zero(owner(f)));
    for (CoeffPP& term: CoeffPPList)
      swap(ans[StdDeg(term.myPP)], term.myCoeff); // swap avoids a wasteful copy!
//    for (vector<CoeffPP>::iterator it=CoeffList.begin(); it != CoeffList.end(); ++it)
//      swap(ans[StdDeg(it->myPP)], it->myCoeff); // swap avoids a wasteful copy!
    return ans;
  }


  // Assumes f is in span of the basis; BasisPoly is a poly whose PPs are the actual basis.
  // Redmine 1210: do we want a version which also accepts PPvector???
  std::vector<RingElem> CoeffVecWRTSupport(ConstRefRingElem f, ConstRefRingElem BasisPoly)
  {
    const PolyRing& P = owner(f);
    if (P != owner(BasisPoly))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    const ring& K = CoeffRing(P);
//???    if (!IsField(K)) CoCoA_THROW_ERROR1(ERR::ReqField);

    const long n = NumTerms(BasisPoly);
    vector<RingElem> C(n, zero(K));
    if (IsZero(f)) return C;
    SparsePolyIter itf = BeginIter(f);
    PPMonoidElem PPf = PP(itf);

    long i=n-1;
    for (SparsePolyIter it=BeginIter(BasisPoly); !IsEnded(it); ++it, --i)
    {
      if (PPf < PP(it)) continue;
      if (PPf != PP(it))  CoCoA_THROW_ERROR2(ERR::BadArg, "not in span");
      C[i] = coeff(itf);
      ++itf;
      if (IsEnded(itf)) return C; // NORMAL EXIT POINT!!
      PPf = PP(itf);
    }
    // Get here only if (!IsEnded(itf))
  CoCoA_THROW_ERROR2(ERR::BadArg, "not in span");
    return C; // NEVER EXECUTED; just to keep compiler quiet
  }
  
  
} // end of namespace CoCoA
