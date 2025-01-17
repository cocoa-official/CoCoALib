//   Copyright (c)  2018  John Abbott and Anna M. Bigatti
//   Authors:  2018  Anna M. Bigatti
//             2017  Alice Moallemy first translation CoCoA5-->CoCoALib
//             2017  Elisa Palezzato original code in CoCoA5

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


// Source code for ideals in SparsePolyRing (functions and member functions)

#include "CoCoA/SparsePolyOps-ideal.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixOps.H" // for "-"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyOps-MinPoly.H" // for MinPolyDef
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"  // for myGcd
#include "CoCoA/random.H" // for RandomLongStream
#include "CoCoA/time.H"
#include "CoCoA/verbose.H"

#include <algorithm>
using std::max;
#include <iostream>
using std::endl;
#include <list>
using std::vector;

namespace CoCoA
{

  namespace  //anonymous  for QuotientBasis
  {

    inline bool IsDivisible(ConstRefPPMonoidElem pp, const std::list<PPMonoidElem>& ByL)
    {
      return any_of(ByL.begin(), ByL.end(), [&pp](const PPMonoidElem& t){return IsDivisible(pp,t);});
///      for (auto i=ByL.begin(); i != ByL.end(); ++i)
///        if (IsDivisible(pp, *i)) return true;
///      return false;
    }
    

    // // Apparently STL bind2nd cannot have (const?) reference arguments.
    // // Apparently BOOST bind would work: what about C++-11?
    // // This is a hack for STL
    // class IsDivHack
    // {
    // public:
    //   IsDivHack (const std::list<PPMonoidElem>& L): myL(L) {}
    //   bool operator () (ConstRefPPMonoidElem pp) {return IsDivisible(pp, myL);}
      
    // private:
    //   const std::list<PPMonoidElem> & myL;
    // };


    void QuotientBasisRec(std::vector<PPMonoidElem>& ans, 
                          const std::list<PPMonoidElem>& L, 
                          ConstRefPPMonoidElem prefix, 
                          long idx)
    {
      const PPMonoid& PPM = owner(L.front());
      const PPMonoidElem& X = indets(PPM)[idx];
      PPMonoidElem prefixXd(prefix);  // prefix * x[idx]^d
      long MaxDeg = 0;
      if (idx == NumIndets(PPM)-1)
      {
        MaxDeg = exponent((L.front()), idx);
        for (int d=0; d < MaxDeg;  ++d, prefixXd *= X)
          ans.push_back(prefixXd);
        return;
      }
      for (const PPMonoidElem& t: L)
        MaxDeg = max(MaxDeg, exponent(t,idx)); // if (exponent(t,idx) > MaxDeg)  MaxDeg = exponent(t,idx);
//      for (auto it=L.begin(); it != L.end() ; ++it)
//        if (exponent((*it),idx) > MaxDeg)  MaxDeg = exponent((*it),idx);
      std::list<PPMonoidElem> CutOff, tmp;
      PPMonoidElem Xd(PPM);  // x[idx]^d
      for (int d=0; d < MaxDeg; ++d, prefixXd *= X, Xd *= X)
      {
        for (const PPMonoidElem& t: L)
          if (exponent(t,idx) == d)  tmp.push_back(t/Xd);
//        for (auto it=L.begin(); it != L.end(); ++it)
//          if (exponent(*it,idx) == d)  tmp.push_back(((*it)/Xd));
        CutOff.remove_if([&tmp](const PPMonoidElem& t){return IsDivisible(t,tmp);}); //        CutOff.remove_if(IsDivHack(tmp));
        CutOff.splice(CutOff.end(), tmp);
        QuotientBasisRec(ans, CutOff, prefixXd, idx+1);
      }
    }
  } // anonymous namespace


  std::vector<PPMonoidElem> QuotientBasis(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (!IsZeroDim(I))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
    const vector<RingElem> GB = GBasis(I);
    std::list<PPMonoidElem> LeadingPPs;
    vector<PPMonoidElem> ans;
    for (long i=0; i < len(GB); ++i)  LeadingPPs.push_back(LPP(GB[i]));
    QuotientBasisRec(ans, LeadingPPs, PPMonoidElem(PPM(RingOf(I))), 0);
    return ans;
  }


  std::vector<PPMonoidElem> QuotientBasisSorted(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (!IsZeroDim(I))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
    vector<PPMonoidElem> QB = QuotientBasis(I);
    std::sort(QB.begin(), QB.end());
    return QB;
  }


//   namespace {  // anonymous ------------------------------
//     RingElem CoeffOfTermSparse(ConstRefRingElem f, ConstRefPPMonoidElem pp)
//     {
//       for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
//         if (PP(itf) == pp) return coeff(itf);
//       return zero(CoeffRing(owner(f)));
//     }

//     matrix MultiplicationMat(ConstRefRingElem f, const ideal& I)
//     {
//       std::vector<PPMonoidElem> QB = QuotientBasis(I);
//       SparsePolyRing P = owner(f);
//       matrix Mf = NewDenseMat(CoeffRing(P), len(QB), len(QB));
//       for (long j=0; j<len(QB); ++j)
//       {
//         RingElem tmpf = f;
//         P->myMulByPP(raw(tmpf), raw(QB[j]));
//         tmpf = NF(tmpf, I);
//         for (long i=0; i<len(QB); ++i)
//           SetEntry(Mf,i,j, CoeffOfTermSparse(tmpf,QB[i]));
//       }
//       return Mf;
//     }

//   } // anonymous end -------------------------------------------


  //--- 0-dimensional ideals: functions using MinPolyQuot ------------------
  // Functions for radical, IsRadical, IsPrimary
  //   were translated from CoCoA-5 by Alice Moallemy.

  void SparsePolyRingBase::IdealImpl::myTestIsMaximal_0dim() const
  {
    VerboseLog VERBOSE("myTestIsMaximal_0dim");
    CoCoA_ASSERT(IamZeroDim());
    const ideal I(const_cast<IdealImpl*>(this));
    SparsePolyRing P = RingOf(I);
    long d = len(QuotientBasis(I)); // should be MultiplicityQuot(I);
    RingElem MP(P);
    for (long i=NumIndets(P)-1; i>=0; --i)
    //    for (long i=0; i<NumIndets(P); ++i)
    {
      VERBOSE(20) << "trying " << indet(P,i) << std::endl;
      MP = MinPolyQuot(indet(P,i), I, indet(P,i));
      if (!IsIrred(MP)) { myAssignPrimeFlag(false); return; }
      if (deg(MP) == d) { myAssignMaximalFlag(true); return; }
    }
    // now we know  I  is radical
    BigInt CharP = characteristic(P);
    //  if (!IsZero(CharP) && CharP<2*d) // small characteristic: heuristic bound    
    if (!IsZero(CharP)) // finite characteristic    
    {
      if (!IsFiniteField(CoeffRing(P)))
        CoCoA_THROW_ERROR2(ERR::NYI, "infinite field, finite characteristic");
      VERBOSE(20) << "Calling FrobeniusMat... " << std::endl;
      if (NumCols(LinKer(FrobeniusMat(I)-IdentityMat(CoeffRing(P),d))) != 1)
        myAssignMaximalFlag(false);
      else
        myAssignMaximalFlag(true);
      return;
    }
    // characteristic 0
    for (long n=2*NumIndets(P); true; n += NumIndets(P))
    {
      RingElem L(RandomLinearForm(P, n));
      VERBOSE(20) << "trying " << L << std::endl;
      MP = MinPolyQuot(L, I, indet(P,0));
      if (!IsIrred(MP)) { myAssignPrimeFlag(false); return; }
      if (deg(MP) == d) { myAssignMaximalFlag(true); return; }
    }
  }


  namespace  // anonymous
  {

    bool IsPrimary_0dimFin(const ideal& I)
    {
      VerboseLog VERBOSE("IsPrimary_0dimFin");
      const PolyRing P(RingOf(I), CoCoA_ERROR_CONTEXT);
      const ring& K = CoeffRing(P);
      if (!IsFiniteField(K))
        CoCoA_THROW_ERROR2(ERR::BadRing, "CoeffRing must be finite field");
      if (!IsZeroDim(I))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);

      VERBOSE(40) << "Calling FrobeniusMat... " << std::endl;
      const auto F = FrobeniusMat(I);
      const int d = NumRows(F);
      const int DimKer = NumCols(LinKer(F -IdentityMat(K, d)));
      return (DimKer == 1);
    }
  }


  void SparsePolyRingBase::IdealImpl::myTestIsPrimary_0dim() const
  {
    VerboseLog VERBOSE("myTestIsPrimary_0dim");
    CoCoA_ASSERT(IamZeroDim());
    const ideal I(const_cast<IdealImpl*>(this));
    const SparsePolyRing P(RingOf(I), CoCoA_ERROR_CONTEXT);;
    long d = len(QuotientBasis(I)); // should be MultiplicityQuot(I)
    ideal RadPartial_withGB = I;
    ideal RadCurrent = I;
    bool IsOriginalIdeal = true;
    RingElem MP(P);
    double t0;
    double timeout_GB=0;
    double time_MinPoly=0;
    double total_time_MinPoly=0;
    for (long i=NumIndets(P)-1; i>=0; --i)
      //for (long i=0; i < NumIndets(P); ++i)
    {
      VERBOSE(20) << "trying " << indet(P,i) << endl;
      t0 = CpuTime();
      MP = MinPolyQuot(indet(P,i), RadPartial_withGB, indet(P,i));
      time_MinPoly = CpuTime()-t0;
      VERBOSE(90) << "time_MinPoly " << time_MinPoly << "s" << endl;
      total_time_MinPoly += time_MinPoly;
      timeout_GB = (total_time_MinPoly*i)/(NumIndets(P)-i) +10; // see redmine 1375
      const factorization<RingElem> F = factor(MP);
      if (len(F.myFactors()) != 1) { myAssignPrimaryFlag(false); return; }
      if (deg(MP) == d)
      {
        if (IsOriginalIdeal && F.myMultiplicities()[0]==1)
          myAssignMaximalFlag(true);
        myAssignPrimaryFlag(true);
        return;
      }
      if (F.myMultiplicities()[0] != 1) // if (!IsSqFree(MP))
      {
        VERBOSE(40) << "MinPoly(" << indet(P,i)
                    << ") not square-free; deg = " << deg(MP) << endl;
        myAssignRadicalFlag(false);
        MP = radical(MP);
        IsOriginalIdeal = false;
        RadCurrent = RadCurrent + ideal(MP);
        // if (i==NumIndets(P)-1) break;  // Seidenberg
        if (i==0) break;  // Seidenberg
        try
        {
          VERBOSE(40) << "-->>>  computing GBASIS with time limit " << timeout_GB << "s  <<<--" << endl;
          GBasisByHomog(RadCurrent, CpuTimeLimit(timeout_GB));
          VERBOSE(80) << "-->>>  GBasisByHomog completed  <<<--" << endl;
          d = len(QuotientBasis(RadCurrent)); // should be MultiplicityQuot(I)
          VERBOSE(40) << "new multiplicity = " << d
                      << ";  deg(rad(minpoly)) = " << deg(MP) << endl;
          RadPartial_withGB = IdealOfGBasis(RadCurrent);
          RadCurrent = RadPartial_withGB;
          if (deg(MP) == d) { myAssignPrimaryFlag(true); return; }
        }
        catch (const CoCoA::TimeoutException&)
        {
          VERBOSE(40) << ">>>  GBASIS timed out  <<<" << endl;
        }
      }
    }//end for
    const ideal J = RadCurrent;
    d = len(QuotientBasis(J)); // should be MultiplicityQuot(I)
    VERBOSE(30) << "now ideal J is radical" << endl;
    const long CharP = ConvertTo<long>(characteristic(P));
    if (CharP != 0)
    {
      if (!IsFiniteField(CoeffRing(P)))  CoCoA_THROW_ERROR2(ERR::NYI, "only for finite fields");
      // and if q = p^2 // ??? Anna 2020-02-12
      myAssignPrimaryFlag(IsPrimary_0dimFin(J));
      return;
    }
    long n = NumIndets(P);
    RingElem z = indet(P, 0); //indet for minpoly
    while (true)  //instead of repeat
    {
      RingElem L = RandomLinearForm(P, n);
      VERBOSE(20) << "-- trying " << L << endl;
      const RingElem MP = MinPolyQuot(L, J, z);
      if (!IsIrred(MP)) { myAssignPrimaryFlag(false); return; }
      if (deg(MP) == d) { myAssignPrimaryFlag(true); return; }
      n *= 2;
    }
  }
  

  void SparsePolyRingBase::IdealImpl::myTestIsRadical_0dim() const
  {
    VerboseLog VERBOSE("myTestIsRadical_0dim");
    CoCoA_ASSERT(IamZeroDim());
    const ideal I(const_cast<IdealImpl*>(this));
    SparsePolyRing P = RingOf(I);
    long d = len(QuotientBasis(I)); // should be MultiplicityQuot(I);
    RingElem MP(P);
    for (long i=NumIndets(P)-1; i>=0; --i)
      //    for (long i=0; i<NumIndets(P); ++i)
    {
      VERBOSE(20) << "trying " << indet(P,i) << endl;
      MP = MinPolyQuot(indet(P,i), I, indet(P,i));
      if (!IsSqFree(MP)) { myAssignRadicalFlag(false); return; }
      if (!IsIrred(MP))  myAssignPrimeFlag(false); // just assign, no "return"!
      else if (deg(MP) == d)  { myAssignMaximalFlag(true); return; }
      if (deg(MP) == d)  { myAssignRadicalFlag(true); return; }
    }
    myAssignRadicalFlag(true);
  }


  ideal SparsePolyRingBase::IdealImpl::myRadical_0dimDRL() const
  {
    VerboseLog VERBOSE("myRadical_0dimDRL");
    CoCoA_ASSERT(IamZeroDim());
    const ideal I(const_cast<IdealImpl*>(this));
    SparsePolyRing P = RingOf(I);
    long d = len(QuotientBasis(I)); // should be MultiplicityQuot(I);
    ideal RadCurrent = I;
    ideal RadPartial_withGB = I;
    //    bool IsOriginalIdeal = true;
    RingElem MP(P);
    double t0;
    double timeout_GB=0;
    double time_MinPoly=0;
    double max_time_MinPoly=0;
    // fixed_time is used in heuristic estimate for Timeout when computing intermediate GBs
    constexpr double fixed_time = 30; // min time-out when computing GB (see timeout_GB 10 lines below)
    for (long i=NumIndets(P)-1; i>=0; --i)
    {
      VERBOSE(20) << "trying " << indet(P,i) << endl;
      t0 = CpuTime();
      MP = MinPolyQuot(indet(P,i), RadPartial_withGB, indet(P,i));
      time_MinPoly = CpuTime()-t0;
      VERBOSE(90) << "time_MinPoly " << time_MinPoly << "s" << endl;
      if (time_MinPoly>max_time_MinPoly)  max_time_MinPoly = time_MinPoly; // see redmine 1375
      timeout_GB = max_time_MinPoly*i + fixed_time;                        // see redmine 1375
//       If (IsOriginalIdeal)
//       {
//         if (!IsIrred(MP))  myAssignPrimeFlag(false);
//         else if (deg(MP) == d)  myAssignMaximalFlag(true);
//       }
      if (!IsSqFree(MP))
      {
        VERBOSE(40) << "MinPoly(" << indet(P,i)
                    << ") not square-free; deg = " << deg(MP) << endl;
        myAssignRadicalFlag(false);
        MP = radical(MP);
        // IsOriginalIdeal = false;
        RadCurrent = RadCurrent + ideal(MP);
        //        if (deg(original MP) == d) break;????????? teoria
///JAA???        if (i==0) break;  // Seidenberg
        VERBOSE(40) << "-->>>  computing GBASIS with time limit " << timeout_GB << "s  <<<--" << endl;
        try 
        {
//           VERBOSE(80) << "computing GBasis ..." << endl;
//           GBasis(RadCurrent, CpuTimeLimit(time_MinPoly));
          GBasisByHomog(RadCurrent, CpuTimeLimit(timeout_GB));
          VERBOSE(80) << "-->>>  GBasisByHomog completed  <<<--" << endl;
          d = len(QuotientBasis(RadCurrent)); // should be MultiplicityQuot(RC);
          VERBOSE(40) << "new multiplicity = " << d
                      << ";  deg(rad(minpoly)) = " << deg(MP) << endl;
          //          RadPartial_withGB = RadCurrent;
          RadPartial_withGB = IdealOfGBasis(RadCurrent);
          RadCurrent = RadPartial_withGB;
        }
        catch (const CoCoA::TimeoutException&)
        {
          VERBOSE(40) << ">>>  GBASIS timed out  <<<" << endl;
        }
      }
      if (deg(MP) == d)  break;
    }
    ourGetPtr(RadCurrent)->myAssignRadicalFlag(true);
    return RadCurrent;
  }


  

//   bool IsPrimary_0dim(const ideal& I)   // paper P.20 Algo 7.16
//   {
//     VerboseLog VERBOSE("IsPrimary_0dim");
//     ideal J = I;
//     PolyRing P = RingOf(J);
//     if (!IsZeroDim(I))
//       CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
//     int d = len(QuotientBasis(I)); // should be MultiplicityQuot(I)

//     for (long i = 0; i < NumIndets(P); ++i)
//     {
//       VERBOSE(20) << "trying " << indet(P,i) << endl;
//       RingElem MP = MinPolyQuot(indet(P,i), J, indet(P,i));
//       const factorization<RingElem> F = factor(MP);
//       const vector<RingElem>& facF = F.myFactors();
//       const vector<long>& multF = F.myMultiplicities();
//       if (len(facF) != 1) return false;
//       if (multF[0] != 1)
//       {
//         VERBOSE(30) << "MinPoly not square-free" << endl;
//         J = IdealOfGBasis(J) + ideal(facF[0]);
//         d = len(QuotientBasis(J)); // should be MultiplicityQuot(I)
//         if (deg(facF[0]) == d)  return true;
//       }   
//     }//end for
//     VERBOSE(30) << "now ideal J is radical" << endl;
//     long CharP = ConvertTo<long>(characteristic(P));
//     if (CharP != 0)
//     {
//       if (!IsFiniteField(CoeffRing(P)))
//         CoCoA_THROW_ERROR2(ERR::NYI, "only for finite field of characteristic p");
//       return IsPrimary_0dimFin(J);
//     }
//     long n = NumIndets(P);
//     RingElem z = indet(P, 0); //indet for minpoly
//     while (true)  //instead of repeat
//     {
//       RingElem L = RandomLinearForm(P, n);
//       VERBOSE(20) << "-- trying " << L << endl;
//       RingElem MP = MinPolyQuot(L, J, z);
//       if (!IsIrred(MP))  return false;
//       if (deg(MP) == d)  return true;
//       n *= 2;
//     }
//   } // IsPrimary_0Dim(I)





  //---- PrimaryDecomposition0 ---------------------------------------------

// -- Splitting functions

  namespace  {// anonymous

    void PDSplittingRandomLinear(RingElem& poly, bool& IsTotalSplit, factorization<RingElem>& facs, const ideal& I)
    {
      VerboseLog VERBOSE("PDSplittingRandomLinear");
      SparsePolyRing P = RingOf(I);
      long d = len(QuotientBasis(I)); // should be MultiplicityQuot(I);
      RingElem MP;
      for (long n=2*NumIndets(P); true; n += NumIndets(P))
      {
        RingElem L(RandomLinearForm(P, n));
        VERBOSE(20) << "L = " << L << std::endl;
        MP = MinPolyQuot(L, I, indet(P,0));
        const factorization<RingElem> F = factor(MP);
        VERBOSE(20) << "deg(mp(" << L << ")) =" << deg(MP)
                    << "  mult =" << d
                    << "  len(fact) =" << len(F.myFactors()) << endl;
        if (deg(MP)==d)
        {
          swap(poly, L);
          IsTotalSplit = true;
          facs = F;
          return;
        }
        if (len(F.myFactors()) > 1)
        {
          swap(poly, L);
          IsTotalSplit = false;
          facs = F;
          return;
        }
      }
    }


//     matrix FrobeniusSpace(const ideal& I)
//     {
//       VerboseLog VERBOSE("FrobeniusSpace");
//       SparsePolyRing  P = RingOf(I);
//       double  t0 = CpuTime();
//       matrix  C = FrobeniusMat(I);
//       VERBOSE(90) << "Time FrobeniusMat: " << CpuTime()-(t0) << endl;
//       t0 = CpuTime();
//       std::vector<PPMonoidElem> QB = QuotientBasisSorted(I);
//       VERBOSE(90) << "Time QuotientBasis: " << CpuTime()-(t0) << endl;
//       t0 = CpuTime();
//       matrix LK = LinKer(C-IdentityMat(CoeffRing(P), len(QB)));
//       VERBOSE(90) << "Time LinKer: " << CpuTime()-(t0) << endl;
//       //        
//       matrix M(DenseMatrix(P, 1, len(QB)));
//       for (long i=0; i<len(QB); ++i)  SetEntry(M,0,i, monomial(P, QB[i]));
//       return M * matrix(P, LK);
//     }


    void PDSplittingCharp(RingElem& poly, bool& IsTotalSplit, factorization<RingElem>& facs, const ideal& I)
    {
      VerboseLog VERBOSE("PDSplittingCharp");
      SparsePolyRing P = RingOf(I);
      std::vector<PPMonoidElem> QB = QuotientBasisSorted(I);
      long d = len(QB); // should be MultiplicityQuot
      matrix FrobSp = LinKer(FrobeniusMat(I) -IdentityMat(CoeffRing(P), d));
      long NumComps = NumCols(FrobSp);
      if (NumComps==1)
      {
        poly = zero(P);
        //		  factorization :=record[factors:=[one(P)],multiplicities:=[1]],
        //		  IsTotalSplit := true];
        return;
      }
      double prob=0.0;
      BigInt q = characteristic(P);
      if (q >= NumComps)
        prob = ConvertTo<double>(factorial(q))/ConvertTo<double>(factorial(q-NumComps)*power(q,NumComps));
      VERBOSE(90) << "NumComps =" << NumComps << "  prob =" << prob << endl;
      RingElem poly_tmp(P);
      for (long i=0; i<len(QB); ++i)
        poly_tmp += monomial(P, FrobSp(i,0),QB[i]);
      if (IsConstant(poly_tmp))
        for (long i=0; i<len(QB); ++i)
          poly_tmp += monomial(P, FrobSp(i,1),QB[i]);
      VERBOSE(99) << "Frobenius: poly =" << poly_tmp << endl;
      RingElem MP = MinPolyQuot(poly_tmp, I, indet(P,1));
      factorization<RingElem> F = factor(MP);
      VERBOSE(20) << "deg(MP) =" << deg(MP)
                  << "  BestDeg =" << (q<NumComps ? q : BigInt(NumComps))
                  << "  len(facs) =" << len(F.myFactors()) << endl;
      poly = poly_tmp;
      IsTotalSplit = (deg(MP)==NumComps);
      facs = F;
    }


    void PDSplitting(RingElem& poly, bool& IsTotalSplit, factorization<RingElem>& facs, const ideal& I)
    {
      VerboseLog VERBOSE("PDSplitting");
      SparsePolyRing P = RingOf(I);
      long d = len(QuotientBasis(I)); // should be MultiplicityQuot(I);
      RingElem MP(P);
      RingElem poly_tmp(P);
      factorization<RingElem> facs_tmp(one(P));  // myfactor for alg ext
      for (long i=NumIndets(P)-1; i>=0; --i)
      {
        poly_tmp = indet(P,i);
        VERBOSE(80) << "trying " << poly_tmp << std::endl;
        MP = MinPolyQuot(poly_tmp, I, poly_tmp);
        facs_tmp = factor(MP);
        VERBOSE(80) << "deg(mp(" << poly_tmp << ")) =" << deg(MP)
                    << "  mult =" << d
                    << "  len(facs) =" << len(facs_tmp.myFactors()) << endl;
        if (deg(MP)==d)
        {
          swap(poly, poly_tmp);
          facs = facs_tmp;
          IsTotalSplit = true;
          return;
        }
        if (len(facs_tmp.myFactors()) > 1)
        {
          swap(poly, poly_tmp);
          facs = facs_tmp;
          IsTotalSplit = false;
          return;
        }
      }
      if (characteristic(P)>0)
      {
        PDSplittingCharp(poly, IsTotalSplit, facs, I);
        return;
      }
      if (IsPrimary(I))
      {
        poly = zero(P);
        return;
      }
      //      CoCoA_THROW_ERROR2(ERR::NYI, "PDSplittingRandomLinear");
      PDSplittingRandomLinear(poly, IsTotalSplit, facs, I);
    }



// GBasisByHomog implemented in SparsePolyOps_ideal:
  
// define GBasisByHomog(ref I)
//   P := RingOf(I);
//   TmpP := NewPolyRing(CoeffRing(P), SymbolRange("Z",1,NumIndets(P)+1));
//   phi := PolyAlgebraHom(P, TmpP, first(indets(TmpP), NumIndets(P)));
//   psi := PolyAlgebraHom(TmpP, P, concat(indets(P), [one(P)]));
//   GensQi := gens(I);
//   I := ideal(psi( GBasis(ideal(homog(phi(GensQi), last(indets(TmpP))))) ));
//   return GBasis(I);
// enddefine; -- GBasisByHomog

  }
  

  namespace // anonymous
  {

    // similar to IsZeroEvalUniPolyMod in SparsePolyOps-MinPoly
    RingElem EvalUniPolyMod(ConstRefRingElem f, ConstRefRingElem g, const ideal& I)
    {
      const ring& Pf = owner(f);
      const ring& Pg = owner(g);
      if (CoeffRing(Pf)!=CoeffRing(Pg))  CoCoA_THROW_ERROR1(ERR::MixedCoeffRings);
      if (Pg!=RingOf(I))  CoCoA_THROW_ERROR2(ERR::MixedRings, "arg 2 and 3");
      const RingElem z = indet(Pf, UnivariateIndetIndex(f));
      const std::vector<RingElem> c = CoeffVecWRT(f, z);
      const long NF_freq = (NumIndets(Pg)>2)?1:2;
      RingElem eval_f = zero(Pg);
      for (long d=len(c)-1; d>=0; --d)
      {
        eval_f = eval_f*g + c[d];
        if (d%NF_freq == 0)  // NF every NF_freq steps?  Worth it?
          eval_f = NF(eval_f, I);
      }
      return eval_f;
    }
  }  // end of anon namespace


  void SparsePolyRingBase::IdealImpl::myPrimaryDecompositionCore_0dim(bool& IsCertified, std::vector<ideal>& Q) const
  {
    VerboseLog VERBOSE("myPrimaryDecompositionCore_0dim");
    SparsePolyRing P = myRing();
    RingElem SplitPoly(P);
    factorization<RingElem> facs(one(P));
    bool IsTotalSplit;
    double T_split = CpuTime();
    const ideal I(const_cast<IdealImpl*>(this));
    PDSplitting(SplitPoly, IsTotalSplit, facs, I);
    VERBOSE(90) << "Time PDSplitting: " << CpuTime()-(T_split) << endl;
    IsCertified = IsTotalSplit;
    if (IsZero(SplitPoly) || len(facs.myFactors()) == 1)
    {
      IamPrimary3Flag = true3;
      Q = std::vector<ideal>(1, ideal(const_cast<IdealImpl*>(this)));
      return;
    }
    IamPrimary3Flag = false3;
    double T_eval=0, T_power=0, T;
    vector<ideal> Qi;
    RingElem fac;
    for (long i=0; i<len(facs.myFactors()); ++i)
    {
      T = CpuTime();
      fac = power(facs.myFactors()[i], facs.myMultiplicities()[i]);
      T_power += CpuTime()-(T);
      if (!IsIndet(SplitPoly))
      {
        T = CpuTime();
        fac = EvalUniPolyMod(fac, SplitPoly, I);
        T_eval += CpuTime()-(T);
      }
      Qi.push_back(IdealOfGBasis(I) + ideal(fac));
    }
    VERBOSE(90) << "Time power:           " << T_power << endl;
    VERBOSE(90) << "Time EvalUniPolyMod:  " << T_eval << endl;
    Q = Qi;
    //    swap(Q, Qi);
  }

  
  std::vector<ideal> SparsePolyRingBase::IdealImpl::myPrimaryDecomposition_0dim() const
  {
    VerboseLog VERBOSE("myPrimaryDecomposition_0dim");
    if (IsTrue3(IamPrimary3()))
      return std::vector<ideal>(1, ideal(const_cast<IdealImpl*>(this)));
    std::vector<ideal> Q;
    bool IsCertified;
    myPrimaryDecompositionCore_0dim(IsCertified, Q);
    if (IsCertified)
    {
      VERBOSE(50) << "certified: " << len(Q) << " component(s)" << endl;
      //      VERBOSE(99) << "Q = " << Q << endl;
      return Q;
    }
    else
    {
      BigInt c = characteristic(myRing());
      double T;
      std::vector<ideal> PD;
      bool IsPrimaryQi = false;
      for (long i=0; i<len(Q); ++i)
      {
        VERBOSE(50) << "--checking: IsPrimary(Q[" << i << "])" << endl;
        T = CpuTime();
        vector<RingElem> GB = GBasisByHomog(Q[i]);
        // vector<RingElem> GB = GBasis(Q[i]); // may be very slow
        VERBOSE(90) << "Time GBasis(Q[i]):    " << CpuTime()-(T) << endl;
        if (c==0)
        {
          T = CpuTime();
          IsPrimaryQi = IsPrimary(Q[i]); // generally faster than rec call
          VERBOSE(90) << "Time IsPrimary(Q[i]): " << CpuTime()-T << endl;
        }
        if (IsPrimaryQi)
          PD.push_back(Q[i]);
        else
        {
          std::vector<ideal> PDQi=SparsePolyRingBase::IdealImpl::ourGetPtr(Q[i])->myPrimaryDecomposition_0dim();
          PD.reserve(len(PD)+len(PDQi));
          PD.insert(PD.end(), PDQi.begin(), PDQi.end());
        }
      }
      return PD;
    }
  }


} // end of namespace CoCoA
