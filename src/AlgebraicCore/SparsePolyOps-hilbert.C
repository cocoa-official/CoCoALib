//   Copyright (c)  2006-2014  John Abbott and Anna Bigatti
//   Author: 2006-2014 Anna Bigatti

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


// Source code for Hilbert-Poincare Series

#include "CoCoA/SparsePolyOps-hilbert.H"

#include "CoCoA/BigInt.H"
#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/degree.H"
#include "CoCoA/factorization.H"
#include "CoCoA/ideal.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"
#include "CoCoA/verbose.H"  // for VerboseLog
#include "TmpHilbertDir/AnnaUtils.h"
#include "TmpHilbertDir/IVectors.h"
#include "TmpHilbertDir/TmpPoincareCPP.H"
#include "TmpHilbertDir/eterms.h"
#include "TmpHilbertDir/poincare.h"
#include "TmpHilbertDir/unipoly.h"

#include <vector>
using std::vector;
// No longer needed: unique_ptr (was auto_ptr)
// #include <memory>
// using std::unique_ptr;

namespace CoCoA
{

  HPSeries::HPSeries(ConstRefRingElem num, const factorization<RingElem>& den):
    myNum(num),
    myDenFactors(den)
//???    myDenFactors(den.myFactors, den.myMultiplicities, den.myRemainingFactor)
  {
    // check consistency of factorization ???
    // check consistency of num and den
    if (!den.myFactors().empty() && (owner(num) != owner(den.myRemainingFactor())))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
  }


  HPSeries::HPSeries(const vector<BigInt>& DenseRepr, const vector<long>& DenExponents, long shift):
    myNum(RingQQt(1)),
    myDenFactors(one(RingQQt(1)))
///???    myDenFactors(vector<RingElem>(), vector<long>(), one(RingQQt(1)))
  {
      PolyRing QQt = RingQQt(1);
      const RingElem t = indet(QQt,0);
      const RingElem ONE = one(QQt);
      // create the numerator polynomial from the dense vector of coefficients
      RingElem tpower = ONE;

      if (shift > 0) tpower = power(t,shift); // multiply num by t^shift
      for (long i=0; i<len(DenseRepr); ++i)
      {
        myNum += DenseRepr[i] * tpower;   // NB tpower = t^i
        tpower *= t;
      }

      // fill the denominator factorization
      const long n = len(DenExponents);
      long i = 0;
      while (i < n)
      {
        const long e = DenExponents[i];
        ++i;
        long mult = 1;
        while (i < n && DenExponents[i] == e)
        {
          ++i;
          ++mult;
        }
        myDenFactors.myAppend(ONE - power(t,e), mult);
      }
      // for negative shift we add an extra t^(-shift) factor to den
      // alternatively the numerator could become a Laurent poly
      if (shift < 0) myDenFactors.myAppend(power(t,-shift), 1);
  }


  std::ostream& operator<<(std::ostream& out, const HPSeries& S)
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << "HPSeries(myNum = " << S.myNum
               << ", myDenFactors = " << S.myDenFactors << ")";
  }

  //----------------------------------------------------------------------
  namespace  // anonymous
  {
    eterm NewEterm(ConstRefPPMonoidElem pp);
    TermList NewLPPTList(const vector<RingElem>& GB);
    RingElem NewPoly(PolyRing HSRing, unipoly f);
  }
  //----------------------------------------------------------------------

  RingElem HilbertNumQuot_C(const ideal& I)
  {
    VerboseLog VERBOSE("HilbertNumQuot_C");
    VERBOSE(90) << "original C code " << std::endl;

    PolyRing QQt = RingQQt(1);
    if (IsZero(I)) return one(QQt);
    const vector<RingElem>& GB = TidyGens(I);
    for (const auto& g: GB)
      if (!IsHomog(g))  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    ::StartPoincare((int)NumIndets(RingOf(I)));
    TermList TL = NewLPPTList(GB);
    unipoly PN = TLPoincareNumerator(TL); // TL is freed (?)
    RingElem ans = NewPoly(QQt, PN);
    FreeUnipoly(PN);
    return ans;
  }


  void EndPoincare_C()
  {
    ::EndPoincare(PoincareMaxPower); // clear *C* global variable
  }
  

  RingElem HilbertNumQuot(const ideal& I)
  {
    VerboseLog VERBOSE("HilbertNumQuot");
    VERBOSE(90) << "using CoCoALib polynomial arithmetic " << std::endl;
    
    if (!IsHomog(I))  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    const PolyRing QQt = RingQQt(1);
    if (IsZero(I)) return one(QQt);
    DenseUPolyRing HSRing = StartPoincareQQt((int)NumIndets(RingOf(I)));
    TermList TL = NewLPPTList(GBasis(I));
    RingElem PN = TLPoincareNumeratorCPP(HSRing, TL); // TL is freed (?)
    return PolyAlgebraHom(HSRing, QQt, indets(QQt)) (PN);
  }


  //  RingElem MGHilbertNumQuot(const SparsePolyRing& HSRing, const ideal& I)
  RingElem MGHilbertNumQuot(const ideal& I)
  {
    //  if (!HasPositiveGrading(RingOf(I))) CoCoA_THROW_ERROR1(ERR::ReqPositiveGrading);
    if (!IsHomog(I))  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    const SparsePolyRing RingOfI = RingOf(I);
    const SparsePolyRing QQt(RingQQt(GradingDim(RingOfI)));
    if (IsZero(I)) return one(QQt);
    ::StartPoincare((int)NumIndets(RingOfI));
    TermList TL = NewLPPTList(GBasis(I));
    RingElem PN = TLPoincareNumeratorCPP(QQt, PPM(RingOfI), TL); // TL is freed (?)
    return PN;
  }


  namespace  // anonymous
  {
    factorization<RingElem> DenFactors(long d)
    {
      SparsePolyRing QQt = RingQQt(1);
      if (d==0) return factorization<RingElem>(std::vector<RingElem>(0), std::vector<long>(0), one(QQt));
      return factorization<RingElem>(std::vector<RingElem>(1,one(QQt)-indet(QQt,0)), std::vector<long>(1,d), one(QQt));
    }

    factorization<RingElem> DenFactors(const SparsePolyRing& P)
    {
      if (IsStdGraded(P)) return DenFactors(NumIndets(P));
      long GrDim = GradingDim(P);
      SparsePolyRing QQt = RingQQt(GrDim);
      std::vector<RingElem> facs;
      std::vector<BigInt> v(GrDim);
      for (long i=0; i<NumIndets(P); ++i)
      {
        degree d(wdeg(indet(P,i)));
        for (long j=0; j<GradingDim(P); ++j)  v[j] = d[j];
        facs.push_back(one(QQt) + monomial(QQt, -1, PPMonoidElem(PPM(QQt),v)));
      }
      return factorization<RingElem>(facs, std::vector<long>(NumIndets(P),1), one(QQt));
    }

  } // anonymous namespace


  HPSeries HilbertSeriesQuot(const ideal& I)
  {
    const SparsePolyRing P(RingOf(I));
    if (IsStdGraded(P)) return HPSeries(HilbertNumQuot(I), DenFactors(P));
    if (!HasPositiveGrading(RingOf(I))) CoCoA_THROW_ERROR1(ERR::ReqPositiveGrading);
    return HPSeries(MGHilbertNumQuot(I), DenFactors(P));
  }


  HPSeries HSSimplified(const HPSeries& PSer)
  {
//     if (!IsStandard(PSer))  CoCoA_THROW_ERROR2(ERR::BadArg, "HPSeries must be standard");
    //    if IsZero(num(PSer)) return 0/1;
    if (DenFactors(PSer).myMultiplicities().empty()) return PSer;
    long HPDenExp = sum(DenFactors(PSer).myMultiplicities());
    PolyRing QQt = RingQQt(1);
    
    RingElem t = indet(QQt, 0);
    RingElem HPNum = num(PSer);
    while (IsDivisible(HPNum, HPNum, 1-t))       // HPNum /= (1-t);
      --HPDenExp;
    return HPSeries(HPNum, DenFactors(HPDenExp));
  }


  vector<BigInt> HVector(const HPSeries& HPS)
  {
    if (IsZero(num(HPS))) return  vector<BigInt>(0);
    HPSeries SimplHPS = HSSimplified(HPS);
    vector<BigInt> HV(deg(num(SimplHPS))+1);
    for (SparsePolyIter it=BeginIter(num(SimplHPS)); !IsEnded(it); ++it)
      HV[deg(PP(it))] = ConvertTo<BigInt>(coeff(it));
    return HV;
  }

  
  long dim(const HPSeries& HPS)
  {
    if (IsZero(num(HPS))) return sum(DenFactors(HPS).myMultiplicities());
    HPSeries SimplHPS = HSSimplified(HPS);
    if (DenFactors(SimplHPS).myMultiplicities().empty()) return 0;
    return sum(DenFactors(SimplHPS).myMultiplicities());
  }
  

  BigInt multiplicity(const HPSeries& HPS)
  {
    HPSeries SimplHPS = HSSimplified(HPS);
    vector<BigInt> HV = HVector(SimplHPS);
    return sum(HV);
  }
  

  long DimQuot(const ideal& I)
  {
    if (!IsStdGraded(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::BadRing, "must be standard graded");
    if (AreGensMonomial(I)) return dim(HilbertSeriesQuot(radical(I)));
    return dim(HilbertSeriesQuot(I));
  }
  

  BigInt MultiplicityQuot(const ideal& I)
  {
    if (!IsStdGraded(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::BadRing, "must be standard graded");
    return multiplicity(HilbertSeriesQuot(I));
  }

  
  RingElem HilbertPoly(const HPSeries& HPS)
  {
    // check in HilbertPoly(PModI)
    //     if (!IsStdGraded(P)) ...
    if (IsZero(num(HPS))) return num(HPS);
    HPSeries SimplHPS = HSSimplified(HPS);
    vector<BigInt> HV = HVector(SimplHPS);
    if (DenFactors(SimplHPS).myMultiplicities().empty())
      return zero(owner(num(SimplHPS)));
    long DIM = sum(DenFactors(SimplHPS).myMultiplicities());
    RingElem t = indet(RingQQt(1),0);  
    RingElem HP(RingQQt(1));
    for (long i=len(HV)-1; i>=0; --i) HP += HV[i]*binomial(t+DIM-i-1, DIM-1);
    return HP;
  }
  

  HPSeries HilbertSeries(const QuotientRing& PModI)
  { return HilbertSeriesQuot(DefiningIdeal(PModI)); }

  
  RingElem HilbertPoly(const QuotientRing& PModI)
  {
    if (!IsStdGraded(BaseRing(PModI)))
      CoCoA_THROW_ERROR2(ERR::BadRing, "must be standard graded");
    return HilbertPoly(HilbertSeries(PModI));
  }

  
  vector<BigInt> HVector(const QuotientRing& PModI)
  { return HVector(HilbertSeries(PModI)); }
    
  
  //----------------------------------------------------------------------
  namespace  // anonymous
  {
    eterm NewEterm(ConstRefPPMonoidElem pp)
    {
      PPMonoid PPM = owner(pp);
      eterm aux = eterm_init(NumIndets(PPM));
      ints OccInd = Indets(aux);
      const vector<long> exps = exponents(pp);
      for ( long i=0 ; i<NumIndets(PPM) ; ++i )
        if ( exps[i] != 0 )
        {
          eterm_put_nth(aux, (unsigned long)i+1, exps[i]); // range from 1
          IntsPutLast(OccInd, i+1);
        }
      return aux;
    }


    TermList NewLPPTList(const vector<RingElem>& GB)
    {
      SparsePolyRing P = owner(GB[0]);
      TermList aux = NewTList(len(GB), NumIndets(P));;
      int NewLen =0;

      for ( long i=0 ; i<len(GB) ; ++i )
      {
        eterm t = NewEterm(LPP(GB[i]));
        InsInTList(aux, t, &NewLen);
      }
      TListReduceSize (aux, NewLen);

      return aux;
    }


    RingElem NewPoly(PolyRing HSRing, unipoly f)
    {
      RingElem ans(HSRing);

      if ( HasMPZCoeffs(f) )
        for ( long d=0 ; d<=UPDeg(f) ; ++d )
          ans += BigIntFromMPZ(f[d].z) * IndetPower(HSRing,0,d);
      else
        for ( long d=0 ; d<=UPDeg(f) ; ++d )
          ans += f[d].i * IndetPower(HSRing,0,d);
      return ans;
    }

  } // anonymous namespace
  
// the following code is for using a Hilbert-server

//   void SendPPListToHilbert(std::ostream& out, const vector<RingElem>& GB)
//   {
//     SparsePolyRing P = owner(GB[0]);

//     // "1 P" stands for unused truncation,
//     out << len(GB) << NumIndets(P) << "1 P" << std::endl << std::endl;
//     for ( long i=0 ; i<len(GB) ; ++i )
//       out << LPP(GB[i]) << "," << std::endl;
//   }


//   RingElem ReceivePoincareNumerator(std::istream& in)
//   {
//     PolyRing P = NewPolyRing(NewRingZZ(), 1);
//     RingElem ans(P);

//     size_t d;
//     in >> d;
//     BigInt tmpcoeff;
//     for ( size_t i=0 ; i<=d ; ++i )
//     {
//       in >> tmpcoeff;
//       ans += tmpcoeff * IndetPower(P,0,i);
//     }
//     return ans;
//   }


//   RingElem TmpUniPoincareNumeratorMod(ideal I)
//   {
//     // (UsingSocket)
//     const unsigned short HilbertPort = 123456;
//     unique_ptr<SocketStream> SockPtr;
//     SockPtr.reset(new SocketStream(HilbertPort));
//     SetGlobalInput(*SockPtr);
//     SetGlobalOutput(*SockPtr);

//     const vector<RingElem>& GB = TidyGens(I);
//     SendPPListToHilbert(*SockPtr, GB);
//     RingElem ans = ReceivePoincareNumerator(*SockPtr);
//     GlobalLogput() << "poincare numerator = " << ans << std::endl;
//     return ans;
//   }


} // end of namespace CoCoA
