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
    //  if (!HasPositiveGrading(RingOf(I))) CoCoA_THROW_ERROR(ERR::NotPositiveGrading, "MGHilbertNumQuot");
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
//     if (!IsStandard(PSer))
//       CoCoA_THROW_ERROR("HPSeries must be standard");
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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-hilbert.C,v 1.20 2024/07/25 14:58:46 bigatti Exp $
// $Log: SparsePolyOps-hilbert.C,v $
// Revision 1.20  2024/07/25 14:58:46  bigatti
// Summary: using CoCoA_THROW_ERROR1/2
//
// Revision 1.19  2024/07/02 15:36:15  bigatti
// Summary:Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.18  2024/03/15 18:57:32  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.17  2024/02/05 12:35:36  bigatti
// Summary: radical_tmp  -->  radical
//
// Revision 1.16  2023/07/06 16:34:06  bigatti
// Summary: added verbosity
//
// Revision 1.15  2023/03/23 04:44:20  bigatti
// Summary: fixed HilbertPoly of 0-dim ideals Redmine #1729
// also added check for standard grading (otherwise there is no single HilbertPoly)
//
// Revision 1.14  2023/03/09 22:41:07  abbott
// Summary: Tidied use of exponents (redmine 1694)
//
// Revision 1.13  2022/10/11 15:44:00  bigatti
// Summary: using new interface of exponents returning vector
//
// Revision 1.12  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.11  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.10  2020/02/13 14:00:07  bigatti
// -- now using auto in for loop
//
// Revision 1.9  2019/10/04 09:42:15  bigatti
// Summary: fixed bug linux/mac with different evaluation order
//
// Revision 1.8  2019/10/03 13:33:56  bigatti
// -- implemented radical for monomial ideals (and used where useful)
//
// Revision 1.7  2019/10/02 10:38:47  bigatti
// -- fixed multivariate HilbertSeriesQuot: check for HasPositiveGrading
//
// Revision 1.6  2019/03/04 11:16:39  abbott
// Summary: Commented out apparently useless include (for auto_ptr/unique_ptr)
//
// Revision 1.5  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.4  2018/05/17 15:48:33  bigatti
// -- added include SparsePolyIter
// -- renamed VectorOperations --> VectorOps
// -- sorted #includes
//
// Revision 1.3  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.2  2018/04/10 15:45:20  bigatti
// -- fixed includes
//
// Revision 1.1  2018/04/06 15:15:36  bigatti
// -- renamed TmpHilbert.C
//
// Revision 1.38  2018/03/15 14:18:06  bigatti
// -- added files SparsePolyOps-ideal.H and SparsePolyOps-involutive.H
//
// Revision 1.37  2017/04/13 15:53:33  bigatti
// -- updated Copyright
//
// Revision 1.36  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.35  2016/10/27 13:20:07  bigatti
// -- added MultiplicityQuot, DimQuot
//
// Revision 1.34  2016/02/16 16:48:32  bigatti
// Summary: undid check for IsPositiveGrading (needs to be cleaned)
//
// Revision 1.33  2016/02/16 16:29:16  bigatti
// -- added check for IsPositiveGrading
// -- added commented out code for dim and multiplicity
//
// Revision 1.32  2015/11/04 12:14:21  abbott
// Summary: Added a const
//
// Revision 1.31  2015/09/03 10:18:41  bigatti
// -- changes by Christof Soeger (in Aarhus)
//
// Revision 1.30  2014/11/11 09:23:20  bigatti
// -- fixes and improvements
//
// Revision 1.29  2014/10/24 12:47:22  bigatti
// -- added HilbertPoly, HVector (to be refined)
//
// Revision 1.28  2014/07/30 14:11:06  abbott
// Summary: Changed name AmbientRing --> RingOf
// Author: JAA
//
// Revision 1.27  2014/07/07 13:08:40  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.26  2014/06/17 10:15:01  abbott
// Summary: Removed pointless call to AsPolyRing (on RingQQt)
// Author: JAA
//
// Revision 1.25  2014/04/30 16:16:24  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.24  2014/03/24 12:09:21  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.23  2014/01/29 18:45:20  bigatti
// -- minor changes to HPSeries ctor (size --> len,  SparsePolyRing --> PolyRing)
//
// Revision 1.22  2014/01/29 17:32:50  bigatti
// -- added HPSeries ctor for DenseRepr (by Christof Soeger)
//
// Revision 1.21  2013/07/30 16:45:31  bigatti
// -- fixed HilbertNumQuot
//
// Revision 1.20  2013/07/30 15:04:03  bigatti
// -- added HPSeries class
// -- added HilbertSeriesQuot
// -- simplified arguments: all results in RingQQt
//
// Revision 1.19  2013/06/18 12:27:36  abbott
// Renamed HibertSeriesPolyRing to RingQQt.
//
// Revision 1.18  2013/06/17 08:54:38  abbott
// Added HilbertSeriesPolyRing (untested).
//
// Revision 1.17  2013/02/04 17:33:26  bigatti
// -- only one poincare_init for leak control (but useless unipoly for
//    some cases)
//
// Revision 1.16  2013/01/31 10:18:03  bigatti
// -- improve comment
//
// Revision 1.15  2012/02/08 16:15:18  bigatti
// -- just a comment fix
//
// Revision 1.14  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.13  2011/05/09 13:41:53  bigatti
// -- moved conversion code into anonymous namespace
//
// Revision 1.12  2011/04/26 10:12:21  bigatti
// -- univariate hilbert creates its own (ZZ[t]) ring
//
// Revision 1.11  2011/04/12 13:35:07  bigatti
// -- fixed for input ideal(0)
//
// Revision 1.10  2011/04/08 14:04:45  bigatti
// -- renamed HilbertNumeratorMod into HilbertNumQuot
//
// Revision 1.9  2011/03/10 18:02:07  bigatti
// -- added cast to unsigned long for call to old eterm_put_nth
//
// Revision 1.8  2011/03/10 17:57:29  bigatti
// -- using long instead of size_t
// -- using len instead of size()
//
// Revision 1.7  2010/10/29 09:40:36  bigatti
// -- Globals for C++ Poincare are now in GlobalManager
// -- Globals for C Poincare are to be freed manually with EndPoincare_C
//
// Revision 1.6  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.5  2009/10/26 15:45:10  bigatti
// -- added CopyFromMPZ in ZZ ctor
//
// Revision 1.4  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.3  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/10/10 14:36:34  bigatti
// new: TmpHilbertDir/TmpPoincareCPP.[CH]: poincare code using
// 	C++ univariate polynomials (DenseUPolyRing)
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.3  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.2  2007/01/17 17:54:24  bigatti
// -- fixed upper case in includes names
//
// Revision 1.1  2007/01/17 17:38:11  bigatti
// -- moved all cocoa-4 code for hilbert into src/TmpHilbertDir
//
// Revision 1.6  2006/11/27 15:30:04  cocoa
// -- fixed: minor memory leak in Hilbert
//
// Revision 1.5  2006/11/24 17:12:05  cocoa
// -- reorganized includes of header files
//
// Revision 1.4  2006/11/20 15:17:23  cocoa
// -- added check for homogeneity
//
// Revision 1.3  2006/11/17 18:14:01  cocoa
// -- fixed: now Hilbert computes TidyGens instead of cheating
// -- myGBasis is much more efficient on monomial input
//
// Revision 1.2  2006/11/16 17:39:51  cocoa
// -- fixed some calls to rum: lots of code could be deleted, but this is
//    all much more fragile than I would like
//
// Revision 1.1  2006/10/09 16:48:58  cocoa
// -- first import
//
