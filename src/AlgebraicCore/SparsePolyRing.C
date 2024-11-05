//   Copyright (c)  2005-2018  John Abbott and Anna M. Bigatti
//   Authors:  2005-2018  John Abbott and Anna M. Bigatti

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


// Source code for abstract class SparsePolyRing and friends

#include "CoCoA/SparsePolyRing.H"

#include "CoCoA/MatrixOps.H"             // for rk
#include "CoCoA/MatrixForOrdering.H"     // for IsPositiveGrading
#include "CoCoA/OpenMath.H"
#include "CoCoA/RingDistrMPolyClean.H"   // for NewPolyRing_DMP
#include "CoCoA/RingDistrMPolyInlFpPP.H" // for NewPolyRing_DMPII
#include "CoCoA/RingDistrMPolyInlPP.H"   // for NewPolyRing_DMPI
#include "CoCoA/RingFp.H"                // for IsRingFp
#include "CoCoA/matrix.H"                // for ConstMatrixView ????


#include <algorithm>
using std::sort;    // for AreGoodIndetNames
#include <iostream>
// using std::ostream in SparsePolyRingBase::myOutput
//#include <vector>
using std::vector;


namespace CoCoA
{

//   namespace // anonymous
//   {
//     // This fn is needed in a call to std::transform
//     CoeffPP CoeffPPCtor(const pair<PPMonoidElem, RingElem>& arg)
//     {
//       return CoeffPP(arg.second, arg.first);
//     }
//   } // end of namespace anonymous


  SparsePolyRing NewPolyRing(const ring& CoeffRing, const PPMonoid& PPM)
  {
//     if (IsPPMonoidOv(PPM))
//     {
//       if (IsRingFp(CoeffRing))
//         CoCoA_THROW_ERROR2(ERR::NYI, "NewPolyRing DMPII with PPM");
//       return NewPolyRing_DMPI(CoeffRing, PPM);
//     }
    return NewPolyRing_DMP(CoeffRing, PPM);// the only one for DMP
    // !!! ANNA: but if PPM is PPMonoidOv we could be clever!!!
  }

  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetSyms, const PPOrdering& ord)
  {
    //    std::cout << "------NewPolyRing ord" << std::endl;
    if (IsRingFp(CoeffRing))
      return NewPolyRing_DMPII(CoeffRing, IndetSyms, ord);
    return NewPolyRing_DMPI(CoeffRing, IndetSyms, ord);
  }

  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetSyms, const PPOrderingCtor& OrdCtor)
  {
    //    std::cout << "------NewPolyRing OrdCtor" << std::endl;
    if (IsRingFp(CoeffRing))
      return NewPolyRing_DMPII(CoeffRing, IndetSyms, OrdCtor);
    return NewPolyRing_DMPI(CoeffRing, IndetSyms, OrdCtor);
  }


  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetSyms)
  {
    return NewPolyRing(CoeffRing, IndetSyms, StdDegRevLex);
  }


  SparsePolyRing NewPolyRingWeights(const ring& CoeffRing, const std::vector<symbol>& IndetSyms, ConstMatrixView Ws)
  {
    if (rk(Ws)!=NumRows(Ws))  CoCoA_THROW_ERROR1(ERR::ReqFullRank);
    for (long i=0; i < NumRows(Ws); ++i)
      for (long j=0; j < NumCols(Ws); ++j)
        if (sign(Ws(i,j)) < 0)
          CoCoA_THROW_ERROR2(ERR::NYI, "Temporarily requiring weights to be non-negative");
    PPOrdering ord = NewMatrixOrdering(MakeTermOrdMat(Ws), NumRows(Ws));
    if (IsRingFp(CoeffRing))
      return NewPolyRing_DMPII(CoeffRing, IndetSyms, ord);
    return NewPolyRing_DMPI(CoeffRing, IndetSyms, ord);
  }


  bool AreGoodIndetNames(const ring& CoeffRing, const std::vector<symbol>& IndetNames)
  {
    // inefficient: we might know that SymbolsCR and IndetNames are good
    // does it matter?
    vector<symbol> syms = symbols(CoeffRing);
    syms.insert(syms.end(), IndetNames.begin(), IndetNames.end());
    sort(syms.begin(), syms.end());
    const long NumSyms = len(syms);
    for (long i=0; i < NumSyms-1; ++i)
    {
      if (syms[i] == syms[i+1]) return false;
      if (head(syms[i]) == head(syms[i+1]) &&
          NumSubscripts(syms[i]) != NumSubscripts(syms[i+1]))
        return false;
    }
    return true;
  }


  ConstMatrixView OrdMat(const SparsePolyRing& Rx)
  { return OrdMat(PPM(Rx)); }

  ConstMatrixView GradingMat(const SparsePolyRing& Rx)
  { return GradingMat(PPM(Rx)); }

  bool HasPositiveGrading(const SparsePolyRing& Rx)
  { return IsPositiveGrading(GradingMat(Rx)); }


//---- Functions for RingElem moved to SparsePolyOps-RingElem.C

  /*----------------------------------------------------------------------
    Member functions inherited from ring with a single implementation
    for all SparsePolyRing implementations
    ----------------------------------------------------------------------*/

  void SparsePolyRingBase::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("polyd", "poly_ring_d");
    OMOut << myCoeffRing();
    OMOut << myNumIndets(); //???? losing the ordering and grading info here!!!
    OMOut->mySendApplyEnd();
  }


  void SparsePolyRingBase::mySymbols(std::vector<symbol>& SymList) const
  {
    myCoeffRing()->mySymbols(SymList);
//    myPPM()->mySymbols(SymList);
    const vector<symbol>& PPs = symbols(myPPM());
    SymList.insert(SymList.end(), PPs.begin(), PPs.end());
  }


//   void SparsePolyRingBase::mySymbolsAndValues(std::vector<symbol>& SymList, std::vector<RingElem>& v) const
//   {
//     std::vector<RingElem>& TmpV;
//     myCoeffRing()->mySymbolsAndValues(SymList, TmpV);
//     myPPM()->mySymbols(SymList);
//     v.clear();
//     for (vector<RingElem>::const_iterator it=TmpV.begin(); it!=TmpV.end() ; ++it)
//       v.push_back(myCoeffEmbeddingHomCtor()(*it));
//     for (vector<PPMonoidElem>::const_iterator it=indets(myPPM()).begin(); it!=indets(myPPM()).end() ; ++it)
//       v.push_back(myMonomial(1, raw(*it)));
//   }


//   void SparsePolyRingBase::myOutput(std::ostream& out, ConstRawPtr rawf) const
//   {
//     if (myIsZero(rawf)) { out << '0'; return; }


//     const ring& R = myCoeffRing();
//     const PPMonoid PPM = myPPM();
//     const bool IsQuotientOfZZ = IsQuotientRing(R) && IsZZ(BaseRing(R));
//     //const bool IsQuotientOfZ = IsRingFp(R);
//     const bool IsNumberRing = IsZZ(R) || IsQuotientOfZZ || IsRingTwinFloat(R); // || IsQQ(R) 

//     bool IsFirstCoeff = true;
//     for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf) ; ++itf, IsFirstCoeff = false)
//     {
//       bool PrintStar = true;
//       RingElemAlias c = coeff(itf);
//       ConstRefPPMonoidElem pp = PP(itf);
//       const bool IsWithMinus = R->myIsPrintedWithMinus(raw(c));

//       // ---- coefficient ----
//       if (!IsFirstCoeff)
//       {
//         out << ' ';
//         if (!IsWithMinus)  out << '+';
//       }
//       if (IsOne(pp)) { out << c; continue; }
//       // ---------- PP != 1 ----------
//       if (IsOne(c))
//       { // Do not print "1 * "...
//         PrintStar = false;
//         goto PrintPP;
//       }
//       if ( IsWithMinus && IsMinusOne(c) ) // in some Z/(n) "-1" prints as n-1
//       { // Do not print "-1 * "...
//         out << '-';
//         PrintStar = false;
//         goto PrintPP;
//       }
//       // General case: coeff is neither +1 nor -1
//       if (IsNumberRing || R->myIsPrintAtom(raw(c)) ||
//           (IsWithMinus && R->myIsPrintAtom(raw(-c))) )
//       {
//         out << c;
//         goto PrintPP;
//       }
//       if (!IsFirstCoeff && IsWithMinus) out << '+'; // not printed before
//       out << '(' << c << ')';

//     PrintPP:      // ---- PP ----
//       if (PrintStar)  out << '*';
//       out << pp;
//     }
//   }




//   //-- HomImpl ----------------------------------------

//   //-- CoeffEmbeddingHomImpl ----------------------------------------

} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyRing.C,v 1.230 2024/08/28 15:24:55 bigatti Exp $
// $Log: SparsePolyRing.C,v $
// Revision 1.230  2024/08/28 15:24:55  bigatti
// Summary: now using CoCoA_THROW_ERROR1/2
//
// Revision 1.229  2024/03/15 16:00:15  bigatti
// Summary: renamed  NewPolyRing(K, X, Ws)  into NewPolyRingWeights(K, X, Ws)
//
// Revision 1.228  2024/03/14 11:04:08  bigatti
// Summary: added check for negative entries in NewPolyRing(K,X,Weights)
//
// Revision 1.227  2024/02/05 16:20:01  bigatti
// Summary: added check for full rank in NewPolyRing with weights
//
// Revision 1.226  2023/07/06 16:32:34  bigatti
// Summary: added NewPolyRing with weight matrix
//
// Revision 1.225  2023/03/13 19:59:29  abbott
// Summary: Removed useless include directives
//
// Revision 1.224  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.223  2022/02/08 20:18:55  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.222  2021/01/12 13:27:59  abbott
// Summary: Removed many superfluous includes
//
// Revision 1.221  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.220  2020/01/26 14:37:51  abbott
// Summary: Removed useless include (of NumTheory)
//
// Revision 1.219  2019/09/25 14:28:04  bigatti
// -- added HasPositiveGrading
// -- removed old code for RingElem (now is in SparsePolyOps-RingElem)
//
// Revision 1.218  2018/10/02 09:45:30  abbott
// Summary: Moved pseudo-ctors NewPolyRing(CoeffRing, NumIndets, ...) to obsolescent
//
// Revision 1.217  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.216  2018/06/15 08:48:16  abbott
// Summary: Added NewPolyRing with args CoeffRing, NumIndets, PPOrd
//
// Revision 1.215  2018/05/18 16:35:52  bigatti
// -- split SparsePolyOps-RingElem from SparsePolyRing
//
// Revision 1.214  2018/05/18 12:23:50  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.213  2018/05/17 15:51:20  bigatti
// -- added include SparsePolyIter
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.212  2018/05/16 14:35:58  abbott
// Summary: Reversed order of loop in RandomLinearForm
//
// Revision 1.211  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.210  2018/04/19 13:40:03  bigatti
// -- minor change in CoeffVecWRT (using indices instead of x)
//
// Revision 1.209  2018/04/16 21:50:29  bigatti
// -- moved out code on RingHom
// -- improved RandomLinearForm
//
// Revision 1.208  2018/04/09 16:36:09  bigatti
// -- fixed includes
//
// Revision 1.207  2018/04/06 17:41:37  bigatti
// -- fixed includes
//
// Revision 1.206  2018/03/20 11:43:31  bigatti
// -- added RandomLinearForm
// -- removed commented code for ideals
//
// Revision 1.205  2018/03/15 14:38:56  bigatti
// -- new file SparsePolyOps-ideal.H
//
// Revision 1.204  2018/03/13 14:24:23  bigatti
// -- fixed CRTPoly (swapped arguments)
//
// Revision 1.203  2018/02/27 16:36:49  bigatti
// -- split SparsePolyRing.C part on ideals into SparsePolyRing_ideal.C
//
// Revision 1.202  2018/02/22 17:07:12  abbott
// Summary: Minor tidying
//
// Revision 1.201  2018/02/21 10:55:11  abbott
// Summary: Updated ctor call for RatReconstructByContFrac
//
// Revision 1.200  2018/01/17 10:31:49  abbott
// Summary: Added PowerOverflowCheck to mySquare
//
// Revision 1.199  2017/12/18 12:07:26  abbott
// Summary: Correct handling for polys with coeffs in rings with zero-divisors
//
// Revision 1.198  2017/11/29 17:39:20  bigatti
// -- added GBasisRealSolveCore
// -- modified GBasisSelfSatCore: now DOES NOT store GBasis
//
// Revision 1.197  2017/11/24 17:46:39  bigatti
// -- renamed GBasisSelfSat --> GBasisSelfSatCore
// -- added GBasisSelfSat in cpkg5
//
// Revision 1.196  2017/11/23 12:33:01  bigatti
// -- added GBasisSelfSat (still buggy)
//
// Revision 1.195  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.194  2017/11/06 15:37:44  bigatti
// -- added GBasisSelfSat (hidden, because buggy)
//
// Revision 1.193  2017/09/25 12:38:12  abbott
// Summary: Added QuotientBasisSorted
//
// Revision 1.192  2017/09/06 14:09:39  abbott
// Summary: Changed ERR::SERIOUS into ERR:ShouldNeverGetHere
//
// Revision 1.191  2017/07/23 15:31:33  abbott
// Summary: Added GBasisTimeout (just for ideals)
//
// Revision 1.190  2017/07/14 09:35:37  bigatti
// -- just some comments for debugging ring construction
//
// Revision 1.189  2017/07/04 14:51:22  abbott
// Summary: CutLf now gives error is arg is zero
//
// Revision 1.188  2017/07/03 19:53:38  abbott
// Summary: Added CutLF; minor cleaning to LF
//
// Revision 1.187  2017/05/11 08:46:46  bigatti
// -- added error ERR::CannotReconstruct
// -- cleaned up code accordingly
//
// Revision 1.186  2017/04/27 16:17:13  bigatti
// -- fixed IdealOfGBasis(ideal(zero(P)))
//
// Revision 1.185  2017/04/26 15:56:00  bigatti
// -- added IdealOfGBasis, IdealOfMinGens
//
// Revision 1.184  2017/04/18 09:53:25  bigatti
// -- removed two comments in NewPolyRing
//
// Revision 1.183  2017/04/07 14:18:21  bigatti
// -- improved LF (PushBack & CmpWDeg)
// -- fixed verbosity for IsMaximalSPR
//
// Revision 1.182  2017/03/08 15:33:55  bigatti
// -- added ReducedGBasis
//
// Revision 1.181  2017/02/14 17:32:49  bigatti
// -- minor changes in IamMaximalSPRTest_aux
//
// Revision 1.180  2017/02/06 16:14:40  bigatti
// -- IsMaximalSPR now is as in paper
//
// Revision 1.179  2016/12/15 14:39:00  bigatti
// -- maximality: moved check Frobenius after indeterminates
//
// Revision 1.178  2016/12/06 16:11:31  bigatti
// -- improved maximal test (as in paper)
//
// Revision 1.177  2016/12/02 06:38:15  bigatti
// -- added verbosity to IsMaximalSPR
//
// Revision 1.176  2016/11/25 15:46:21  abbott
// Summary: Fixed bug redmine #981
//
// Revision 1.175  2016/11/18 07:37:30  bigatti
// -- just comments
//
// Revision 1.174  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.173  2016/11/07 13:54:21  bigatti
// ++ added AreGensSqFreeMonomial
// ++ replaced all SquareFree into SqFree
//
// Revision 1.172  2016/11/07 12:11:21  bigatti
// -- added HasGBasis / IhaveGBasis
//
// Revision 1.171  2016/11/03 12:25:25  abbott
// Summary: Changed IsRadical (for PPMonoidElem) into IsSqFree
//
// Revision 1.170  2016/09/08 14:12:31  bigatti
// -- mySetMaximalFlag --> myAssignMaximalFlag
// -- mySetPrimeFlag --> myAssignPrimeFlag
// -- updated the related code
// -- (still "old design", but better aligned)
//
// Revision 1.169  2016/06/24 14:27:41  bigatti
// -- renamed CRT_poly --> CRTPoly
//
// Revision 1.168  2016/06/20 15:27:07  bigatti
// -- renamed MinPolyXX --> MinPolyQuotXX
// -- added default MinPolyQuot
//
// Revision 1.167  2016/06/14 08:26:23  bigatti
// -- commented out debugging printouts
//
// Revision 1.166  2016/06/10 15:49:59  bigatti
// -- added IsRadicalSPR, IsMaximalSPR (and then updated myMaximalTest)
//
// Revision 1.165  2016/04/27 07:01:55  bigatti
// -- new: myReset()  (and using it in all functions modifying ideal)
// -- removed: myClearGBasis()
//
// Revision 1.164  2015/12/08 13:56:09  abbott
// Summary: Updated Mario's code!  Very many changes!
//
// Revision 1.163  2015/12/04 15:22:24  bigatti
// -- renamed ComputeSSaturation into ComputeSaturation
//
// Revision 1.162  2015/12/01 13:11:01  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.161  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.160  2015/10/01 10:14:24  bigatti
// -- added CRT_poly, RatReconstructPoly
//
// Revision 1.159  2015/06/30 12:51:10  abbott
// Summary: Added new fn IndetsIn
// Author: JAA
//
// Revision 1.158  2015/06/11 16:55:16  bigatti
// -- added monomial(ring, pp) and monomial(ring, expv)
//
// Revision 1.157  2015/05/20 12:59:31  bigatti
// -- renamed ComputeCColon --> ComputeColon
//
// Revision 1.156  2015/04/30 08:50:10  bigatti
// -- added myIsZeroDivisor
// -- some fixes in myDivMod
//
// Revision 1.155  2015/04/24 15:40:58  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.154  2015/04/16 16:36:33  abbott
// Summary: Cleaned impls of myPowerOverflowCheck
// Author: JAA
//
// Revision 1.153  2015/04/13 14:42:58  abbott
// Summary: Added call to PowerOverflowCheck
// Author: JAA
//
// Revision 1.152  2015/03/09 12:12:32  bigatti
// -- just two empty lines
//
// Revision 1.151  2014/12/04 08:56:26  bigatti
// -- added check for empty reducer list in NR
//
// Revision 1.150  2014/07/31 13:10:45  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.149  2014/07/30 14:09:59  abbott
// Summary: Changed name AmbientRing --> RingOf
// Author: JAA
//
// Revision 1.148  2014/07/28 15:52:10  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble); added general defn of myCoeffEmbeddingHomCtor
// Author: JAA
//
// Revision 1.147  2014/07/11 15:37:16  bigatti
// -- default implementation of myOutputSelf(OpenMathOutput& OMOut)
//
// Revision 1.146  2014/07/09 11:42:41  abbott
// Summary: Removed some cruft
// Author: JAA
//
// Revision 1.145  2014/07/08 15:23:34  abbott
// Summary: Updated comment
// Author: JAA
//
// Revision 1.144  2014/07/08 13:14:41  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.143  2014/07/08 09:02:09  abbott
// Summary: Corrected silly typo
// Author: JAA
//
// Revision 1.142  2014/07/08 08:38:39  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.141  2014/07/07 17:11:57  abbott
// Summary: [MAJOR CHANGE] Removed AsSparsePolyRing; added SparsePolyRingPtr
// Author: JAA
//
// Revision 1.140  2014/07/04 12:55:30  abbott
// Summary: Revised following redefinition of symbols for a PPM
// Author: JAA
//
// Revision 1.139  2014/05/16 13:24:42  bigatti
// -- fixed problem with elim(.., ideal(0)) (?)
//
// Revision 1.138  2014/05/06 16:01:26  abbott
// Summary: Commented out unused "suspect" function MonomialMulTrick
// Author: JAA
//
// Revision 1.137  2014/05/06 13:20:41  abbott
// Summary: Changed names (my)MaxExponents into (my)Deg
// Author: JAA
//
// Revision 1.136  2014/04/30 16:13:56  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.135  2014/04/22 12:45:16  abbott
// Summary: Added commente out "general soln" for myDivMod
// Author: JAA
//
// Revision 1.134  2014/04/15 13:28:55  abbott
// Summary: Added new ClearDenom fn (with 2 args)
// Author: JAA
//
// Revision 1.133  2014/04/11 15:07:39  abbott
// Summary: Renamed TmpFactor to factor in include
// Author: JAA
//
// Revision 1.132  2014/03/27 14:57:23  bigatti
// -- added myMinimalize
//
// Revision 1.131  2014/03/21 13:07:51  bigatti
// -- added check in myIntersect
// -- cosmetics
//
// Revision 1.130  2014/02/26 15:35:57  abbott
// Summary: Changed some 1 char strings into chars
// Author: JAA
//
// Revision 1.129  2014/01/28 13:29:25  bigatti
// -- moved back check for IsInvertible in gcd (so that gcd(inv,0) --> 1)
//
// Revision 1.128  2014/01/28 13:09:06  bigatti
// -- added LF for ideal
//
// Revision 1.127  2014/01/28 09:47:56  abbott
// Improved impls of PushFront & PushBack (avoids double creation of PP).
//
// Revision 1.126  2013/10/28 13:21:58  bigatti
// -- fixed error message
//
// Revision 1.125  2013/06/28 17:03:51  abbott
// Modified semantics of IdealBase::myDivMod;
// it now returns a boolean.
// Several consequential changes.
//
// Revision 1.124  2013/06/27 16:56:09  abbott
// Changed arg checking fn to expect a "const char*" const instead of a "std::string"
// (as the ctor cost for a string is paid whenever the checking fn is called, even if
// there is no error to report).
//
// Revision 1.123  2013/06/12 08:49:17  bigatti
// -- added MinGens
//
// Revision 1.122  2013/05/30 13:13:58  bigatti
// -- added SetGBasisAsGens
//
// Revision 1.121  2013/05/30 07:24:10  bigatti
// -- improved myPrimeTest (waiting for IsMaximal3)
//
// Revision 1.120  2013/05/28 13:29:54  abbott
// Added printing for struct CoeffPP.
//
// Revision 1.119  2013/05/27 17:16:37  abbott
// Added missing const.
//
// Revision 1.118  2013/05/27 16:25:56  abbott
// Fixed "typo" (forgot to rename a couple of variables).
// Renamed args and local variables in myGCD.
//
// Revision 1.117  2013/05/27 15:00:11  abbott
// Added some arg checks to CoefficientsWRT and ContentWRT.
//
// Revision 1.116  2013/05/27 07:07:49  bigatti
// -- fixed all cases with 0's
//
// Revision 1.115  2013/05/24 17:55:19  bigatti
// -- fixed myColon for input with 0
//
// Revision 1.114  2013/03/25 17:04:19  abbott
// Major clean-up of interface to SmallFpImpl/SmallFpLogImpl/SmallFpDoubleImpl
// (underlying impl remains much the same).  Removed lots of cruft.
// Consequential changes to RingFp* classes; small change to SparsePolyRing.
//
// Revision 1.113  2013/02/21 17:35:00  bigatti
// -- simplified gcd (call to SyzOfGens instead of ComputeSyz)
//
// Revision 1.112  2013/02/21 12:51:42  abbott
// Added new fn UnivariateIndetIndex.
//
// Revision 1.111  2013/02/12 16:29:15  bigatti
// -- added ERR::ZeroGradingDim in IsHomog
//
// Revision 1.110  2013/01/18 18:16:02  abbott
// Modified CoefficientsWRT so that it now gives result guaranteeing that the entries are in decsreasing order of PP.
//
// Revision 1.109  2013/01/17 15:11:19  abbott
// Added new fn CoeffVecWRT.
// Added new virt mem fns myImageLiesInSubfield & IamPartial.
//
// Revision 1.108  2012/10/24 13:38:52  abbott
// Changed ConstRefRingElem into RingElemAlias in ctor calls and some local variable types.
//
// Revision 1.107  2012/10/17 09:40:15  abbott
// Replaced  RefRingElem  by  RingElem&
// (plus a few consequential changes)
//
// Revision 1.106  2012/10/05 10:21:39  bigatti
// -- added LF (leading form)
//
// Revision 1.105  2012/10/03 12:28:46  bigatti
// -- fixed error message (homog)
//
// Revision 1.104  2012/10/02 16:44:04  bigatti
// -- added homog for ideal
//
// Revision 1.103  2012/05/30 13:44:11  bigatti
// -- renamed IhaveMonomialGensB3Value --> IhaveMonomialGens3Value
//
// Revision 1.102  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.101  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.100  2012/05/24 14:49:23  bigatti
// -- changed symbol "index" into "subscripts"
//
// Revision 1.99  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.98  2012/04/27 15:06:48  abbott
// Corrected myMaximalTest
//
// Revision 1.97  2012/04/15 20:10:13  abbott
// Added quick hacks to use old "faster" code for polynomial GCD in some cases.
//
// Revision 1.96  2012/04/13 16:27:51  abbott
// Added a working simple case in myDivMod; enough for Laura to proceed.
//
// Revision 1.95  2012/04/12 12:22:40  abbott
// Corrected handling of 0 in myNormalizeFracNoGcd;
// previously  0*x/(-1)  produced the result 0/(-1).
//
// Revision 1.94  2012/04/11 10:47:20  abbott
// Added special handling for the case 0/1 in myNormalizeFracNoGcd
// (without it the code wrongly divides by zero).
//
// Revision 1.93  2012/04/04 08:47:53  bigatti
// -- improved error messages in myGcd
//
// Revision 1.92  2012/04/03 16:12:46  abbott
// Changed slightly the fn signatures of CoefficientsWRT.
// Added some things to *.C file -- will complete later.
//
// Revision 1.91  2012/03/02 14:21:46  bigatti
// -- added ContentWRT, CoefficientsWRT(f, x)
//
// Revision 1.90  2012/02/24 13:10:16  abbott
// Removed cruft (commented out useless code).
//
// Revision 1.89  2012/02/14 15:16:14  bigatti
// -- commented out unused code
//
// Revision 1.88  2012/02/10 17:09:38  abbott
// Added new fns  indets, CoefficientsWRT, ContentWRT.
//
// Revision 1.87  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.86  2012/02/08 15:14:13  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.85  2012/01/26 16:47:16  bigatti
// -- changed back_inserter into insert
//
// Revision 1.84  2011/12/05 16:32:11  bigatti
// -- fixed bug about saturation (by non-principal ideal)
//
// Revision 1.83  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.82  2011/11/07 10:55:23  bigatti
// -- AreMonomials is now public
//
// Revision 1.81  2011/08/24 10:29:55  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.80  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.79  2011/08/12 15:22:49  abbott
// Improved a comment.
//
// Revision 1.78  2011/07/07 10:03:54  abbott
// Removed "inline" before ourGetPtr (as compilation fails with -O2).
//
// Revision 1.77  2011/07/05 15:02:17  bigatti
// -- added AlexanderDual
// -- added ad-hoc functions for colon, elim on monomial ideals
//
// Revision 1.76  2011/06/27 13:30:37  bigatti
// -- new file for monomial ideals:
// -- some functions moved there, but some have to be declared in SparsePolyRing.H
//
// Revision 1.75  2011/06/23 16:04:46  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.74  2011/05/26 16:33:47  bigatti
// -- improved implementation of IamZero()
// -- check for IamZero before computing GBasis
//
// Revision 1.73  2011/05/20 16:12:59  bigatti
// -- added QuotientBasis
//
// Revision 1.72  2011/04/27 08:21:48  bigatti
// -- added gcd with coefficients in GCDDomain
//
// Revision 1.71  2011/04/19 13:59:49  bigatti
// -- added AreGensMonomial
//
// Revision 1.70  2011/04/12 09:52:54  bigatti
// -- added IsHomog(ideal), LT(ideal)
//
// Revision 1.69  2011/03/22 15:28:20  bigatti
// -- fixed IsIndetPosPower
//
// Revision 1.68  2011/03/16 15:40:22  bigatti
// -- added myIsIndetPosPower(f), IsIndetPosPower(f)
//
// Revision 1.67  2011/03/16 13:21:23  abbott
// Added comments for myIsHomogPartial & myCmpWDegPartial.
// Cleaned up impls of myIsHomog & myIsHomogPartial.
// Corrected typo in impl of myElim.
//
// Revision 1.66  2011/03/11 10:54:03  bigatti
// -- added mySaturate
//
// Revision 1.65  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.64  2011/03/01 14:10:47  bigatti
// -- added ClearDenom/myClearDenom
//
// Revision 1.63  2011/02/28 14:16:08  bigatti
// -- added GBasis(ideal)  -- only for SparsePolyRing
// -- error for myGcd when CoeffRing is not a field
//
// Revision 1.62  2011/02/23 16:04:50  bigatti
// -- more fixes to get gcd==1 when coprime
//
// Revision 1.61  2011/02/23 15:02:22  bigatti
// -- gcd for polys now returns 1 when coprime (i.e. not any invertible constant)
//
// Revision 1.60  2011/02/04 11:05:06  bigatti
// -- changed printing for coefficients in QQ: added parentheses
//
// Revision 1.59  2011/01/28 17:58:07  bigatti
// -- added myElim
//
// Revision 1.58  2011/01/28 11:41:42  bigatti
// -- added IsPrintedWithMinus
// -- improved myOutput
// -- fixed bug in IsMinusOne
//
// Revision 1.57  2011/01/18 14:38:43  bigatti
// -- moved **_forC5 functions into CoCoA-5/CoCoALibSupplement:
//    myMonomials_forC5, mySupport_forC5, monomials_forC5, support_forC5,
//    LPP_forC5, LT_forC5, LM_forC5
//
// Revision 1.56  2010/11/30 11:28:27  bigatti
// -- moved IndetsCalled into unique implementation in PolyRing
//
// Revision 1.55  2010/11/25 12:31:22  bigatti
// -- added myIndetsCalled
//
// Revision 1.54  2010/11/05 15:59:31  bigatti
// -- added myMonomials_forC5, mySupport_forC5
//
// Revision 1.53  2010/11/02 15:33:23  bigatti
// -- fixed myIsPrintAtom
//
// Revision 1.52  2010/10/08 11:39:52  abbott
// Renamed DistrMPoly to DistrMPolyClean.
//
// Revision 1.51  2010/10/01 15:52:23  bigatti
// -- added mySymbolValue
//
// Revision 1.50  2010/06/10 08:00:02  bigatti
// -- fixed naming conventions
//
// Revision 1.49  2010/03/30 16:06:10  bigatti
// -- fixed call to characteristic
//
// Revision 1.48  2010/03/30 15:17:51  bigatti
// -- added nasty trick for simplifying multiplication of monomial ideals
//
// Revision 1.47  2010/03/18 13:55:56  abbott
// Added pseudo-ctors for monomials with QQ coeffs.
//
// Revision 1.46  2010/03/05 18:43:48  abbott
// Added pseudo-ctors allowing polynomial rings to be created specifying
// the ordering using a PPOrderingCtor object.
//
// Revision 1.45  2010/02/04 10:14:38  bigatti
// -- fixed "mul" for ideal
//
// Revision 1.44  2010/02/04 09:57:11  bigatti
// -- added "mul" for ideals.  Implemented only for SparsePolyRing
//
// Revision 1.43  2010/02/03 14:15:58  bigatti
// -- moved IsMonomial(vector<RingElem>) into anomymous namespace
//
// Revision 1.42  2010/01/21 13:17:08  bigatti
// -- commented out unnecessary default definition for bool3
//
// Revision 1.41  2010/01/20 16:49:52  bigatti
// -- fixed IhaveMonomialGens3Value in some functions modifying the ideal
//
// Revision 1.40  2009/11/26 17:21:33  bigatti
// -- added PushFront/PushBack(f, c, pp)
// -- added in .C inline functions:
// --   CheckCompatible, CheckElemSparsePolyRing, CheckCoeffExpv, CheckCoeffPP
//
// Revision 1.39  2009/10/29 18:43:23  abbott
// Added necessary include directive for iostream
// (previously was wrongly commented out).
//
// Revision 1.38  2009/10/02 13:27:26  bigatti
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.37  2009/09/28 16:19:43  bigatti
// -- unique implementation for myDeriv
//
// Revision 1.36  2009/09/25 13:09:27  bigatti
// -- fixed: reinserted #include QuotientRing
//
// Revision 1.35  2009/09/25 13:02:09  bigatti
// -- myDiv with one implementation in SparsePolyRing
//
// Revision 1.34  2009/09/24 16:24:32  abbott
// Added include directive (after removing it from RingFp.H).
//
// Revision 1.33  2009/09/22 14:07:33  bigatti
// -- added CmpWDegPartial and IsHomogPartial
//
// Revision 1.32  2009/07/30 15:43:56  bigatti
// -- new implementation of IsMonomial(vector<RingElem>) using STL functions
//
// Revision 1.31  2009/07/24 12:26:42  abbott
// Added CommonDenom function for polynomials.
//
// Revision 1.30  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.29  2009/05/22 10:27:08  bigatti
// -- fixed myOutput and myIsPrintAtom
//
// Revision 1.28  2009/05/21 14:49:33  abbott
// Cleaned up myOutput: now handles -1 in small finite fields correctly
// (if residues are symmetric).
//
// Revision 1.27  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.26  2008/11/18 15:20:09  bigatti
// -- added const to myGBasis return value
// -- added myIdealCtor to RingWeyl for proper inheritance
//
// Revision 1.25  2008/10/09 15:08:42  bigatti
// -- added MonomialIntersection
//
// Revision 1.24  2008/07/04 09:11:04  bigatti
// -- new PPVector class
//
// Revision 1.23  2008/05/29 16:02:32  bigatti
// -- improved myMul for the case [coeff * poly]
//
// Revision 1.22  2008/05/27 16:43:35  bigatti
// -- minor change
//
// Revision 1.21  2008/05/27 16:25:59  bigatti
// -- moved code for monomial ideals into TmpPPVector
// -- small improvement for myMul
//
// Revision 1.20  2008/04/15 15:53:07  bigatti
// -- new for SparsePolyRing: mySquare, first draft
//
// Revision 1.19  2008/04/10 16:36:03  bigatti
// -- prototype for squaring polynomials
//
// Revision 1.18  2008/04/03 13:07:16  bigatti
// -- improved "myPowerSmallExp" for monomial input
//
// Revision 1.17  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.16  2007/12/07 15:27:01  bigatti
// -- default implementation of "IamOne" in ideal.C
//
// Revision 1.15  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.14  2007/12/04 14:27:06  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.13  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.12  2007/10/05 14:37:42  bigatti
// -- just fixed a comment
//
// Revision 1.11  2007/09/24 14:25:14  abbott
// Changed IsIndetPower to IsIndetPosPower.
// Added parens to shut up gcc-4.3.
//
// Revision 1.10  2007/06/25 11:38:15  bigatti
// -- reverted IsIndetPosPower --> IsIndetPower
//
// Revision 1.9  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.8  2007/05/31 16:34:37  bigatti
// -- Changed IsValid (now returns true of false and does not throw an error)
// -- using IsValid for sanity check in PushBack
//
// Revision 1.7  2007/05/31 15:43:56  bigatti
// -- added mySymbols and AreGoodIndetNames
//
// Revision 1.5  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.4  2007/05/21 14:46:34  bigatti
// -- PushFront now accepts zero coefficient
//
// Revision 1.3  2007/03/12 16:00:29  bigatti
// -- moved myLog(F, index) into unique implementation in SparsePolyRing
//
// Revision 1.2  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.32  2007/03/08 18:49:20  bigatti
// -- removed TmpGcd.C: moved code into SparsePolyRing.C
//
// Revision 1.31  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.30  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.29  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.28  2007/03/07 14:08:52  bigatti
// -- minor: commented argument names for -Wextra
//
// Revision 1.27  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.26  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.25  2007/02/28 13:51:59  bigatti
// -- added function IsMonomial
//
// Revision 1.24  2007/02/26 15:00:01  bigatti
// -- just a comment fix
//
// Revision 1.23  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.22  2007/01/20 14:07:25  bigatti
// -- moved code for homomorphism into common implementation in SparsePolyRing
//
// Revision 1.21  2006/12/07 17:36:19  cocoa
// -- migrated  myRemoveBigContent myContent myPowerSmallExp  into
//    single implementation in SparsePolyRing
// -- removed  content  from DistrMPoly(..)
//
// Revision 1.20  2006/11/24 17:06:10  cocoa
// -- reorganized includes of header files
//
// Revision 1.19  2006/11/23 17:48:43  cocoa
// -- minor change
//
// Revision 1.18  2006/11/22 17:51:31  cocoa
// -- moved printing functions into unified implementation in SparsePolyRing
//
// Revision 1.17  2006/11/22 15:12:48  cocoa
// -- minor cleaning (indicated by Intel compiler)
//
// Revision 1.16  2006/11/21 18:09:23  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPoly(..) and RingDistrMPoly(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.15  2006/11/17 18:14:01  cocoa
// -- fixed: now Hilbert computes TidyGens instead of cheating
// -- myGBasis is much more efficient on monomial input
//
// Revision 1.14  2006/11/09 17:46:58  cocoa
// -- version 0.9712:
// --   IdealImpl moved to SparsePolyRing from concrete rings
// -- PolyList in GTypes is now vector<RingElem> (was list)
// -- "my" coding convention applied to DistrMPoly
//
// Revision 1.13  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.12  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.11  2006/09/27 14:30:45  cocoa
// -- changed f into rawf for ConstRawPtr and RawPtr variables
//
// Revision 1.10  2006/08/17 09:45:07  cocoa
// -- added: homogenization
//
// Revision 1.9  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.8  2006/07/20 17:06:08  cocoa
// -- moved myStdDeg into SparsePolyRing
//
// Revision 1.7  2006/07/20 16:51:38  cocoa
// -- added common implementation of myStdDeg
//
// Revision 1.6  2006/07/20 14:24:12  cocoa
// -- added special cases for myMul (0) and myGcd (unit)
//
// Revision 1.5  2006/07/17 19:32:55  cocoa
// Added arg check to PushFront.
//
// Revision 1.4  2006/07/17 11:05:53  cocoa
// -- added: myIsValid, myIsHomogeneous, IsHomogeneous
//
// Revision 1.3  2006/06/20 17:25:27  cocoa
// -- added function geobucket::myAddMul(monom, g, gLen);   [without SkipLMFlag]
//
// Revision 1.2  2006/06/08 16:45:27  cocoa
// -- RingDistrMPoly*.H  have been "moved" into RingDistrMPoly*.C
// -- some coding conventions fixed in DistrMPoly*
// -- functions wdeg and CmpWDeg have a common implementation in SparsePolyRing
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.11  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.10  2006/05/12 17:00:34  cocoa
// -- added some function whose implementation in RingDistr*** were identical
// -- added multivariate gcd
//
// Revision 1.9  2006/05/12 13:16:30  cocoa
// Added functions for preprocessing approximate points.
//
// Revision 1.8  2006/04/27 14:08:24  cocoa
// -- sorry: missed parenthesis
//
// Revision 1.7  2006/04/27 14:05:29  cocoa
// -- changed CmpWDeg: I think it should throw an error if polys are 0
//
// Revision 1.6  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPoly* and DistrMPoly* have been disabled
//
// Revision 1.5  2006/04/21 14:58:04  cocoa
// Removed myWDeg member function: it is no longer needed.
//
// Revision 1.4  2006/04/11 14:16:29  cocoa
// -- reorganization of fns between reduce,SparsePolyRing,GPoly
// -- added optional "len" argument to myAssignReset in ReductionCog
//
// Revision 1.3  2006/02/13 13:19:01  cocoa
// -- fixed: "const PPMonoidElem&" --> "ConstRefPPMonoidElem"
//
// Revision 1.2  2006/01/19 16:34:42  cocoa
// -- added NF, myReduceMod functions (not yet tested)
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.2  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.1  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
