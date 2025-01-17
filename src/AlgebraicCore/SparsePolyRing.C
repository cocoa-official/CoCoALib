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
