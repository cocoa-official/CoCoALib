//   Copyright (c)  2024  John Abbott and Anna M. Bigatti
//   Main authors: Anna M Bigatti, Evelina Lanteri

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

#include "CoCoA/SparsePolyOps-ideal-FGLM.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/LinDepMill.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/QBGenerator.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"
#include "CoCoA/verbose.H" // for VerboseLog
#include <algorithm> // for swap

// #include <algorithm> // already included in apply.H
// using std::swap;   in  SOI, NBM
#include <iostream>
using std::ostream;
using std::endl;
//#include <vector>  --- already included by ApproxPts2.H
using std::vector;
#include <list>
using std::list;


namespace CoCoA
{

  namespace // anonymous namespace for file local auxiliary functions and definitions.
  {
    
    void AddToQBAndUpdate(QBGenerator& QBG,
                          vector<RingElem>& QB_P,
                          const PPMonoidElem& t)
    {
      VerboseLog VERBOSE("AddToQBAndUpdate");
      SparsePolyRing P = owner(QB_P.front());
      QB_P.push_back(monomial(P,t));
      QBG.myCornerPPIntoQB(t);
      VERBOSE(90) << t << " added to QBG and QB_P" << endl;
    }
    
    
    void AddToGBAndUpdate(QBGenerator& QBG,
                          vector<RingElem>& GB,
                          const PPMonoidElem& t,
                          const vector<RingElem>& v, // in K
                          const vector<RingElem>& QB_P)
    {
      VerboseLog VERBOSE("AddToCornerAndUpdate");
      SparsePolyRing P = owner(QB_P.front());
      RingElem f = monomial(P, t);
      for (long i=0; i < len(QB_P); ++i)      f += v[i]*QB_P[i];
      GB.push_back(f);
      QBG.myCornerPPIntoAvoidSet(t);
      VERBOSE(90) << t << " added to AvoidSet, and poly to GB" << endl;
    }

    //// copiata da SparsePolyOps-MinPoly.C perche' in anonymous namespace

    // >>>ASSUMES<<< that QB is in increasing term-order!
    void coefficients(vector<RingElem>& coeffs, ConstRefRingElem f, const vector<PPMonoidElem>& QB)
    {
      if (len(coeffs)!=len(QB))  CoCoA_THROW_ERROR2(ERR::IncompatDims, "len(coeffs)!=len(QB)");
      SparsePolyIter itf=BeginIter(f);
      long i=len(QB)-1;
      for (; i>=0 && !IsEnded(itf); --i)
      {
        CoCoA_ASSERT(PP(itf) <= QB[i]);
        if (PP(itf) == QB[i])
        {
          coeffs[i] = coeff(itf);
          ++itf;
        }
        else
          coeffs[i] = 0;
      }
      for (; i>=0; --i)  coeffs[i] = 0;
    }


  }  // anonymous namespace -- end
  
//--------------  FGLM  ---------------------------------------------

  std::vector<RingElem> FGLM(const SparsePolyRing& P_new, const ideal& I)
  {
    std::vector<PPMonoidElem> QB;
    std::vector<RingElem> GB;
    FGLM(QB, GB, P_new, I);
    return GB;
  }

  
  void FGLM(std::vector<PPMonoidElem>& QB,
            std::vector<RingElem>& GB,
            const SparsePolyRing& P_new,
            const ideal& I)   // in P_old
  {
    VerboseLog VERBOSE("FGLM");
    //    CheckInput(P_new, I);  // da fare: compatibilita', HasGBasis(I), I=(1), ..

    const SparsePolyRing& P_old = RingOf(I);
    ring K = CoeffRing(P_new);
      
    //   --------------     INITIALIZE   ---------------
    long IterNum = 1;
    vector<RingElem> GB_new, QB_new;

    QBGenerator QBG(PPM(P_new));         // to generate the new QB
    vector<PPMonoidElem> QB_old = QuotientBasisSorted(I);
    LinDepMill LDM(K, len(QB_old)); // to search for lin dep
    vector<RingElem> coeffs(len(QB_old), zero(K));

    // add pp 1
    QBG.myCornerPPIntoQB(one(PPM(P_new)));
    QB_new.push_back(one(P_new));
    coeffs[0] = one(K);    LDM.myAppendVec(coeffs);

    // ------------   MAIN CYCLE   ----------------- verbosity 80/90/95
    while (!QBG.myCorners().empty())
    {
      VERBOSE(95) << "---------------------------------------" << endl;
      VERBOSE(95) << "New QB in progress = " << QBG.myQB() << endl;
      VERBOSE(95) << "List of PP yet to test = " << QBG.myCorners() << endl;
      const PPMonoidElem t = QBG.myCorners().front(); // PP considered in the current step
      VERBOSE(80) << "Iteration " << IterNum << " PP = " << t << endl;
      
      coefficients(coeffs, NF(monomial(P_old, exponents(t)), I), QB_old);
      LDM.myAppendVec(coeffs);
      if (LDM.myLinReln().empty()) // independent
        AddToQBAndUpdate(QBG, QB_new, t);
      else
        AddToGBAndUpdate(QBG, GB_new, t, LDM.myLinReln(), QB_new);

      ++IterNum;
    }
    VERBOSE(95) << "==========================================" << endl
                << "All power-products have been considered" << endl;
    QB = QBG.myQB();  
    std::swap(GB, GB_new);
  }


} // end of namespace CoCoA
