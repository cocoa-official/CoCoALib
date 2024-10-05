//   Copyright (c)  2018  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/SparsePolyOps-EvalInterval.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/BigRatInterval.H"
#include "CoCoA/utils.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/verbose.H"

#include <algorithm>
using std::max;
using std::min;
#include <iostream>
#include <utility>
using std::pair;
using std::make_pair;
//#include <vector>
using std::vector;

namespace CoCoA
{

  std::pair<BigRatInterval, BigRatInterval> SplitIntvl(const BigRatInterval& I)
  {
    const BigRat MidPoint = (min(I) + max(I))/2;
    return make_pair(BigRatInterval(min(I), MidPoint),
                     BigRatInterval(MidPoint, max(I)));
  }

  BigRatInterval EvalByHorner(const std::vector<BigRat>& CoeffVec, const BigRatInterval& I)
  {
    if (IsZeroInside(I))
    {
      BigRatInterval L(min(I), BigRat(0));
      BigRatInterval R(BigRat(0), max(I));
      return merge(EvalByHorner(CoeffVec, L), EvalByHorner(CoeffVec, R));
    }
    // Here the interval I does not contain 0 in its interior
    const int deg = len(CoeffVec)-1;
    BigRatInterval ans(CoeffVec[deg], CoeffVec[deg]);
    for (int k=deg-1; k >= 0; --k)
    {
      ans = ans*I + CoeffVec[k];
    }
    return ans;
  }


  BigRatInterval eval(const vector<BigRat>& f, const BigRatInterval& I)
  {
    VerboseLog VERBOSE("eval(poly,intvl)");
    vector<BigRatInterval> ListXIntvl;
    if (!IsZeroInside(I))
      ListXIntvl.push_back(I);
    else
    {
      VERBOSE(20) << "Intvl straddles 0; splitting" << std::endl;
      BigRatInterval L(min(I), BigRat(0));
      BigRatInterval R(BigRat(0), max(I));
      ListXIntvl.push_back(L);
      ListXIntvl.push_back(R);
    }

    vector<BigRatInterval> ListValIntvl;
    int n = len(ListXIntvl);
    for (int k=0; k < n; ++k)
      ListValIntvl.push_back(EvalByHorner(f, ListXIntvl[k]));

    BigRat UPB = max(ListValIntvl[0]);
    BigRat LWB = min(ListValIntvl[0]);
    for (int k=1; k < n; ++k)
    {
      UPB = max(UPB, max(ListValIntvl[k]));
      LWB = min(LWB, min(ListValIntvl[k]));
    }

    // main loop
    while (true)
    {
      const BigRat ValWidth = UPB - LWB;

      VERBOSE(20) << "MAIN LOOP:" << std::endl;
      VERBOSE(20) << "X intvls: " << ListXIntvl << std::endl;
      VERBOSE(20) << "Val intvls: " << ListValIntvl << std::endl;
      VERBOSE(20) << "LWB = " << LWB << std::endl;
      VERBOSE(20) << "UPB = " << UPB << std::endl;
      VERBOSE(20) << "ValWidth = " << ValWidth << std::endl;
      
      vector<long> ToRefineUPB;
      vector<long> ToRefineLWB;
      bool RecomputeLWB = false;

      for (int k=0; k < n; ++k)
      {
        const bool ReachesUPB = (max(ListValIntvl[k]) == UPB);
        const bool ReachesLWB = (min(ListValIntvl[k]) == LWB);
        if (ReachesUPB && ReachesLWB) RecomputeLWB = true;
        if (ReachesUPB)
          ToRefineUPB.push_back(k);
        if (ReachesLWB)
          ToRefineLWB.push_back(k);
      }          

      VERBOSE(20) << "ToRefineUPB = " << ToRefineUPB << std::endl;
      const int NumRefineUPB = len(ToRefineUPB);
      BigRat LocalUPB = LWB;
      for (int j=0; j < NumRefineUPB; ++j)
      {
        const int i=ToRefineUPB[j];
        pair<BigRatInterval,BigRatInterval> LR = SplitIntvl(ListXIntvl[i]);
        const BigRatInterval ValL = EvalByHorner(f, LR.first);
        const BigRatInterval ValR = EvalByHorner(f, LR.second);
        LocalUPB = max(LocalUPB, max(max(ValL), max(ValR)));
        ListXIntvl[i] = LR.first;
        ListValIntvl[i] = ValL;
        ListXIntvl.push_back(LR.second);
        ListValIntvl.push_back(ValR);
      }

      VERBOSE(20) << std::endl;
      VERBOSE(20) << "AFTER REFINE UPB:" << std::endl;
      VERBOSE(20) << "X intvls: " << ListXIntvl << std::endl;
      VERBOSE(20) << "Val intvls: " << ListValIntvl << std::endl;

      if (RecomputeLWB)
      {
        n = len(ListXIntvl);
        BigRat LWB = min(ListValIntvl[0]);
        for (int k=1; k < n; ++k)
        {
          LWB = min(LWB, min(ListValIntvl[k]));
        }
        VERBOSE(20) << "Recomputing LWB..." << LWB << std::endl;

        ToRefineLWB.clear();
        for (int k=0; k < n; ++k)
        {
//        const bool ReachesUPB = (max(ListValIntvl[k]) == UPB);
          const bool ReachesLWB = (min(ListValIntvl[k]) == LWB);
//        if (ReachesUPB && ReachesLWB) RecomputeLWB = true;
//        if (ReachesUPB)
//          ToRefineUPB.push_back(k);
          if (ReachesLWB)
            ToRefineLWB.push_back(k);
        }          

      }

      VERBOSE(20) << "ToRefineLWB = " << ToRefineLWB << std::endl;
      const int NumRefineLWB = len(ToRefineLWB);
      BigRat LocalLWB = UPB;
      for (int j=0; j < NumRefineLWB; ++j)
      {
        const int i=ToRefineLWB[j];
        pair<BigRatInterval,BigRatInterval> LR = SplitIntvl(ListXIntvl[i]);
        const BigRatInterval ValL = EvalByHorner(f, LR.first);
        const BigRatInterval ValR = EvalByHorner(f, LR.second);
        LocalLWB = min(LocalLWB, min(min(ValL), min(ValR)));
        ListXIntvl[i] = LR.first;
        ListValIntvl[i] = ValL;
        ListXIntvl.push_back(LR.second);
        ListValIntvl.push_back(ValR);
      }

      if (LocalUPB-LocalLWB > BigRat(15,16)*ValWidth) break;

      n = len(ListXIntvl);
    UPB = max(ListValIntvl[0]);
    LWB = min(ListValIntvl[0]);
    for (int k=1; k < n; ++k)
    {
      UPB = max(UPB, max(ListValIntvl[k]));
      LWB = min(LWB, min(ListValIntvl[k]));
    }
    VERBOSE(20) << "New LWB = " << LWB << std::endl;
    VERBOSE(20) << "New UPB = " << UPB << std::endl;
    VERBOSE(20) << "END ITER " << std::endl;

//    if (UPB-LWB > 65535*ValWidth/65536) break;

  }
    VERBOSE(10) << "END  len = " << len(ListXIntvl) << std::endl;
    return BigRatInterval(LWB,UPB);
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-EvalInterval.C,v 1.6 2022/02/18 14:11:58 abbott Exp $
// $Log: SparsePolyOps-EvalInterval.C,v $
// Revision 1.6  2022/02/18 14:11:58  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2021/01/15 16:59:34  abbott
// Summary: Redmine 1563: new ctor for BigRat direct from integer
//
// Revision 1.4  2018/05/17 15:45:25  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.3  2018/05/05 15:24:43  abbott
// Summary: Decreased "reduction factor" in main loop; faster, but looser result
//
// Revision 1.2  2018/04/23 09:00:13  abbott
// Summary: Minor improvements
//
// Revision 1.1  2018/04/20 13:12:15  abbott
// Summary: Added new fn eval(poly, interval) -- still only a prototype
//
//
