//   Copyright (c)  2008-2009  Anna Bigatti and Eduardo Saenz-de-Cabezon

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

#include "CoCoA/TmpMonomialFunctions.H"

#include "CoCoA/PPMonoid.H"  // for MultidegreeMap
#include "CoCoA/TmpPPVector.H"
#include "CoCoA/BigInt.H"  // for position_t

using std::vector;

/**********************************************************/

namespace CoCoA
{

  void support(std::vector<long>& sup, const PPMonoidElem& monomial)
  {
    const std::vector<long> expv = exponents(monomial); 
    sup.clear();
    const long size = len(expv);
    for (long i=0; i<size; ++i)
      if (expv[i] != 0)
        sup.push_back(i);
  }

  bool IsIrreducible(const PPVector& ideal)
  {
    long numgens=len(ideal);
    for (long cont=0; cont<numgens; ++cont)
      if (!IsIndetPosPower(PP(ideal[cont])))
        return false;
    return true;
  }

  bool IsPrime(const PPVector& ideal)
  {
    long numgens=len(ideal);
    long var;
    for (long cont=0; cont<numgens; ++cont)
      if (!IsIndet(var, PP(ideal[cont])))
        return false;
    return true;
  }

  bool IsPrimary(const PPVector& ideal)
  {
    long numgens=len(ideal);
    long var;
    BigInt expo;
    long N=NumIndets(PPM(ideal));
    vector<bool> pures(N,false);
    vector<bool> other(N,false);
    for (long cont=0; cont<numgens; ++cont)
    {
      if (IsIndetPosPower(var, expo, PP(ideal[cont])))
        pures[var] = true;
      else //I copied part of "support"
      {
        const std::vector<long> expv = exponents(PP(ideal[cont]));
        for (long i=0; i<N; ++i)
          if (expv[i] != 0)
            other[i]=true;
      }
    }
    for (long conta=0; conta<N; ++conta)
      if (pures[conta]==false && other[conta]==true)
        return false;
    return true;
  }


  // done! 2011-07-04
//   void ColonIdeal(PPVector& PPs, const PPVector& PPs1, const PPVector& PPs2)
//   {
//     //is it relevant which one has more generators?
//     PushBack(PPs,one(PPM(PPs)));
//     for(long i=0; i<len(PPs2); ++i)
//     {
//       PPVector current(PPM(PPs),DMR(PPs));
//       for (long j=0; j<len(PPs1); ++j)
//         PushBack(current, colon(PP(PPs1[j]),PP(PPs2[i])));
//       interreduce(current); //is it better to just add and then interreduce at the end?
//       lcms(PPs, PPs, current);
//     }
//     interreduce(PPs);
//   }
  


}  // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpMonomialFunctions.C,v 1.16 2023/03/09 22:41:07 abbott Exp $
// $Log: TmpMonomialFunctions.C,v $
// Revision 1.16  2023/03/09 22:41:07  abbott
// Summary: Tidied use of exponents (redmine 1694)
//
// Revision 1.15  2022/10/11 15:39:35  bigatti
// Summary: using new interface of exponents returning vector
//
// Revision 1.14  2022/02/18 14:12:01  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.13  2021/07/19 11:26:08  abbott
// Summary: Layout
//
// Revision 1.12  2014/04/30 16:24:27  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.11  2012/02/01 13:31:39  abbott
// Removed several unnecessary includes.
//
// Revision 1.10  2011/12/23 14:57:11  bigatti
// -- changed and --> && , or --> ||
//
// Revision 1.9  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.8  2011/07/05 15:02:17  bigatti
// -- added AlexanderDual
// -- added ad-hoc functions for colon, elim on monomial ideals
//
// Revision 1.7  2011/03/22 15:30:07  bigatti
// -- added cvs log
//
