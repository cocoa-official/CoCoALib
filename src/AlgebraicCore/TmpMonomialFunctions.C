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
