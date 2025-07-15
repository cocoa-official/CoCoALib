//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti
//   Original author: 2015  Mario Albert

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

#include "CoCoA/TmpPBMill.H"
#include "CoCoA/QuotientRing.H"


namespace CoCoA
{
  namespace Involutive
  {
    using std::list;
    using std::map;
    using std::multimap;
    using std::make_pair;
    using std::pair;
    using std::unique_ptr;
    using std::vector;
//    using namespace std;

    PBMill::PBMill(const PBMill::Builder& builder) : JBMill(builder.build())
    {
      if (!IamPommaretBasis())
      {
        CoCoA_THROW_ERROR1("This is not a Pommaret Basis");
      }
    }

    long PBMill::myCountElementsWithDegree(const vector<RingElem>& vec, long degree) const
    {
      long counter(0);
      for (const RingElem& g: vec)
      {
        if (StdDeg(g) == degree)
        {
          ++counter;
        }
      }
      return counter;
    }

    long PBMill::myDepth() const
    {
      long res(0);
      if (IamHomogenous())
      {
        map<PPMonoidElem, vector<bool> > MultVars = myMultVars();
        long depth(NumIndets(myPolyRing));
        for (const auto& entry: MultVars)
        {
          long CurrentMultVars(myCountTrues(entry.second));
          if (depth > CurrentMultVars)
          {
            depth = CurrentMultVars;
          }
        }
        res = depth;
        // --res;
      }
      else
      {
        CoCoA_THROW_ERROR1("Janet basis isn't a Pommaret basis or the ideal isn't homogenous");
      }
      return res;
    }

    long PBMill::myProjDim() const
    {
      return NumIndets(myPolyRing) - myDepth(); //depth throws the error if not pommaret or homogenous
    }

    vector<RingElem> PBMill::mySocle() const
    {
      if (!(IamHomogenous()))
      {
        CoCoA_THROW_ERROR1(ERR::ReqHomog);
      }

      if (!(IamCohenMacaulay()))
      {
        CoCoA_THROW_ERROR1("Ideal isn't CohenMacaulay");
      }
      if (!(IsStdDegRevLex(ordering(myPPMValue))))
      {
        CoCoA_THROW_ERROR2(ERR::PPOrder, "Require degrevlex ordering");
      }

      long min = myMinCls();
      vector<RingElem> socle = myElementsWithClass(min);
      vector<RingElem> ResultWithoutResidueClass;
      for (const RingElem& g: socle)
      {
        // RingElem rest = myStandardRepresentation(*i).second;
        bool NoConstant(true);
        for (int iter = 0; iter < (min - 1); ++iter)
        {
          RingElem rest = myStandardRepresentation(g * indet(myPolyRing, iter)).second;
          while(!IsZero(rest))
          {
            if (IsConstant(rest))
            {
              NoConstant = false;
              break;
            }
            rest = rest - monomial(myPolyRing,LC(rest),LPP(rest));
          }
          if (!NoConstant)
          {
            break;
          }
        }
        if (NoConstant)
        {
          ResultWithoutResidueClass.push_back(g / indet(myPolyRing, min - 1));
        }
      }
      long dimension = myDim();
      vector<RingElem> x = indets(myPolyRing);
      vector<RingElem> GenSet;
      for (long i = 0; i < dimension; ++i)
      {
        GenSet.push_back(x[NumIndets(myPolyRing) - 1 - i]);
      }
      ring quotientRing = NewQuotientRing(myPolyRing, ideal(myPolyRing, GenSet) );
      RingHom phi = QuotientingHom(quotientRing);
      vector<RingElem> res;
      for (const RingElem& g: ResultWithoutResidueClass)
      {
        res.push_back(phi(g));
      }
      return res;
    }

    map<pair<long, long>, long> PBMill::myExtremalBettiNumbers() const
    {
      if (!(IamHomogenous()))
      {
        CoCoA_THROW_ERROR1("Janet basis isn't homogenous");
      }
      if (!(IsStdDegRevLex(ordering(myPPMValue))))
      {
        CoCoA_THROW_ERROR2(ERR::PPOrder, "Require degrevlex ordering");
      }

      long CurPDim = myProjDim();
      long n(NumIndets(myPolyRing));
      map<long, long> DegOfPommMultVars(myDegSeperatedByNumMultVar());
      long deg(-1);
      long numMultVars(-1);
      for(long i = 1; i != n + 1; ++i)
      {
        if (deg < DegOfPommMultVars[i])
        {
          deg = DegOfPommMultVars[i];
          numMultVars = i;
        }
      }
      map<pair<long, long>, long> result;
      vector<RingElem>  ElementsWithNumMultVars(myElementsWithClass(n - numMultVars));
      result.insert(pair<pair<long, long>, long>(pair<long,long>(n - numMultVars, deg + n - numMultVars), myCountElementsWithDegree(ElementsWithNumMultVars, deg)));
      while( n - numMultVars < CurPDim)
      {
        deg = 0;
        long OldCls(numMultVars);
        for(long i = 1; i != OldCls; ++i)
        {
          if (deg < DegOfPommMultVars[i])
          {
            deg = DegOfPommMultVars[i];
            numMultVars = i;
          }
        }
        ElementsWithNumMultVars = myElementsWithClass(numMultVars);
        result.insert(pair<pair<long, long>, long>(pair<long,long>(n - numMultVars, deg + n - numMultVars), myCountElementsWithDegree(ElementsWithNumMultVars, deg)));
      }
      return result;
    }

    vector<RingElem> PBMill::myRegSeq() const//mod I!!!
    {
      vector<RingElem> x = indets(myPolyRing);
      vector<RingElem> res;
      for(long iter = 0; iter !=  myDepth(); ++iter)//throws an error if not a pommaret basis or homogenous
      {
        res.push_back(x[iter]);
      }
      return res;
    }

    vector<RingElem> PBMill::myMaxStronglyIndependentSet() const
    {
      vector<RingElem> res;
      if (IsOne(LPP(myBasis.myListBegin()->myGetPol())))
      {
        return res;
      }
      ring quotientRing = NewQuotientRing(myPolyRing, ideal(myPolyRing, myReturnGB()) );
      RingHom phi = QuotientingHom(quotientRing);
      long CountIndets = NumIndets(myPolyRing);
      vector<RingElem> x = indets(myPolyRing);
      for(long iter = CountIndets - myDim(); iter != CountIndets; ++iter)
      {
        res.push_back(phi(x[iter]));
      }
      return res;
    }

    bool PBMill::IamCohenMacaulay() const
    {
      if (myDepth() == myDim()) //throws an error if it is not a pommaret basis or not homogenous
      {
        return true;
      }

      return false;
    }

    long PBMill::myRegularity() const
    {
      if (!(IsStdDegRevLex(ordering(myPPMValue))))
      {
        CoCoA_THROW_ERROR2(ERR::PPOrder, "Require degrevlex ordering");
      }
      if (!IamHomogenous())
      {
        CoCoA_THROW_ERROR1(ERR::ReqHomog);
      }

      long MaxDeg(0);
      for(list<JanetTriple>::const_iterator iter = myBasis.myListBegin(); iter != myBasis.myListEnd(); ++iter)
      {
        const long d = StdDeg(iter->myGetPol());
        if (d > MaxDeg) MaxDeg = d;
      }
      return MaxDeg;
    }

    vector<RingElem> PBMill::mySaturation() const
    {
      if (!IamHomogenous())
      {
        CoCoA_THROW_ERROR1("Janet basis isn't isn't homogenous");
      }
      if (!(IsStdDegRevLex(ordering(myPPMValue))))
      {
        CoCoA_THROW_ERROR2(ERR::PPOrder, "Require degrevlex ordering");
      }
      // initialization
      multimap<long, RingElem> SepBasis = myBasisSeperatedByNumberOfMultVar();
      long n = NumIndets(myPolyRing);
      long LastVarIndex = n - 1;
      vector<RingElem> res;
      //handling H_1
      for (multimap<long, RingElem>::iterator iter = SepBasis.lower_bound(1); iter != SepBasis.upper_bound(1); ++iter)
      {
        // remove x_n part
        res.push_back(iter->second / IndetPower(myPolyRing, LastVarIndex, exponent(LPP(iter->second), LastVarIndex)));
      }
      //handling rest of H_n
      for (int i = 2; i <= n; ++i)
      {
        for (multimap<long, RingElem>::iterator iter = SepBasis.lower_bound(i); iter != SepBasis.upper_bound(i); ++iter)
        {
          res.push_back(iter->second);
        }
      }
      // not minimized!
      return res;
    }

    long PBMill::mySatiety() const
    {
      if (!IamHomogenous())
      {
        CoCoA_THROW_ERROR1("Janet basis isn't isn't homogenous");
      }
      if (!(IsStdDegRevLex(ordering(myPPMValue))))
      {
        CoCoA_THROW_ERROR2(ERR::PPOrder, "Require degrevlex ordering");
      }
      // maximal degree of elements with only one mult var
      return myDegForNumMultVar(1);
    }

    multimap<long, RingElem> PBMill::myBasisSeperatedByNumberOfMultVar() const
    {
      vector<RingElem> pb(myReturnPB());
      multimap<long, RingElem> res;
      long n = NumIndets(myPPMValue);
      for (const RingElem& g: pb)
      {
        // the coordinates are delta regular, hence the number of mult vars
        // of the given element r is NumIndets - cls(r)
        res.insert(make_pair(n - myCls(LPP(g)), g));
      }

      return res;
    }

    long PBMill::myDegForNumMultVar(long index) const
    {
      multimap<long, RingElem> MultVarSeperation(myBasisSeperatedByNumberOfMultVar());
      long res(0);
      for(multimap<long, RingElem>::iterator iter = MultVarSeperation.lower_bound(index); iter != MultVarSeperation.upper_bound(index); ++iter)
      {
        if (res < StdDeg(iter->second))
        {
          res = StdDeg(iter->second);
        }
      }
      return res;
    }

    map<long, long> PBMill::myDegSeperatedByNumMultVar() const
    {
      long CountIndets(NumIndets(myPolyRing));
      map<long, long> res;
      for(long iter = 1; iter <= CountIndets; ++iter)
      {
        res.insert(pair<long, long>(iter, myDegForNumMultVar(iter)));
      }
      return res;
    }

    unique_ptr<StabilityAlgorithm> PBMill::DeltaRegularTransformator::myChooseStrategy() const
    {
      unique_ptr<StabilityAlgorithm> res;
      switch(myStrategy)
      {
      case  SingleWithPermutation:
        res = unique_ptr<StabilityAlgorithm>(new DeltaRegular(myPolyRing(), DeltaRegular::UsePermutations));
        break;
      case  SingleWithoutPermutation:
        res = unique_ptr<StabilityAlgorithm>(new DeltaRegular(myPolyRing(), DeltaRegular::NotUsePermutations));
        break;
      case  AllWithPermutation:
        res = unique_ptr<StabilityAlgorithm>(new DeltaRegularAll(myPolyRing(), DeltaRegular::UsePermutations));
        break;
      case  AllWithoutPermutation:
        res = unique_ptr<StabilityAlgorithm>(new DeltaRegularAll(myPolyRing(), DeltaRegular::NotUsePermutations));
        break;
      }
      return res;
    }


    unique_ptr<StabilityAlgorithm> PBMill::StableLTITransformator::myChooseStrategy() const
    {

      unique_ptr<StabilityAlgorithm> res;
      switch(myStrategy)
      {
      case  Single:
        res = unique_ptr<StabilityAlgorithm>(new StableLTI(myPolyRing()));
        break;
      case  All:
        res = unique_ptr<StabilityAlgorithm>(new StableLTIAll(myPolyRing()));
        break;
      }
      return res;
    }


    unique_ptr<StabilityAlgorithm> PBMill::StronglyStableLTITransformator::myChooseStrategy() const
    {

      unique_ptr<StabilityAlgorithm> res;
      switch(myStrategy)
      {
      case  Single:
        res = unique_ptr<StabilityAlgorithm>(new StronglyStableLTI(myPolyRing()));
        break;
      case  All:
        res = unique_ptr<StabilityAlgorithm>(new StronglyStableLTIAll(myPolyRing()));
        break;
      }
      return res;
    }


    JanetContainer PBMill::Transformator::myUseStrategy() const
    {
      unique_ptr<StabilityAlgorithm> strategy(myChooseStrategy());
      strategy->SetStatisticLevel(StatLevel);
      strategy->SetJanetStrategy(myJanetStrategy);
      strategy->myComputer(myInput.begin(), myInput.end());
      return strategy->myOutputResult();
    }

    pair<JanetContainer, StabilityAlgorithm::Statistics> PBMill::Transformator::myUseStrategyWithStatistics()
    {
      unique_ptr<StabilityAlgorithm> strategy(myChooseStrategy());
      strategy->SetStatisticLevel(StatLevel);
      strategy->SetJanetStrategy(myJanetStrategy);
      strategy->myComputer(myInput.begin(), myInput.end());
      return make_pair(strategy->myOutputResult(), strategy->myOutputStatistics());
    }


  } // end of namespace Involutive
} // end of namespace CoCoA
