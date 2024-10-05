//   Copyright (c)  2014-2015  John Abbott,  Anna M. Bigatti
//   Original authors: 2014-2015  Pierre Pytlik, Mario Albert

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

#include "CoCoA/TmpStabilityAlgorithm.H"

#include "CoCoA/SparsePolyOps-ideal.H"


namespace CoCoA
{
  namespace Involutive
  {
    using std::map;
    using std::set;
    using std::vector;

    
    void StabilityAlgorithm::myApplyHomToVector(std::vector<RingElem>::iterator beg, std::vector<RingElem>::iterator end, RingHom& hom) const
    {
      for (std::vector<RingElem>::iterator i = beg; i != end; ++i)
      {
        *i = hom(*i);
      }
    }

    void DeltaRegular::myComputer(vector<RingElem>::const_iterator beginInput,
                                  vector<RingElem>::const_iterator endInput)
    {
      JBMill::Builder builder;
      builder.setStrategy(myJanetStrategy);
      vector<vector<RingElem> > result;
      vector<RingElem> generators(beginInput, endInput);
      while(true)
      {
        vector<RingElem> images(indets(myPolyRing));
        if (myUsagePermutations == UsePermutations)
        {
          images = ReturnPermutation(generators);
          myTrackTransformation(images);
          myPerformTransformation(images, generators.begin(), generators.end());
        }
        JBMill mill(builder.setInput(generators));

        myTrackJBMill(mill);
        if (mill.IamPommaretBasis())
        {
          myContainer.myChangeBasis(mill.myGetJanetContainer().myListBegin(),
                                    mill.myGetJanetContainer().myListEnd());
          break;
        }
        images = DoComputeImage(mill);

        generators = mill.myReturnGB();

        myTrackTransformation(images);
        myPerformTransformation(images, generators.begin(), generators.end());
      }
    }

    void DeltaRegular::myPerformTransformation(const vector<RingElem>& images,
                                               vector<RingElem>::iterator beginGens,
                                               vector<RingElem>::iterator endGens) const
    {
      if (images != indets(myPolyRing))
      {
        RingHom phi = PolyAlgebraHom(myPolyRing, myPolyRing, images);
        myApplyHomToVector(beginGens, endGens, phi);
      }
    }


    vector<RingElem> DeltaRegular::ReturnPermutation(const vector<RingElem>& generators) const
    {
      vector<RingElem> gens(generators);
      vector<RingElem>::iterator gensBegin(gens.begin());
      vector<RingElem>::iterator gensEnd(gens.end());

      gensEnd = remove_if(gensBegin, gensEnd, myIsNotIndetPower);
      sort(gensBegin, gensEnd, mySortByPowerAndIndex);

      vector<RingElem> result;
      for (vector<RingElem>::iterator i = gensBegin; i != gensEnd; ++i)
      {
        result.push_back(indet(owner(*i), UnivariateIndetIndex(*i)));
      }
      vector<RingElem> x = indets(myPolyRing);
      const long n = len(x);
      result.insert(result.end(), x.begin(), x.end());
      if(len(result) != n)
      {
        for (std::vector<RingElem>::iterator i = result.begin(); i != result.end(); ++i)
        {
          std::vector<RingElem>::iterator j = i;
          ++j;
          while(j != result.end())
          {
            if (*i == *j)
            {
              CoCoA_ASSERT(j != i);
              j = result.erase(j);
            } else {
              ++j;
            }
          }
        }
      }

      return result;
    }

    vector<RingElem> DeltaRegular::DoComputeImage(const JBMill& mill) const
    {
      vector<RingElem> x = indets(myPolyRing);
      const long nvars = NumIndets(myPolyRing);
      vector<RingElem> images = x;
      vector<long> MinJMultVarNotPMult(nvars, nvars);

      for (const auto& gen: mill.myMultVars())
      {
        const long cls(mill.myCls(gen.first));
        for (long var = 0; var < cls; ++var)
        {
          if ((gen.second)[var] && MinJMultVarNotPMult[cls] > var)
          {
            MinJMultVarNotPMult[cls] = var;
          }
        }
      }

      for (int var = 0; var < nvars; ++var)
      {
        if (MinJMultVarNotPMult[var] < nvars)
        {
          images[var] += x[MinJMultVarNotPMult[var]];
        }
      }
      return images;
    }

    bool DeltaRegular::mySortByPowerAndIndex(ConstRefRingElem r1, ConstRefRingElem r2)
    {
      long i1(UnivariateIndetIndex(r1));
      long i2(UnivariateIndetIndex(r2));
      long n1(deg(r1));
      long n2(deg(r2));

      CoCoA_ASSERT(i1 != -1 && i2 != -1);

      return (n1 < n2) || ((n1 == n2) && i1 < i2);
    }

    vector<RingElem> DeltaRegularAll::DoComputeImage(const JBMill& mill) const
    {
      vector<RingElem> x = indets(myPolyRing);
      const long nvars = NumIndets(myPolyRing);
      vector<RingElem> images = x;

      vector<set<long> > JMultVarsNotPMult(nvars);

      for (const auto& gen: mill.myMultVars())
      {
        const long cls(mill.myCls(gen.first));
        for (long var = 0; var < cls; ++var)
        {
          if ((gen.second)[var])
          {
            JMultVarsNotPMult[cls].insert(var);
          }
        }
      }

      for (int var = 0; var < nvars; ++var)
      {
        for (long i: JMultVarsNotPMult[var])
        {
          images[var] += x[i];
        }
      }

      return images;
    }

    bool StableLTI::myIsLTIStable(const vector<RingElem>& generators) const
    {
      const vector<RingElem>& x = indets(myPolyRing);
      const ideal LTIdeal = LT(ideal(generators));
      const vector<RingElem> GensOfLTIdeal = myLeadingMonomials(GBasis(LTIdeal));
      CoCoA_ASSERT(AreMonomials(GensOfLTIdeal));
      for (const RingElem& m: GensOfLTIdeal)
      {
        const long cls(myCls(m));
        const RingElem g(m / x[cls]);
        for (long i = 0; i < cls; i++)
        {
          if (!IsElem(g * x[i], LTIdeal))
            return false;
        }
      }
      return true;
    }


    bool StronglyStableLTI::myIsLTIStronglyStable(const vector<RingElem>& generators) const
    {
      const vector<RingElem>& x = indets(myPolyRing);
      const ideal LTIdeal = LT(ideal(generators));
      const vector<RingElem> GensLTIdeal = myLeadingMonomials(GBasis(LTIdeal));
      for (const RingElem& m: GensLTIdeal)
      {
        CoCoA_ASSERT(IsMonomial(m));
        const long cls(myCls(m));
        RingElem g(myPolyRing);
        for (long i = 1; i <= cls; i++)
        {
          if (!IsDivisible(g, m, x[i])) continue;
          for (long j = 0; j < i; j++)
          {
            if (!IsElem(g * x[j], LTIdeal))
              return false;
          }
        }
      }
      return true;
    } // bool myIsLTIStronglyStable

    long StableLTI::myCls(ConstRefRingElem r) const
    {
      vector<long> exps;
      exponents(exps, LPP(r));
      long cls = len(exps) - 1;
      for (vector<long>::const_reverse_iterator i = exps.rbegin(); i != exps.rend(); ++i)
      {
        if (*i != 0)
        {
          return cls;
        } else {
          --cls;
        }
      }
      return cls;
    }

    void StableLTI::myComputer(vector<RingElem>::const_iterator beginInput,
                               vector<RingElem>::const_iterator endInput)
    {
      vector<RingElem> NewGenerators(beginInput, endInput);
      while (!IamInPosition(NewGenerators))
      {
        ideal LTIdeal = ideal(myLeadingMonomials(GBasis(ideal(NewGenerators))));
        vector<RingElem> GensLTIdeal = myLeadingMonomials(GBasis(LTIdeal));
        RingHom phi = PolyAlgebraHom(myPolyRing, myPolyRing, DoComputeImage(GensLTIdeal, LTIdeal));
        myApplyHomToVector(NewGenerators.begin(), NewGenerators.end(), phi);
      }
      vector<vector<RingElem> > result;
      JBMill::Builder builder;
      builder.setStrategy(myJanetStrategy);
      JBMill mill(builder.setInput(NewGenerators));
      myContainer.myChangeBasis(mill.myGetJanetContainer().myListBegin(),
                                mill.myGetJanetContainer().myListEnd());
    }

    vector<RingElem> StableLTI::DoComputeImage(const vector<RingElem>& GensLTIdeal, const ideal& LTIdeal) const
    {
      const vector<RingElem>& x(indets(myPolyRing));
      vector<RingElem> images(x);
      for (const RingElem& m: GensLTIdeal)
      {
        CoCoA_ASSERT(IsMonomial(m));
        const long cls(myCls(m));
        const RingElem g(m / x[cls]);
        for (long l = 0; l < cls; l++)
        {
          if (!IsElem(g * x[l], LTIdeal))
          {
            images[cls] += images[l];
            return images;
          }
        }
      }
      return vector<RingElem>();
    }

    vector<RingElem> StableLTI::myLeadingMonomials(const vector<RingElem>& v) const
    {
      vector<RingElem> res; res.reserve(v.size());
      for (const RingElem& g: v)
        res.push_back(monomial(myPolyRing, one(CoeffRing(myPolyRing)), LPP(g)));

      return res;
    }

    vector<RingElem> StableLTIAll::DoComputeImage(const vector<RingElem>& GensLTIdeal, const ideal& LTIdeal) const
    {
      const vector<RingElem>& x(indets(myPolyRing));
      vector<RingElem> images(x);
      vector<vector<bool> > done(NumIndets(myPolyRing), vector<bool>(NumIndets(myPolyRing), false));

      for (const RingElem& m: GensLTIdeal)
      {
        CoCoA_ASSERT(IsMonomial(m));
        const long cls(myCls(m));
        const RingElem g(m / x[cls]);
        for (long l = 0; l < cls; l++)
        {
          if (!done[cls][l] && !IsElem(g * x[l], LTIdeal))
          {
            images[cls] += images[l];
            done[cls][l]  = true;
          }
        }
      }

      return images;
    }

    vector<RingElem> StronglyStableLTIAll::DoComputeImage(const vector<RingElem>& GensLTIdeal, const ideal& LTIdeal) const
    {
      const vector<RingElem>& x(indets(myPolyRing));
      vector<RingElem> images(x);
      vector<vector<bool> > done(NumIndets(myPolyRing), vector<bool>(NumIndets(myPolyRing), false));
      for (const RingElem& m: GensLTIdeal)
      {
        CoCoA_ASSERT(IsMonomial(m));
        const long cls(myCls(m));
        RingElem g(myPolyRing);
        for (long k = 1; k <= cls; k++)
        {
          if (!IsDivisible(g, m, x[k])) continue;
          for (long l = 0; l < k; l++)
          {
            if (!done[k][l] && !IsElem(g * x[l], LTIdeal))
            {
              images[k] += images[l];
              done[k][l] = true;
            }
          }
        }
      }

      return images;
    }

    vector<RingElem> StronglyStableLTI::DoComputeImage(const vector<RingElem>& GensLTIdeal, const ideal& LTIdeal) const
    {
      const vector<RingElem>& x(indets(myPolyRing));
      vector<RingElem> images(x);

      for (const RingElem& m: GensLTIdeal)
      {
        CoCoA_ASSERT(IsMonomial(m));
        const long cls(myCls(m));
        RingElem g(myPolyRing);
        for (long k = 1; k <= cls; k++)
        {
          if (!IsDivisible(g, m, x[k])) continue;
          for (long l = 0; l < k; l++)
          {
            if (!IsElem(g * x[l], LTIdeal))
            {
              images[k] += images[l];
              return images;
            }
          }
        }
      }
      return vector<RingElem>();
    }

    void StabilityAlgorithm::myPrintJBMill(const JBMill& mill) const
    {
      if (StatLevel == Logging)
      {
        std::cout << "==============================================================" << std::endl;
        std::cout << "====================== Begin JBMill ==========================" << std::endl;
        std::cout <<  "mill.myReturnJB() = " << mill.myReturnJB() << std::endl;
        std::cout <<  "len(mill.myReturnJB()) = " << len(mill.myReturnJB()) << std::endl;
        std::cout << "====================== End JBMill ============================" << std::endl;
        std::cout << "==============================================================" << std::endl;
      }
    }

    void StabilityAlgorithm::myPrintTransformation(const std::vector<RingElem>& transformation) const
    {
      if (StatLevel == Logging)
      {
        std::cout << "==============================================================" << std::endl;
        std::cout << "======================== Begin Transformation ================" << std::endl;
        std::cout <<  "transformation = " << transformation << std::endl;
        std::cout << "======================== End Transformation ==================" << std::endl;
        std::cout << "==============================================================" << std::endl;
      }
    }
  } // end namespace Involutive
} // end namespace CoCoA
