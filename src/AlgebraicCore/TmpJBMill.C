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

#include "CoCoA/TmpJBMill.H"
#include "CoCoA/submodule.H"
#include "CoCoA/RingQQ.H"

namespace CoCoA {
  namespace Involutive {

    // default value for the involutive criteria
    const std::bitset<3> JBMill::Builder::defaultInvCriteria(std::bitset<3>(3));

    std::unique_ptr<JBAlgorithm> JBMill::Builder::myChooseStrategy() const
    {
      typedef std::unique_ptr<JBAlgorithm> ptr;
      switch(myStrategy)
      {
        case  TQDegree:     return ptr(new DegreeTQ(myPolyRing(), myInvolutiveCriteria));
        case  TQBlockHigh:  return ptr(new BlockTQ(myPolyRing(), myInvolutiveCriteria, TQPSets::High));
        case  TQBlockLow:   return ptr(new BlockTQ(myPolyRing(), myInvolutiveCriteria, TQPSets::Low));
        case  GBCompletion: return ptr(new CompletionGB(myPolyRing()));
      }
      CoCoA_THROW_ERROR("Unknown strategy", __FUNCTION__);
      return nullptr;
    }

    JanetContainer JBMill::Builder::myUseStrategy() const
    {
      const char* const FnName = "Involutive::JBMill::Builder";
//      myCheckEmptyInput();
      if (myInput.empty()) CoCoA_THROW_ERROR("Empty or zero input", FnName);
//      myCheckMixedRingElems();
      if (!HasUniqueOwner(myInput)) CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
      std::unique_ptr<JBAlgorithm> strategy(myChooseStrategy());
      strategy->myComputer(myInput.begin(), myInput.end());
      return strategy->myOutputResult();
    }

    // void JBMill::Builder::myCheckEmptyInput() const
    // {
    //   if (myInput.empty())
    //   {
    //     CoCoA_THROW_ERROR("Empty or zero input", "Involutive::JBMill::Builder");
    //   }
    // }

    void JBMill::Builder::myRemoveZeroInput()
    {
      std::vector<RingElem> input;
      swap(input, myInput);
//      myInput.clear();
//      for(std::vector<RingElem>::const_iterator iter = input.begin(); iter != input.end(); ++iter)
      for (const RingElem& g: input)
      {
        if (!IsZero(g))
          myInput.push_back(g);
      }
    }

    // This is just HasUniqueOwner...
    // void JBMill::Builder::myCheckMixedRingElems() const
    // {
    //   CoCoA_ASSERT(!myInput.empty());
    //   std::vector<RingElem>::const_iterator iter(myInput.begin());
    //   SparsePolyRing ring(owner(*iter));
    //   for(std::vector<RingElem>::const_iterator iter = myInput.begin(); iter != myInput.end(); ++iter)
    //   {
    //     if(ring != owner(*iter))
    //     {
    //       CoCoA_THROW_ERROR(ERR::MixedRings, "Involutive::JBMill::Builder");
    //     }
    //   }
    // }

    JBMill::JBMill(const JBMill::Builder& builder) : myBasis(builder.myUseStrategy()),
                                                     myPolyRing(builder.myPolyRing()),
                                                     myPPMValue(PPM(myPolyRing))
    {}


    std::vector<RingElem> JBMill::myReturnJB() const
    {
      std::vector<RingElem> output;
        for(std::list<JanetTriple>::const_iterator it = myBasis.myListBegin(); it != myBasis.myListEnd(); ++it)
        {
          output.push_back(it->myGetPol());
        }
      return output;
    }


    std::vector<RingElem> JBMill::myReturnGB() const
    {
      std::vector<RingElem> output;
      // an element is part of the reduced GB if anc = lpp
      for(std::list<JanetTriple>::const_iterator it = myBasis.myListBegin(); it != myBasis.myListEnd(); ++it)
      {
        if(LPP(it->myGetPol()) == it->myGetAnc())
        {
          output.push_back(it->myGetPol());
        }
      }
      return output;
    }


    void JBMill::myPrintMultVar() const
    {
      std::cout << "computation of multiplicative variables:" << std::endl;
      std::map<PPMonoidElem,std::vector<bool> > MultVars = myMultVars();
      myOutputVar(MultVars, true);
    }

    void JBMill::myPrintNonMultVar() const
    {
      std::cout << "computation of nonmultiplicative variables:" << std::endl;
      std::map<PPMonoidElem,std::vector<bool> > MultVars = myMultVars();
      myOutputVar(MultVars, false);
    }


    void JBMill::myOutputVar(const std::map<PPMonoidElem,std::vector<bool> >& MultVars, bool OutputMultVar) const
    {
//      for(std::map<PPMonoidElem,std::vector<bool> >::const_iterator iter(MultVars.begin()); iter != MultVars.end(); ++iter)
      for (const auto& entry: MultVars)
      {
        std::cout << "LT(pol) = " << entry.first << std::endl;
        if (OutputMultVar)
        {
          std::cout << "multiplicative variables: ";
        }
        else
        {
          std::cout << "nonmultiplicative variables: ";
        }
        bool FirstIteration(true);
        const long VecSize = len(entry.second);
        // printing index of every iter, which is equal to OutputMultVar
        for(long VecIndex = 0; VecIndex != VecSize; ++VecIndex)
        {
          if (entry.second[VecIndex] == OutputMultVar)
          {
            if (!FirstIteration) std::cout << ", ";
            FirstIteration = false;
            std::cout << VecIndex;
          }
        }
        std::cout << std::endl;
        std::cout << "---------------------------------------" << std::endl;
      }
    }


    std::map<PPMonoidElem,std::vector<bool> > JBMill::myMultVars() const
    {
      //initialization
      JanetIterator iter(myBasis.myGetTree());
      std::vector<int> CurrentNonMultVars;
      const long CountIndets = NumIndets(myPolyRing);
      std::map<PPMonoidElem,std::vector<bool> > MultVars;
      for(std::list<JanetTriple>::const_iterator VecIter(myBasis.myListBegin()); VecIter != myBasis.myListEnd(); ++VecIter)
      {
        MultVars.insert(std::pair<PPMonoidElem, std::vector<bool> >(LPP(VecIter->myGetPol()), std::vector<bool>(CountIndets, true)));
      }
      //computation
      myRekComputeMultVar(MultVars, iter, CurrentNonMultVars, 0);
      //output
      return MultVars;
    }

    std::map<PPMonoidElem, std::vector<bool> > JBMill::myNonMultVars() const
    {
      std::map<PPMonoidElem, std::vector<bool> > MultVars = myMultVars();
      //reversing multiplicative variables -> nonMultVars
//      for(std::map<PPMonoidElem, std::vector<bool> >::iterator i = MultVars.begin(); i != MultVars.end(); ++i)
      for (auto& entry: MultVars)
      {
        entry.second = myReverseBoolVec(entry.second);
      }
      return MultVars;
    }

    std::vector<bool> JBMill::myNonMultVarsOf(ConstRefRingElem elem) const
    {
      std::map<PPMonoidElem, std::vector<bool> > NonMultVars = myNonMultVars();
      std::map<PPMonoidElem, std::vector<bool> >::iterator i(NonMultVars.find(LPP(elem)));
      if(i != NonMultVars.end())
      {
        return i->second;
      }
      return std::vector<bool>(NumIndets(myPolyRing), true);
    }


    std::vector< std::pair<RingElem, std::vector<bool> > > JBMill::myNonMultVarsWithRingElem() const
    {
      std::map<PPMonoidElem, std::vector<bool> > MultVarsWithPPM = myNonMultVars();
      std::vector< std::pair<RingElem, std::vector<bool> > > result;
      if(IsOne(LPP(myBasis.myListBegin()->myGetPol())))
      {
        result.push_back(std::pair<RingElem, std::vector<bool> >(myBasis.myListBegin()->myGetPol(), (MultVarsWithPPM.begin())->second));
        return result;
      }
//      for(std::map<PPMonoidElem, std::vector<bool> >::iterator i = MultVarsWithPPM.begin(); i != MultVarsWithPPM.end(); ++i)
      for (const auto& entry: MultVarsWithPPM)
      {
        JanetTriple* gPtr(myBasis.myGetTree().myJDivisor(entry.first));
        result.push_back(std::pair<RingElem, std::vector<bool> >(gPtr->myGetPol(), entry.second));
      }
      return result;
    }

    std::vector< std::pair<RingElem, std::vector<bool> > > JBMill::myMultVarsWithRingElem() const
    {
      std::map<PPMonoidElem, std::vector<bool> > MultVarsWithPPM = myMultVars();
      std::vector< std::pair<RingElem, std::vector<bool> > > result;
      if(IsOne(LPP(myBasis.myListBegin()->myGetPol())))
      {
        result.push_back(std::pair<RingElem, std::vector<bool> >(myBasis.myListBegin()->myGetPol(), (MultVarsWithPPM.begin())->second));
        return result;
      }
//      for(std::map<PPMonoidElem, std::vector<bool> >::iterator i = MultVarsWithPPM.begin(); i != MultVarsWithPPM.end(); ++i)
      for (const auto& entry: MultVarsWithPPM)
      {
        JanetTriple* gPtr(myBasis.myGetTree().myJDivisor(entry.first));
        result.push_back(std::pair<RingElem, std::vector<bool> >(gPtr->myGetPol(), entry.second));
      }
      return result;
    }



    void JBMill::myRekComputeMultVar(std::map<PPMonoidElem, std::vector<bool> >& MultVars, JanetIterator iter, std::vector<int> CurrentNonMultVars, int CurVar) const
    {
      do //until highest node in degree direction
      {
        // copy the current non mult vars
        std::vector<int> CopyVars = CurrentNonMultVars;
        // if this isn't the highest degree in current var add var to CopyVars
        if(iter.myDisNextDeg() != 0)
        {
          CopyVars.push_back(CurVar);
        }
        // if there is a node in variable-direction call myRekComputeMultVar
        if(iter.myDisNextVar())
        {
          int distance(iter.myDisNextVar());
          JanetIterator TmpIter(iter);
          TmpIter.myNextVar();
          myRekComputeMultVar(MultVars, TmpIter, CopyVars, CurVar + distance);
        }
      }
      while(iter.myNextDeg());

      // highest degree node for current variable -> this variable is multiplicative
      if(!(iter.myDisNextVar()))
      {
        //everything is multiplicative except the CurrentNonMultVars
//        for(std::vector<int>::iterator VecIter(CurrentNonMultVars.begin()); VecIter != CurrentNonMultVars.end(); ++VecIter)
        for (int VarIndex: CurrentNonMultVars)
        {
          MultVars[LPP(iter.myGetPol())][VarIndex]= false;
        }
      }
    }

    std::vector<bool> JBMill::myPommaretMultVar(ConstRefPPMonoidElem pp) const
    {
      std::vector<bool> MultVars(NumIndets(myPolyRing), false);
      long cls(0);
      if (!IsOne(pp))
      {
        cls = myCls(pp);
      }
      for (int i = cls; i < NumIndets(myPolyRing); ++i)
      {
        MultVars[i] = true;
      }
      return MultVars;
    }

    bool JBMill::IamPommaretBasis() const
    {
//      std::map<PPMonoidElem,std::vector<bool> > MultVars = myMultVars();
//      for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter(MultVars.begin()); iter != MultVars.end(); ++iter)
      for (const auto& entry: myMultVars())
      {
//??? Same as if (entry.second != myPommaretMultVar(entry.first)) return false; ???
        const std::vector<bool> IPBPommaretMultVar(myPommaretMultVar(entry.first));
        const long VecSize(len(entry.second));
        for(long i = 0; i != VecSize; ++i)
        {
          if (entry.second[i] != IPBPommaretMultVar[i])
          {
            return false;
          }
        }
      }
      return true;
    }

    bool JBMill::IamHomogenous() const
    {
//      std::vector<RingElem> gb(myReturnGB());
//      for(std::vector<RingElem>::iterator iter = gb.begin(); iter != gb.end(); ++iter)
      for (const RingElem& g: myReturnGB())
      {
        if (!IsHomog(g)) return false;
      }
      return true;
    }

    bool JBMill::IamMonomialIdeal() const
    {
//      std::vector<RingElem> jb(myReturnJB());
//      for(std::vector<RingElem>::iterator iter = jb.begin(); iter != jb.end(); ++iter)
      for (const RingElem& g: myReturnJB())
      {
        if (!IsMonomial(g)) return false;
      }
      return true;
    }


    //TODO: Throw error if r has the wrong type!!!
    std::pair<std::map<PPMonoidElem, RingElem>, RingElem> JBMill::myStandardRepresentation(ConstRefRingElem r) const
    {
      std::map<PPMonoidElem, RingElem> StandardRep;
      //initializing reduction history -> nothing is reduced every factor must be zero
      for(std::list<JanetTriple>::const_iterator iter(myBasis.myListBegin()); iter != myBasis.myListEnd(); ++iter)
      {
        StandardRep.insert(std::pair<PPMonoidElem, RingElem>(LPP(iter->myGetPol()), zero(myPolyRing)));
      }
      RingElem h(r);
      RingElem res(myPolyRing);
      //testing if basis is one
    //  if(IsOne(LPP(myBasis.myListBegin()->myGetPol())))
      if(myPPMValue->myIsOne(raw(myPolyRing->myLPP(raw(myBasis.myListBegin()->myGetPol())))))
      {
    //    StandardRep[LPP(myBasis.myListBegin()->myGetPol())] = h;
        StandardRep[myPolyRing->myLPP(raw(myBasis.myListBegin()->myGetPol()))] = h;
        h = 0;
      }
      //perform involutive full reduction
    //  while(!IsZero(h))
      while(!myPolyRing->myIsZero(raw(h)))
      {
    //    JanetTriple* gPtr(myBasis.myGetTree().myJDivisor(LPP(h)));
        JanetTriple* gPtr(myBasis.myGetTree().myJDivisor(myPolyRing->myLPP(raw(h))));
        if(gPtr == nullptr)
        {
    //      res += monomial(myPolyRing,LC(h),LPP(h));
    //      h -= monomial(myPolyRing,LC(h),LPP(h));

          RingElem tmp(myPolyRing);
          myPolyRing->myMoveLMToFront(raw(tmp), raw(h));
          myPolyRing->myAppendClear(raw(res), raw(tmp));
          // myPolyRing->myAddMul(raw(res), raw(h), raw(one(myPolyRing)));
          // myPolyRing->myDeleteLM(raw(h));
        }
        else
        {
          //if we can reduce something add factor to the according element in the map
          RingElem factor(myPolyRing);
          // factor = LPP(h)/LPP(gPtr->myGetPol())
          myPolyRing->myDivLM(raw(factor), raw(h), raw(gPtr->myGetPol()));
          // h = h - factor * gPtr->myGetPol()
          myPolyRing->myAddMulLM(raw(h), raw(-factor), raw(gPtr->myGetPol()));
          // factor = factor + StandardRep.find(LPP(gPtr->myGetPol()))->second)
          myPolyRing->myAdd(raw(factor), raw(factor), raw(StandardRep.find(LPP(gPtr->myGetPol()))->second));
          // these lines are inefficient!
          StandardRep.erase(LPP(gPtr->myGetPol()));
          StandardRep.insert(std::pair<PPMonoidElem,RingElem>(LPP(gPtr->myGetPol()), factor));
        }
      }
      // store in std::pair, first part is the reduction history, second part is the rest
      std::pair<std::map<PPMonoidElem, RingElem>, RingElem> pair(StandardRep, res);
      return pair;
    }

    std::vector<std::pair<RingElem, RingElem> > JBMill::myStandardRepresentationWithoutRest(ConstRefRingElem r) const
    {
      std::vector<std::pair<RingElem, RingElem> > result;
      const std::map<PPMonoidElem, RingElem> StandardRep(myStandardRepresentation(r).first);
//      for(std::map<PPMonoidElem, RingElem>::iterator i = StandardRep.begin(); i != StandardRep.end(); ++i)
      for (const auto& entry: StandardRep)
      {
        JanetTriple* gPtr(myBasis.myGetTree().myJDivisor(entry.first));
        result.push_back(std::pair<RingElem, RingElem>(gPtr->myGetPol(), entry.second));
      }
      return result;
    }

    std::vector<RingElem> JBMill::myStandardRepresentationWithoutRestShort(ConstRefRingElem r) const
    {
      std::vector<RingElem> result;
      const std::map<PPMonoidElem, RingElem> StandardRep(myStandardRepresentation(r).first);
//      for(std::map<PPMonoidElem, RingElem>::iterator i = StandardRep.begin(); i != StandardRep.end(); ++i)
      for (const auto& entry: StandardRep)
      {
        result.push_back(entry.second);
      }
      return result;
    }


    void JBMill::myOutputStandardRepresentation(ConstRefRingElem r) const
    {
      const std::pair<std::map<PPMonoidElem, RingElem>, RingElem> representation = myStandardRepresentation(r);
      const auto& StandardRep = representation.first;
      const RingElem& res = representation.second;
      std::cout << "Involutive Standard Representation" << std::endl;
      std::cout << r << " =" << std::endl;
      std::cout << std::endl;
//      for(std::map<PPMonoidElem, RingElem>::iterator iter(StandardRep.begin()); iter != StandardRep.end(); ++iter)
      for (const auto& entry: StandardRep)
      {
        std::cout << entry.second <<  std::endl;
      }
      std::cout << "rest = "<< res << std::endl;
      std::cout << "----------------------------------" << std::endl;
    }


    RingElem JBMill::myJNormalForm(const RingElem& elem) const
    {
      RingElem h(elem);
      RingElem res(zero(myPolyRing));
      if(IsOne(LPP(myBasis.myListBegin()->myGetPol())))
      {
        h = 0;
      }
      while(!IsZero(h))
      {
        JanetTriple* gPtr(myBasis.myGetTree().myJDivisor(myPolyRing->myLPP(raw(h))));
        if(gPtr == nullptr)
        {
          res += monomial(myPolyRing, LC(h), LPP(h));
          h -= monomial(myPolyRing, LC(h), LPP(h));
        }
        else
        {
          myPolyRing->myReductionStep(raw(h), raw(gPtr->myGetPol()));
        }
      }
      if(IsFractionField(CoeffRing(myPolyRing)) && !IsZero(h))
      {
        myPolyRing->myDivByCoeff(raw(res),raw(LC(res)));
        res = ClearDenom(res);
      }
      return res;
    }


    RingElem JBMill::myBinLike(PolyRing ring, RingElem PolAbove, long IntBelow) const
    {
      RingElem result = one(ring);
      for (long j = 0; j < IntBelow; ++j)
      {
        result *= (PolAbove  - j)/(j + 1);
      }
      return result;
    }

    long JBMill::myCountTrues(const std::vector<bool>& vec) const
    {
      long result(0);
//      for(std::vector<bool>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter)
      for (bool b: vec)
      {
        if (b) ++result;
      }
      return result;
    }

/*
q_t = deg t
k_t = Anzahl Mult Vars
hp(q) = \summe_{t \in T}((q - q_t + k_t - 1)(q - q_t))

*/
    RingElem JBMill::myHilbertPol(ConstRefRingElem s) const
    {
      const PolyRing& ring = owner(s);
      RingElem result(ring);

      result = myBinLike(ring, s + NumIndets(myPolyRing) - 1, NumIndets(myPolyRing) - 1);

//      std::map<PPMonoidElem, std::vector<bool> > MultVar(myMultVars());

      // for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter = MultVar.begin(); iter != MultVar.end(); ++iter)
      // {
      //   long deg = StdDeg(iter->first);
      //   long NumMultVars = myCountTrues(iter->second);
      //   result -= myBinLike(ring, s - deg + NumMultVars - 1, NumMultVars - 1);
      // }
//      for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter = MultVar.begin(); iter != MultVar.end(); ++iter)
      for (const auto& entry: myMultVars())
      {
        const long deg = StdDeg(entry.first);
        const long NumMultVars = myCountTrues(entry.second);
        result -= myBinLike(ring, s - deg + NumMultVars - 1, NumMultVars - 1);
      }

      return result;
    }

    long JBMill::myDeg() const
    {
      long res(0);
      for(std::list<JanetTriple>::const_iterator it = myBasis.myListBegin(); it != myBasis.myListEnd(); ++it)
      {
        if (res < deg(LPP(it->myGetPol()))) {
          res = deg(LPP(it->myGetPol()));
        }
      }
      return res;
    }

    void JBMill::myHilbertFunc() const
    {
      ring Q = RingQQ();
      PolyRing HilbPolyRing = NewPolyRing(Q, symbols("t"));
      RingElem FuncExpression(myHilbertPol(indet(HilbPolyRing, 0)));
      long maxDeg = myDeg();
      for (int i = 0; i < maxDeg; ++i)
      {
        std::cout << "H(" << i << ") = " << myHilbertFunc(BigInt(i)) << std::endl;
      }
      std::cout << "H(t) = " << FuncExpression << "    for t >= " << maxDeg << std::endl;
    }

    BigInt JBMill::myHilbertFunc(const BigInt& m) const
    {
      BigInt k(NumIndets(myPolyRing));
      BigInt r = m + k - 1;
      BigInt s = k - 1;
      BigInt res;
      if((r < s) || (r < 0) || (s < 0))
      {
        if((r == -1) && (s == -1))
        {
          res = 1;
        }
      }
      else
      {
        res = binomial(r, s);
      }


//      std::map<PPMonoidElem, std::vector<bool> > MultVar(myMultVars());
//      for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter = MultVar.begin(); iter != MultVar.end(); ++iter)
      for (const auto& entry: myMultVars())
      {
        ///??? SLUG why BigInt and not long?  JAA 2020-02-18
        const BigInt deg = BigInt(StdDeg(entry.first));
        const BigInt NumMultVars = BigInt(myCountTrues(entry.second));
        const BigInt r = m - deg + NumMultVars - 1;
        const BigInt s = NumMultVars - 1;
        if((r < s) || (r < 0) || (s < 0))
        {
          if((r == -1) && (s == -1))
          {
            res -= 1;
          }
        }
        else
        {
          res -= binomial(r, s);
        }
      }
      return res;
    }

    RingElem JBMill::myHilbertSeries(ConstRefRingElem s) const
    {
      RingElem res(owner(s));
      if (!IsFractionField(owner(s)))
        CoCoA_THROW_ERROR(ERR::NotElemFrF, "s must be in a Fraction Field!!");

//        std::map<PPMonoidElem, std::vector<bool> > MultVar(myMultVars());
//        for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter = MultVar.begin(); iter != MultVar.end(); ++iter)
      for (const auto& entry: myMultVars())
      {
        const long deg = StdDeg(entry.first);
        const long NumMultVars = myCountTrues(entry.second);
        RingElem numerator(power(s,deg));
        RingElem denominator(power(1 - s, NumMultVars));
        res += numerator/denominator;
      }
      return res;
    }

    FreeModule JBMill::myMakeNewFreeModuleForSyz() const
    {
      std::vector<RingElem> jb = myReturnJB();
      return NewFreeModuleForSyz(jb);
    }

    FGModule JBMill::mySyzygy() const
    {
      std::vector<ModuleElem> GenSyz;
      FreeModule FModule = myMakeNewFreeModuleForSyz(); //generates the module
      std::map<PPMonoidElem, std::vector<bool> > listMultVars = myMultVars(); //computes the multiplicative variables for the janet basis
      std::vector<ModuleElem> GenFModule = gens(FModule); //generators of the module
      std::vector<ModuleElem>::iterator GenFModuleIter = GenFModule.begin();
      for(std::list<JanetTriple>::const_iterator ListIter = myBasis.myListBegin(); ListIter != myBasis.myListEnd(); ++ListIter) //iterates over all elements in the janet basis
      {
        std::vector<bool> MultVars = (listMultVars.find(LPP(ListIter->myGetPol())))->second; //multiplicative variables of the current element
        long VarPos(0); //variable position, used in the for-loop below
//        for(std::vector<bool>::iterator MultVar = MultVars.begin(); MultVar != MultVars.end(); ++MultVar) //iterates over all variables
        for (const bool MultVar: MultVars)
        {
          if (!MultVar) //tests if the variables is nonmultiplicative
          {
            std::map<PPMonoidElem, RingElem> StandardRep = myStandardRepresentation(indet(myPolyRing, VarPos) * ListIter->myGetPol()).first; //computes the standardrepresentation of nm-var * element
            ModuleElem s = indet(myPolyRing, VarPos) * (*GenFModuleIter);
            std::vector<ModuleElem>::const_iterator StandardRepModuleIter = GenFModule.begin();
            for(std::map<PPMonoidElem, RingElem>::const_iterator StandardRepIter = StandardRep.begin(); StandardRepIter != StandardRep.end(); ++StandardRepIter)
            {
              s = s - (StandardRepIter->second) * (*StandardRepModuleIter);
              ++StandardRepModuleIter;
            }
            GenSyz.push_back(s);
          }
          ++VarPos;
        }
        ++GenFModuleIter;
      }
      return submodule(FModule, GenSyz);
    }

    long JBMill::myDim() const
    {
      ring Q = RingQQ();
      PolyRing HilbPolyRing = NewPolyRing(Q, symbols("t"));
      RingElem HilbPol(myHilbertPol(indet(HilbPolyRing, 0)));
      if (IsZero(HilbPol))
      {
        return 0;
      }
      return StdDeg(HilbPol);
    }



    void JBMill::myComplementaryDecompositionRecPart(std::vector< std::pair<PPMonoidElem, std::vector<bool> > >& output, JanetIterator JIter) const
    {
      // recursive part: Compute the decomposition for the next variable
      JanetIterator DegCounter(JIter);
      do
      {
        JanetIterator TmpIter(DegCounter);
        if( 0 != TmpIter.myNextVar())
        {
          myComplementaryDecompositionRecPart(output, TmpIter);
        }
      } while(DegCounter.myNextDeg() != 0);
      // compute the multiplicative variables up to the current variable
      JanetIterator HighestNode = JIter.myGotoHighestNode(); //HighestNode == my
      std::vector<bool> multVars = myMultVar(HighestNode); //mult or nonmult?
      multVars[JIter.myCurrentVar()] = false;
      for(long iter = JIter.myCurrentVar() + 1; iter != NumIndets(myPolyRing); ++iter)
      {
        multVars[iter] = true;
      }
      // If H is the JB, then every (myCurrentVar() exponent times current monomial) belongs to the complementary decomposition,
      // where no element in the JB is in the same class than these element. If there is already a var deg in degree 0 no element
      // of this class belongs to the complementary decomposition (0 = JIter.myCurrentDeg())
      CoCoA_ASSERT(JIter.myCurrentDeg() == 0);
      PPMonoidElem pp = JIter.myGetMonomial();
      if(JIter.myDisNextVar() == 0)
      {
        JIter.myNextDeg();
        CoCoA_ASSERT(JIter.myCurrentDeg() != 0);
      }

      for (long exp = 0; exp < JIter.myCurrentDeg(); ++exp)
      {
        output.push_back(make_pair(pp * IndetPower(owner(pp), JIter.myCurrentVar(), exp), multVars));
      }
    }

    std::vector< std::pair<PPMonoidElem, std::vector<bool> > > JBMill::myComplementaryDecomposition() const
    {
      if(!IamMonomialIdeal())
      {
        CoCoA_THROW_ERROR("ideal isn't monomial", "myComplementaryDecomposition");
      }
      return myComplementaryDecompositionLeadingIdeal();
    }

    std::vector< std::pair<PPMonoidElem, std::vector<bool> > > JBMill::myComplementaryDecompositionPolynomial() const
    {
      return myComplementaryDecompositionLeadingIdeal();
    }

    std::vector< std::pair<PPMonoidElem, std::vector<bool> > > JBMill::myComplementaryDecompositionLeadingIdeal() const
    {
      std::vector< std::pair<PPMonoidElem, std::vector<bool> > > output;
      if(!IsOne(LPP(myBasis.myListBegin()->myGetPol())))
      {
        JanetIterator iter(myBasis.myGetTree());
        myComplementaryDecompositionRecPart(output, iter);
      }
      return output;
    }


    std::vector< std::pair<PPMonoidElem, std::vector<bool> > > JBMill::myStandardPairs() const
    {
      if(!IamMonomialIdeal())
      {
        CoCoA_THROW_ERROR("ideal isn't monomial", "myStandardPairs");
      }
//      std::vector< std::pair<PPMonoidElem, std::vector<bool> > > CompDecomp = myComplementaryDecomposition();
      std::vector< std::pair<PPMonoidElem, std::vector<bool> > > res;
//      for (std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator iter = CompDecomp.begin(); iter != CompDecomp.end(); ++iter)
      for (const auto& entry: myComplementaryDecomposition())
      {
        bool NBarEmpty = true;
        for(long i = 0; i != NumIndets(myPolyRing); ++i)
        {
          NBarEmpty = NBarEmpty && !(entry.second[i] && (0 != exponent(entry.first, i)));
        }
        if(NBarEmpty)
        {
          res.push_back(entry);
        }
        else
        {
          PPMonoid pp = PPM(myPolyRing);
          PPMonoidElem ResultPPMElem(one(pp));
          for(long i = 0; i != NumIndets(myPolyRing); ++i)
          {
            if(!(entry.second)[i])
            {
              ResultPPMElem *= IndetPower(pp, i, exponent(entry.first, i));
            }
          }
          res.push_back(std::pair<PPMonoidElem, std::vector<bool> >(ResultPPMElem, entry.second));
        }
      }
      return res;
    }

/*
    std::pair<RingHom, std::vector<bool> > JBMill::myNoetherNormalization() const
    {
      //construnction of decomposition
      std::vector< std::pair<PPMonoidElem, std::vector<bool> > > CompDecomp = myComplementaryDecompositionPolynomial();
      //construction of nu
      std::vector<bool> nu;
      if(len(CompDecomp) > 1)
      {
        std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator iter = CompDecomp.begin();
        std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator TmpIter(iter);
        ++iter;
        nu = myUnionBoolVectors(TmpIter->second, iter->second);
        ++iter;
        for (std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator i = iter; i != CompDecomp.end(); ++i)
        {
          nu = myUnionBoolVectors(nu, i->second);
        }
      }
      else
      {
        std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator iter = CompDecomp.begin();
        nu = iter->second;
      }
      //construction of d
      long d(nu.size());
      for (std::vector<bool>::reverse_iterator i = nu.rbegin(); i != nu.rend(); ++i)
      {
        if(*i == false)
        {
          break;
        }
        else
        {
          --d;
        }
      }
      if(myDim() == d)
      {
        return std::pair<RingHom, std::vector<bool> >(IdentityHom(myPolyRing), nu);
      }
      RingElem z = myGreatestMultVar(nu);
      std::vector<RingElem> x = indets(myPolyRing);
      //choose p
      RingElem p(myPolyRing);

      long AddZ(1);
      long ChoosePoly(0);
      while(true)
      {
        long TmpChoosePoly(ChoosePoly);
        for(std::list<JanetTriple>::const_iterator i = myBasis.myListBegin(); i != myBasis.myListEnd(); ++i)
        {
          if(IamMultipleOfMultVars(LPP(i->myGetPol()), nu) == true)
          {
            if(TmpChoosePoly == 0)
            {
              p = i->myGetPol();
              --TmpChoosePoly;
              break;
            }
            else
            {
              --TmpChoosePoly;
            }
          }
        }
        if((TmpChoosePoly != -1) && (ChoosePoly != 0))
        {
          CoCoA_THROW_ERROR("Not able to compute noether normalization", "myNoetherNormalization");
        }
        //construction of homomorphism
        std::vector<RingElem> xImages;
        long n = NumIndets(myPolyRing);
        if(IsZero(p) ==  true)
        {
          for (long i = 0; i < n; ++i)
          {
            xImages.push_back(zero(myPolyRing));
          }
          RingHom varphi = PolyAlgebraHom(myPolyRing, myPolyRing, xImages);
          return std::pair<RingHom, std::vector<bool> >(varphi, nu);
        }
        for (long i = 0; i < n; ++i)
        {
          if((indet(myPolyRing, i) == z) || exponent(LPP(p), i) == 0)
          {
            xImages.push_back(indet(myPolyRing, i));
          }
          else
          {
            RingElem poly(indet(myPolyRing, i) - z);
            for (long i = 1; i < AddZ; ++i)
            {
              poly = poly - z;
            }
            if(IsZero(poly) == true)
            {
              xImages.clear();
              for (long i = 0; i < n; ++i)
              {
                xImages.push_back(zero(myPolyRing));
              }
              RingHom varphi = PolyAlgebraHom(myPolyRing, myPolyRing, xImages);
            }
            xImages.push_back(poly);
          }
        }
        RingHom varphi = PolyAlgebraHom(myPolyRing, myPolyRing, xImages);
        //recursion step
        std::vector<RingElem> VarphiImage;
        for(std::list<JanetTriple>::const_iterator i = myBasis.myListBegin(); i != myBasis.myListEnd(); ++i)
        {
          VarphiImage.push_back(varphi(i->myGetPol()));
        }
        CoCoA_THROW_ERROR("FIXME", "JBMill::myNoetherNormalization");
        JBMill TmpMill(*this);
        // JBMill TmpMill = ExtendedJanetBasis(VarphiImage);
        std::pair<RingHom, std::vector<bool> > noether = TmpMill.myNoetherNormalization();
        //output
        if(IsZero(noether.first(indet(myPolyRing, 0))) == false)
        {
          return std::pair<RingHom, std::vector<bool> >(noether.first(varphi), noether.second);
        }
        else
        {
          if(AddZ == 10)
          {
            AddZ = 0;
            ++ChoosePoly;
          }
          else
          {
            ++AddZ;
          }
        }
      }
    }
*/

/*
    bool JBMill::IamMultipleOfMultVars(PPMonoidElem ppm, std::vector<bool> v) const
    {
      long n = NumIndets(myPolyRing);
      for (long i = 0; i < n; ++i)
      {
        if(exponent(ppm, i) > 0)
        {
          if(v[i] == false)
          {
            return false;
          }
        }
      }
      return true;
    }
*/

/*
    RingElem JBMill::myGreatestMultVar(std::vector<bool> v) const
    {
      long n = NumIndets(myPolyRing);
      long MaxVar(-1);
      for(long i = 0; i < n; ++i)
      {
        if(v[i] == true)
        {
          if(MaxVar != -1)
          {
            if(LPP(indet(myPolyRing, i)) > LPP(indet(myPolyRing, MaxVar)))
            {
              MaxVar = i;
            }
          }
          else
          {
            MaxVar = i;
          }
        }
      }
      return indet(myPolyRing, MaxVar);
    }
*/
/*
    std::vector<bool> JBMill::myUnionBoolVectors(std::vector<bool> v1, std::vector<bool> v2) const
    {
      CoCoA_ASSERT(v1.size() == v2.size());
      std::vector<bool> result;
      std::vector<bool>::iterator iter2 = v2.begin();
      for (std::vector<bool>::iterator iter1 = v1.begin(); iter1 != v1.end(); ++iter1)
      {
        if(*iter1 == true)
        {
          result.push_back(true);
        }
        else
        {
          if(*iter2 == true)
          {
            result.push_back(true);
          }
          else
          {
            result.push_back(false);
          }
        }
        ++iter2;
      }
      return result;
    }
*/

    std::vector<bool> JBMill::myMultVar(JanetIterator iter) const
    {
      std::map<PPMonoidElem, std::vector<bool> > MultVars = myMultVars();
      std::map<PPMonoidElem, std::vector<bool> >::iterator MapIter = MultVars.find(iter.myGetMonomial());
      return MapIter->second;
    }

    std::vector<bool> JBMill::myNonMultVar(JanetIterator iter) const
    {
      return myReverseBoolVec(myMultVar(iter));
    }

    std::vector<bool> JBMill::myReverseBoolVec(const std::vector<bool>& vec) const
    {
      std::vector<bool> res(vec);
      res.flip();
//      for(std::vector<bool>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter)
//      {
//        res.push_back(!(*iter));
//      }
      return res;
    }

    long JBMill::myCls(ConstRefPPMonoidElem elem) const
    {
      const std::vector<long> expv = exponents(elem);
      long counter = len(expv);
      for (std::vector<long>::const_reverse_iterator i = expv.rbegin(); i != expv.rend(); ++i)
      {
        --counter;
        if (*i != 0)
        {
          return counter;
        }
      }
      // if we leave the for loop the exponent vector is zero -> elem is constant
      // the class must be -1
      return -1;
    }

    long JBMill::myMinCls() const
    {
      long indets = NumIndets(myPolyRing);
      long min(indets);
      for(std::list<JanetTriple>::const_iterator i = myBasis.myListBegin(); i != myBasis.myListEnd(); ++i)
      {
        if(min > myCls(LPP(i->myGetPol())))
        {
          min = myCls(LPP(i->myGetPol()));
        }
      }
      return min;
    }

    long JBMill::myMaxCls() const
    {
      long indets = NumIndets(myPolyRing);
      long max(indets);
      for(std::list<JanetTriple>::const_iterator i = myBasis.myListBegin(); i != myBasis.myListEnd(); ++i)
      {
        if(max < myCls(LPP(i->myGetPol())))
        {
          max = myCls(LPP(i->myGetPol()));
        }
      }
      return max;
    }

    std::vector<RingElem> JBMill::myElementsWithClass(long InputCls) const
    {
      std::vector<RingElem> res;
      for(std::list<JanetTriple>::const_iterator i = myBasis.myListBegin(); i != myBasis.myListEnd(); ++i)
      {
        if(InputCls == myCls(LPP(i->myGetPol())))
        {
          res.push_back(i->myGetPol());
        }
      }
      return res;
    }

  } // end of namespace Involutive
} // end of namespace CoCoA
