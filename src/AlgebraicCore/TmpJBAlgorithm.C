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

#include "CoCoA/TmpJBAlgorithm.H"
#include "CoCoA/TmpGOperations.H"
#include "CoCoA/FractionField.H"


namespace CoCoA
{
  namespace Involutive
  {
    using std::make_pair;
    using std::pair;
    using std::vector;

    void TQAlgorithm::myTrivialIdeal() {
      myJTree.myDelete();
      myGetSets().myClearSets();
      JanetTriple trivalTriple(one(myPolyRing));
      myGetSets().myInsertSetT(trivalTriple);
    }

    void DegreeTQ::myComputer(std::vector<RingElem>::const_iterator beginInput,
                              std::vector<RingElem>::const_iterator endInput) {
      //initialization of janet-triples
      myInitialization(beginInput, endInput);

      std::list<JanetTriple>::iterator iter(mySets.myMinimizeAndInsertSetT(myJTree));
      //exit if there is a constant input
      // ??????
      if (IsConstant(iter->myGetPol())) {
        myTrivialIdeal();
        return;
      }
      // adding first element to janet tree
      JanetTree tree(JanetTree::OneLeafTree(myPolyRing, iter, 0, one(myPPMValue)));
      myJTree.myAddAtBegin(tree);
      //main loop
      while (mySets.myHeadReduceSetQ(myJTree)) {
        long SizeT(mySets.mySizeSetT());
        //insert element to setT and minimize setT
        std::list<JanetTriple>::iterator iter(mySets.myMinimizeAndInsertSetT(myJTree));
        // exit if we get a constant element
        if (IsConstant(iter->myGetPol())) {
          myTrivialIdeal();
          return;
        }
        //if we 'lost' an element in setT rebuild janet tree
        //otherwise add new element to setT and tailreduce setT
        if (mySets.mySizeSetT() <= SizeT) {
          myJTree.myDelete();
          for (TQSets::MultisetIterator IterT(mySets.myBeginSetT()); IterT != mySets.myEndSetT(); ++IterT) {
            std::list<JanetTriple>::iterator ListIter(*IterT);
            mySets.myJTailNormalForm(myJTree, ListIter);
            mySets.myJInsert(ListIter, myJTree);
          }
        } else {
          mySets.myJInsert(iter, myJTree);
          mySets.myTailReduceSetT(myJTree, iter);
        }
      }
    }

    JanetContainer TQAlgorithm::myOutputResult()
    {
      std::list<JanetTriple> jBasis;
      for (TQSets::MultisetIterator IterT(myGetSets().myBeginSetT()); IterT != myGetSets().myEndSetT(); ++IterT) {
        jBasis.push_back(*(*IterT));
      }
      return JanetContainer(myPolyRing, jBasis);
    }

    void TQAlgorithm::myInitialization(std::vector<RingElem>::const_iterator beginInput,
                                       std::vector<RingElem>::const_iterator endInput) {
      //initialization
      for (std::vector<RingElem>::const_iterator iter = beginInput; iter != endInput; ++iter) {
        if (IsZero(*iter)) {
          continue;
        }
        RingElem elem(*iter);
        if (IsFractionField(CoeffRing(myPolyRing))) {
          //      myPolyRing->myDivByCoeff(raw(*iter),raw(LC(*iter)));
          myPolyRing->myDivByCoeff(raw(elem), raw(LC(elem)));
          elem = ClearDenom(elem);
        }
        JanetTriple triple(elem, myPolyRing->myLPP(raw(elem)));
        myGetSets().myInsertSetQ(triple);
      }
    }

    void BlockTQ::myComputer(std::vector<RingElem>::const_iterator beginInput,
                             std::vector<RingElem>::const_iterator endInput) {
      //initialization
      myInitialization(beginInput, endInput);
      std::list<JanetTriple>::iterator iter(mySets.myMinimizeAndInsertSetT(myJTree));
      //return if constant input
      if (IsConstant(iter->myGetPol())) {
        myTrivialIdeal();
        return;
      }
      JanetTree tree(JanetTree::OneLeafTree(myPolyRing, iter, 0, one(myPPMValue)));
      myJTree.myAddAtBegin(tree);
      //main loop
      while (!mySets.IamEmptySetQ()) {
        //construct mySetP
        mySets.myMoveFromQtoP(myJTree);
        // update mySetP
        bool ModifySetT = mySets.IamJUpdate();
        // if mySetP empty continue
        if (mySets.myBeginSetP() == mySets.myEndSetP()) {
          continue;
        }
        if (ModifySetT) {
          //recompute janettree
          long SizeT(mySets.mySizeSetT());
          mySets.myMinimizeSetT();
          if (mySets.mySizeSetT() < SizeT) {
            myJTree.myDelete();
            for (TQSets::MultisetIterator IterT(mySets.myBeginSetT()); IterT != mySets.myEndSetT(); ++IterT) {
              std::list<JanetTriple>::iterator ListIter(*IterT);
              mySets.myJInsert(ListIter, myJTree);
            }
          }
          // mySets.myTailReduceSetTAll(myJTree);
        }
        TQSets::MultisetIterator IterP = mySets.myEndSetP();
        //constant input -> retun
        if (IsConstant((*(mySets.myBeginSetP()))->myGetPol())) {
          myTrivialIdeal();
          return;
        }
        // insert elements from mySetP into janet tree
        do {
          --IterP;
          std::list<JanetTriple>::iterator ListIter(*IterP);
          // mySets.myJTailNormalForm(myJTree, ListIter);
          mySets.myJInsert(ListIter, myJTree);
        } while (IterP != mySets.myBeginSetP());
        //tail reduce set T if necessary
        // if (!ModifySetT) {
        //   mySets.myTailReduceSetTAll(myJTree);
        // }
        mySets.myInsertSetPInSetT();
        mySets.myTailReduceSetTAll(myJTree);
      }
    }

    void CompletionGB::myBuildInitialTree(std::vector<RingElem>::const_iterator beginInput,
                                          std::vector<RingElem>::const_iterator endInput) {
      std::vector<RingElem> inputVec(beginInput, endInput);
      std::vector<RingElem> groebnerBasis;
      std::vector<RingElem> minGens;
      // computation of GBasis
      ComputeGBasis(groebnerBasis, minGens, inputVec);
      // if GBasis is one our ideal is trivial an we can skip everything else.
      if(len(groebnerBasis) == 1 && IsConstant(groebnerBasis[0]))
      {
        myJBasis.push_front(JanetTriple(one(myPolyRing)));
        return;
      }
      // add GBasis to JanetTree and JBasis
//      for(std::vector<RingElem>::iterator i = groebnerBasis.begin(); i != groebnerBasis.end(); ++i)
      for (const RingElem& g: groebnerBasis)
      {
        myJBasis.push_front(JanetTriple(g));
        myJTree.myInsert(myJBasis.begin());
      }
    }

    bool CompletionGB::myElemIsNonMultForIndex(const PPMonoidElem& m, long index)
    {
      CoCoA_ASSERT(NumIndets(myPolyRing) > index);
      CoCoA_ASSERT(index >= 0);
      JanetIterator iter(myJTree);
      const std::vector<long> expv = exponents(m);

      while(iter.myCurrentVar() < index) {
        while(iter.myCurrentDeg() < expv[iter.myCurrentVar()])
        {
          long degDis(iter.myNextDeg());
          (void) degDis;
          CoCoA_ASSERT(degDis > 0);
        }
        CoCoA_ASSERT(iter.myCurrentDeg() == expv[iter.myCurrentVar()]);
        long varDis(iter.myNextVar());
        if (varDis == 0) {
          return true;
        }
        CoCoA_ASSERT(varDis > 0);
      }
      if(iter.myCurrentVar() > index) {
        CoCoA_ASSERT(expv[index] == 0);
        return true;
      }
      while(iter.myNextDeg() != 0) { /* do nothing */ }
      CoCoA_ASSERT(iter.myCurrentVar() == index);
      CoCoA_ASSERT(iter.myCurrentDeg() >= expv[index]);
      return iter.myCurrentDeg() != expv[index];
    }

    void CompletionGB::myComputer(std::vector<RingElem>::const_iterator beginInput, std::vector<RingElem>::const_iterator endInput) {
      // first build initial JTree out of GBasis
      myBuildInitialTree(beginInput, endInput);
      // skip everything if the ideal is trivial
      if(len(myJBasis) == 1 && IsConstant(myJBasis.front().myGetPol()))
      {
        return;
      }
      // iterating over number of vars
      for(long i = 0; i != NumIndets(myPolyRing); ++i)
      {
        elemPairList prolongSet;
        // identify elements for which i is nonmultiplicative
        // if it is the case add a prolongation
//        for(std::list<JanetTriple>::iterator iter = myJBasis.begin(); iter != myJBasis.end(); ++iter)
        for (const JanetTriple& JT: myJBasis)
        {
          if (myElemIsNonMultForIndex(LPP(JT.myGetPol()), i))
            myPushToProlongList(prolongSet, JT.myGetPol() * indet(myPolyRing, i), LPP(JT.myGetPol()));
        }

        // sort the prolongSet after LPP of elemPair.first
        prolongSet.sort(myLPPComparison);
        // repeat until prolongSet is empty
        while(!prolongSet.empty()) {
          elemPair pair(prolongSet.front());
          // if pair has no involutive divisor add this element to JBasis
          // after that remove the current element from the prolongSet
          // and add new ones which arises from the newly inserted element
          if(myJTree.myJDivisor(LPP(pair.first)) == 0)
          {
            myJBasis.push_front(JanetTriple(pair.first, pair.second));
            // std::cout <<  "new_jb_elem = " << pair.first << std::endl;
            // if the new element would be not a leaf node we have to remove
            // everything behind it.
            if(myJTree.IamInternalNode(LPP(pair.first))) {
              std::vector<nodeData> removedNodes(myJTree.myInsertInTreeStructure(myJBasis.begin()));
//              for (std::vector<nodeData>::iterator remove = removedNodes.begin(); remove != removedNodes.end(); ++remove) {
              for (const nodeData& remove: removedNodes)
                myJBasis.erase(remove);
            } else {
              myJTree.myInsert(myJBasis.begin());
            }
            prolongSet.pop_front();
            if(myElemIsNonMultForIndex(LPP(pair.first), i)) {
              myInsertToProlongList(prolongSet, pair.first * indet(myPolyRing, i), pair.second);
            }
          } else {
            prolongSet.pop_front();
          }
        }
      }
      mySortAndMinimizeJanetBasis();
    }

    void CompletionGB::mySortAndMinimizeJanetBasis()
    {
      myJTree.myDelete();
      std::list<JanetTriple> TmpBasis;
      TmpBasis.swap(myJBasis);
      std::list<JanetTriple>::iterator iter = TmpBasis.end();
      --iter;
      // build initial tree out of Groebner basis elements
      // They are at the end of the list & the ancestor is equal
      // to the leading polynomial
      while (!TmpBasis.empty() && LPP(iter->myGetPol()) == iter->myGetAnc())
      {
        std::list<JanetTriple>::iterator TmpIter(iter);
        --iter;
        myJBasis.splice(myJBasis.begin(), TmpBasis, TmpIter);
        myJTree.myInsert(myJBasis.begin());
      }
      // sort rest of TmpBasis by degree
      TmpBasis.sort(myLPPTripleComparison);
      // Insert the rest of TmpBasis to Janet basis:
      // We take an element out of the sorted list and check j-divisibility.
      // If it is j-divisible we skip this element. If it is not j-divisible
      // we insert it to the Janet basis and return to the beginning of the list
      // to ensure minimality of the Janet basis.
      bool EverythingReducedToZero = true;
      do {
        std::list<JanetTriple>::iterator iter = TmpBasis.begin();
        EverythingReducedToZero = true;
        while (iter != TmpBasis.end())
        {
          std::list<JanetTriple>::iterator TmpIter(iter);
          ++iter;
          JanetTriple* divisor(myJTree.myJDivisor(LPP(TmpIter->myGetPol())));
          if(divisor == 0)
          {
            EverythingReducedToZero = false;
            myJBasis.splice(myJBasis.begin(), TmpBasis, TmpIter);
            myJTree.myInsert(myJBasis.begin());
            myJTree.myJTailNormalForm(myJBasis.begin()->myGetPol());
            break;
          }
        }
      } while (!EverythingReducedToZero);
      // perform tail normal form computation
//      for (std::list<JanetTriple>::iterator i = myJBasis.begin(); i != myJBasis.end(); ++i)
      for (JanetTriple& JT: myJBasis)
      {
        myJTree.myJTailNormalForm(JT.myGetPol());
      }
      // sort myJBasis
      myJBasis.sort(myLPPTripleComparison);
    }

    void CompletionGB::myPushToProlongList(elemPairList& list, const RingElem& elem, const PPMonoidElem& anc)
    {
      list.push_back(make_pair(elem, anc));
    }

    void CompletionGB::myInsertToProlongList(elemPairList& list, const RingElem& elem, const PPMonoidElem& anc)
    {
      elemPair pair(elem, anc);
      elemPairIter posIter(lower_bound(list.begin(), list.end(), pair, myLPPComparison));
      list.insert(posIter, pair);
    }

    JanetContainer CompletionGB::myOutputResult() {
      return JanetContainer(myPolyRing, myJBasis);
    }

  } // end of namespace Involutive
} // end of namespace CoCoA
