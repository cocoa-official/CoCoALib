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

#include "CoCoA/TmpJBDatastructure.H"
#include "CoCoA/ReductionCog.H"
#include "CoCoA/FractionField.H"

namespace CoCoA {
  namespace Involutive {

    void JanetTree::myAddAtBegin(JanetTree& tree) {
      if (myBeginVar == tree.myBeginVar) {
        //if the current tree has the same beginning var as the new tree, then we can
        //connect these trees easily
        JanetIterator iter(*this);
        iter.myConnectJanetTreeDeg(tree);
      } else {
        if (tree.myBeginDeg == 0) {
          //dangerours if tree.myBeginVar < myBeginVar we have a negative var_dis in the tree
          (*myArm.begin())->mySetDisNextVar(tree.myBeginVar - myBeginVar);
          (*myArm.begin())->mySetNextVarArm(*tree.myGetRootArm());
        } else {
          JanetIterator iter(*this);
          //create new var node such that this var matches the beginVar of tree
          iter.mySetNextVar(tree.myBeginVar - myBeginVar);
          iter.myNextVar();
          //then connect
          iter.myConnectJanetTreeDeg(tree);
        }
      }
    }

    void JanetIterator::myChangeToLeafNode(const nodeData triple) {
      (*myCurIter) = JanetHandle(new JanetLeafNodeImpl(triple));
      //if there is a arm which hanging around we must delete this arm:
      myCurArm->erase(++myCurIter, myCurArm->end());
      myCurIter = myCurArm->end();
      --myCurIter;
    }

    long JanetIterator::myNextDeg() {
      // check distance to next degree
      long dis = (*myCurIter)->myGetDisNextDeg();
      // if distance is not 0 then we add the distance to the current distance
      // and moving the vector iterator
      if (dis) {
        myMonomial[myCurVar] += dis;
        ++myCurIter;
        return dis;
      }
      // if not there is no moving -> moving distance is zero
      return 0;
    }

    long JanetIterator::myNextVar() {
      // check distance to next var
      long dis((*myCurIter)->myGetDisNextVar());
      // if distance is not 0 then we 'initilize the degree at the current var'
      // and jumping to the next arm
      if (dis) {
        myCurVar += dis;
        CoCoA_ASSERT(myCurVar < NumIndets(myCurTree->myGetPolyRing()));
        myMonomial[myCurVar] = 0;
        myCurArm = (*myCurIter)->myNextVarArm();
        myCurIter = myCurArm->begin();
        return dis;
      }
      return 0;
    }

    void JanetIterator::mySetNextDeg(long dis) {
      //nothing to do if dis == 0
      if (dis > 0) {
        //dis to next node
#ifdef CoCoA_DEBUG
        const std::list<JanetHandle>::iterator curIter(myCurIter);
#endif
        long deg = (*myCurIter)->myGetDisNextDeg();
        if (deg) {
          // it there is a next deg we have to check if the new deg would be between
          // the current deg and the following one
          if (deg > dis) {
            // if this is the case we add the new node between both nodes
            myCurArm->insert(++myCurIter, JanetHandle(new JanetInternalNodeImpl()));
            (*(--myCurIter))->mySetDisNextDeg(deg - dis);
            (*(--myCurIter))->mySetDisNextDeg(dis);
          } else {
            // if not we erase the deg-direction starting with the node after current node
            // (in deg direction)
            myCurArm->erase((++myCurIter), myCurArm->end());
            myCurIter = myCurArm->end();
            (*(--myCurIter))->mySetDisNextDeg(dis);
            myCurArm->push_back(JanetHandle(new JanetInternalNodeImpl()));
          }
        } else {
          //if there is no next deg we can append the new node easily
          if ((*myCurIter)->IamLeafNode()) {
            *myCurIter = JanetHandle(new JanetInternalNodeImpl());
          }
          (*myCurIter)->mySetDisNextDeg(dis);
          myCurArm->push_back(JanetHandle(new JanetInternalNodeImpl()));
        }
        CoCoA_ASSERT((*myCurIter)->myGetDisNextDeg() == dis);
        CoCoA_ASSERT(myCurIter == curIter);
      }
    }

    void JanetIterator::mySetNextVar(long dis) {
      //same as above but with var.
      if (dis) {
        long var((*myCurIter)->myGetDisNextVar());
        if (var) {
          if (var > dis) {
            JanetHandle NodePtr(new JanetInternalNodeImpl());
            // JanetInternalNodeImpl node;
            std::list<JanetHandle>* varArm((*myCurIter)->myNextVarArm());
            NodePtr->mySetNextVarArm(*varArm);
            // node.mySetNextVarArm(*varArm);
            NodePtr->mySetDisNextVar(var - dis);
            // node.mySetDisNextVar(var - dis);
            varArm->erase(varArm->begin(), varArm->end());
            varArm->push_back(NodePtr);
            (*myCurIter)->mySetDisNextVar(dis);
          } else {
            (*myCurIter)->mySetDisNextVar(dis);
            std::list<JanetHandle>* varArm((*myCurIter)->myNextVarArm());
            varArm->erase(varArm->begin(), varArm->end());
            varArm->push_back(JanetHandle(new JanetInternalNodeImpl()));
          }
        } else {
          if ((*myCurIter)->IamLeafNode()) {
            *myCurIter = JanetHandle(new JanetInternalNodeImpl());
          }
          (*myCurIter)->mySetDisNextVar(dis);
          std::list<JanetHandle>* varArm((*myCurIter)->myNextVarArm());
          varArm->push_back(JanetHandle(new JanetInternalNodeImpl()));
        }
      }
    }

    void JanetIterator::myConnectJanetTreeDeg(JanetTree &tree) {
      // check if current position and the beginning of the tree matches
      if ((tree.myGetBeginDeg() == myMonomial[myCurVar]) && (tree.myGetBeginVar() == myCurVar)) {
        // if yes everything is nice and we can combine the trees (we combine the last node
        // of the current tree and the first node of tree)
        // maybe here is a mistake because we don't check whether the current node is useless
        if ((*myCurIter)->IamLeafNode()) {
          (*myCurIter) = JanetHandle(new JanetInternalNodeImpl());
        }
        std::list<JanetHandle>* rootArm(tree.myGetRootArm());
        std::list<JanetHandle>::iterator tempIter(rootArm->begin());
        (*myCurIter)->mySetNextVarArm(*((*tempIter)->myNextVarArm()));
        (*myCurIter)->mySetDisNextDeg((*tempIter)->myGetDisNextDeg());
        std::list<JanetHandle>::iterator beg(rootArm->begin()), end(rootArm->end());
        for (std::list<JanetHandle>::iterator it = ++beg; it != end; ++it) {
          myCurArm->push_back(*it);
        }
      } else {
        // also simpling combining. But no we have a distande between the curNode and the beginning of
        // the new tree
        long extraDis(0);
        std::list<JanetHandle>* arm(tree.myGetRootArm());
        std::list<JanetHandle>::iterator beg(arm->begin()), end(arm->end());
        // checks if there is at the beginning a node in var direction and if the beginning node is a leaf node
        // if this isn't the case then we ignore the beginning node (to provide the minimality of the tree)
        if (!((*beg)->myGetDisNextVar()) && !(*beg)->IamLeafNode()) {
          extraDis = (*beg)->myGetDisNextDeg();
          ++beg;
        }
        (*myCurIter)->mySetDisNextDeg(tree.myGetBeginDeg() + extraDis - myMonomial[myCurVar]);
        myCurArm->erase(++myCurIter, myCurArm->end());
        myCurIter = myCurArm->end();
        --myCurIter;
        for (std::list<JanetHandle>::iterator it = beg; it != end; ++it) {
          myCurArm->push_back(*it);
        }
      }
    }

    bool JanetIterator::IamSetNM(const long &index) {
      // if this isn't a leaf node return false
      if ((*myCurIter)->IamLeafNode()) {
        JanetLeafNodeImpl* node(dynamic_cast<JanetLeafNodeImpl*>(myCurIter->get()));
        CoCoA_ASSERT(node != 0);
        node->myGetTriple().myVarIsProlonged(index);
        return true;
      }
      return false;
    }

    void JanetIterator::myReturnToBegin() {
      // everything to the beginning
      myCurArm = myCurTree->myGetRootArm();
      myCurIter = myCurArm->begin();
      myCurVar = myCurTree->myGetBeginVar();
      myMonomial = std::vector<long>(NumIndets(myCurTree->myGetPolyRing()));
      myMonomial[myCurVar] = myCurTree->myGetBeginDeg();
    }

    long JanetIterator::myPrevDeg() {
      // returns to the previous deg. If we are already at the beginning do nothing
      if (myCurIter == myCurArm->begin()) {
        return 0;
      } else {
        --myCurIter;
        myMonomial[myCurVar] -= (*myCurIter)->myGetDisNextDeg();
        return (*myCurIter)->myGetDisNextDeg();
      }
    }

    void JanetTriple::myClearProlongedVars() {
      // delete all nonmultiplicative variables -> all entries false
//      for (std::vector<bool>::iterator iter = myAlreadyProlongedVars.begin(); iter != myAlreadyProlongedVars.end(); ++iter) {
//        *iter = false;
      // JAA the following should be functionally equivalent & faster (but possibly less clear?)
      const long n = len(myAlreadyProlongedVars);
      myAlreadyProlongedVars.clear();
      myAlreadyProlongedVars.resize(n);
  }

    // TODOJOHN Put in cost
    JanetTriple* JanetTree::myJDivisor(ConstRefPPMonoidElem w) const {
      // straight forward implementation
      std::vector<long> expv(NumIndets(PPM(myPolyRing)));
      PPM(myPolyRing)->myExponents(expv, raw(w));
      std::list<JanetHandle>::const_iterator iter(myArm.begin());
      long curVar(myBeginVar);
      long curDeg(myBeginDeg);
      long numIndets(NumIndets(owner(w)));
      while (numIndets >= curVar) {
        while (((*iter)->myGetDisNextDeg()) && (curDeg < expv[curVar]))    //exponent(w,curVar)))
        {
          if (curDeg + (*iter)->myGetDisNextDeg() > expv[curVar]) {
            return 0;
          }
          curDeg += (*iter)->myGetDisNextDeg();
          ++iter;
        }
        if ((*iter)->myGetDisNextVar()) {
          curVar += (*iter)->myGetDisNextVar();
          iter = ((*iter)->myNextVarArm())->begin();
          curDeg = 0;
        } else {
          if ((*iter)->myGetDisNextDeg()) {
            return 0;
          } else {
            JanetLeafNodeImpl* node(dynamic_cast<JanetLeafNodeImpl*>(iter->get()));
            if (node == 0) {
              return 0;
            }
            CoCoA_ASSERT(node != 0);
            return node->myGetTriplePtr();
          }
        }
      }
      return 0;  //only of safety. A return command going to be used in the while loop
    }

    JanetIterator JanetIterator::myGotoHighestNode() {
      // goto the highest node: If we are in the highest node go to the next var if possible
      // if not return the current node
      JanetIterator iter(*this);
      while (true) {
        while (iter.myNextDeg() != 0) {
        }
        if (iter.myDisNextVar() == 0) {
          return iter;
        } else {
          iter.myNextVar();
        }
      }
      return iter;
    }

    JanetTree JanetTree::OneLeafTree(const SparsePolyRing& polyRing, const std::list<JanetTriple>::iterator& ps, long i,
                                     PPMonoidElem w) {
    // if we say ps we mean myPolyRing->myLPP(raw(ps->myGetPol()))
      long n = NumIndets(polyRing);
    //navigates to the first index where deg_i(ps) != 0
      while ((i < n - 1) && (exponent(polyRing->myLPP(raw(ps->myGetPol())), i) == 0))  //we don't reach n because the algorithms which use myOneLeafTree check this
      {
        ++i;
      }
    //janet tree starts at deg_i(w) at variable_index i
      JanetTree LeafTree(polyRing, exponent(w, i), i);
      JanetIterator iter(LeafTree);
      JanetIterator TempIter(iter);
      PPMonoidElem u(polyRing->myLPP(raw(ps->myGetPol())));
      while ((u != w) && (i < n)) {
        // add node in deg-direction
        if (exponent(w, i) < exponent(u, i)) {
          iter.mySetNextDeg(exponent(u, i) - exponent(w, i));
          w *= IndetPower(PPM(polyRing), i, exponent(u, i) - exponent(w, i));
          iter.myNextDeg();
          TempIter = iter;
        }
        // add node in var-direction
        if (u != w && (iter.myCurrentVar() < (n - 1))) {
          if (i == n - 1) {
            //should raise an error???????
            return JanetTree(polyRing, exponent(w, i), i);
          }
          // let i = iter.myCurrentVar()
          // deg_i(ps) = 0 -> skipping this variable in the janet tree
          if ((TempIter.myCurrentVar()) != (iter.myCurrentVar())) {
            TempIter.mySetNextVar(iter.myCurrentVar() - TempIter.myCurrentVar() + 1);
            iter = TempIter;
          } else {
            iter.mySetNextVar(1);
          }
          ++i;
          iter.myNextVar();
        }
      }
    // the node was a internal node, but it must be a leaf node
      iter.myChangeToLeafNode(ps);
      return LeafTree;
    }

    void JanetTree::myInsert(const nodeData& ps) {
      //like myJInsert skipping prolongations
      long i(0);
      JanetIterator iter(*this);
      long n(NumIndets(myPolyRing));
      while (i < n) {
        while ((iter.myCurrentDeg() < exponent(LPP(*(ps->myGetPolPtr())), i)) && (iter.myDisNextDeg())) {
          iter.myNextDeg();
        }
        if (iter.myCurrentDeg() < exponent(LPP(*(ps->myGetPolPtr())), i)) {
          JanetTree tree(OneLeafTree(myPolyRing, ps, i, iter.myGetMonomial() * indet(PPM(myPolyRing), i)));
          iter.myConnectJanetTreeDeg(tree);
          break;
        } else {
          if (exponent(LPP(*(ps->myGetPolPtr())), i) == iter.myCurrentDeg() && (iter.myDisNextVar())) {
            bool NewVarNode(false);
            long DisNextVarNode(0);
            for (long j = i + 1; j < (iter.myDisNextVar() + iter.myCurrentVar()); ++j) {
              NewVarNode = (exponent(LPP(*(ps->myGetPolPtr())), j) != 0);
              if (NewVarNode) {
                DisNextVarNode = j - i;
                break;
              }
            }
            if (NewVarNode) {
              iter.mySetNextVar(DisNextVarNode);
            }
            i += iter.myNextVar();
          } else  //exponent(LPP(ps.myGetPol()),i) > iter.myCurrentDeg()
          {
            if (exponent(LPP(*(ps->myGetPolPtr())), i) < iter.myCurrentDeg()) {
              iter.myPrevDeg();
              iter.mySetNextDeg(exponent(LPP(*(ps->myGetPolPtr())), i) - iter.myCurrentDeg());
              iter.myNextDeg();
            }
            JanetTree tree(OneLeafTree(myPolyRing, ps, i + 1, iter.myGetMonomial()));
            iter.mySetNextVar(tree.myGetBeginVar() - iter.myCurrentVar());
            iter.myNextVar();
            iter.myConnectJanetTreeDeg(tree);
            break;
          }
        }
      }
    }

    void JanetTree::myJFullNormalForm(RingElem& elem) {
      myJHeadNormalForm(elem);
      myJTailNormalForm(elem);
    }

    void JanetTree::myJTailNormalForm(RingElem& elem) {
     //  TODO Reduction Cog must be a data member of JanetTree
      if(IsZero(elem)) {
        return;
      }
      ReductionCog myReductionCog(NewRedCogPolyField(myPolyRing));
      //skipping the leading-monomial (we want to perform a !tail!-reduction)
      myReductionCog->myAssignReset(elem);
      myReductionCog->myMoveToNextLM();
      while (!IsActiveZero(myReductionCog)) {
        JanetTriple* g = myJDivisor(ActiveLPP(myReductionCog));
        if (g != 0) {
          myReductionCog->myReduce(*(g->myGetPolPtr()));
        } else {
          myReductionCog->myMoveToNextLM();
        }
      }
      myReductionCog->myRelease(elem);
      if (IsFractionField(CoeffRing(myPolyRing))) {
        myPolyRing->myDivByCoeff(raw(elem), raw(LC(elem)));
        elem = ClearDenom(elem);
      }
    }

    void JanetTree::myJHeadNormalForm(RingElem& elem) {
      if(IsZero(elem)) {
        return;
      }
      ReductionCog myReductionCog(NewRedCogPolyField(myPolyRing));
      //skipping the leading-monomial (we want to perform a !tail!-reduction)
      myReductionCog->myAssignReset(elem);
      myReductionCog->myMoveToNextLM();
      while (!IsActiveZero(myReductionCog)) {
        JanetTriple* g = myJDivisor(ActiveLPP(myReductionCog));
        if (g != 0) {
          myReductionCog->myReduce(*(g->myGetPolPtr()));
        } else {
         break;
        }
      }
      myReductionCog->myRelease(elem);
      if (IsFractionField(CoeffRing(myPolyRing))) {
        myPolyRing->myDivByCoeff(raw(elem), raw(LC(elem)));
        elem = ClearDenom(elem);
      }
    }

    bool JanetTree::IamInternalNode(ConstRefPPMonoidElem w) const
    {
      long i(0);
      JanetIterator iter(*this);
      long n(NumIndets(myPolyRing));
      const std::vector<long> expv = exponents(w);
      while (i < n) {
        while ((iter.myCurrentDeg() < expv[i]) && (iter.myDisNextDeg())) {
          iter.myNextDeg();
        }
        if (iter.myCurrentDeg() < expv[i])
        {
          return false;
        }
        CoCoA_ASSERT(iter.myCurrentDeg() >= expv[i]);
        if (iter.myCurrentDeg() == expv[i]) {
          if (iter.myDisNextVar() == 0) {
            for (long j = i + 1; j != n; ++j) {
              if (expv[j] != 0) {
                return false;
              }
            }
            return true;
          } else {
            long dis(iter.myNextVar());
    //      check if everything between i and dis is zero
            for (long j = i + 1; j != i + dis; ++j)
            {
              if(expv[j] != 0) {
                return false;
              }
            }
            i += dis;
          }
        } else { // iter.myCurrentDeg() > expv[i
          for (long j = i + 1; j != n; ++j) {
            if (expv[j] != 0) {
              return false;
            }
          }
          return true;
        }
      }
      return true;
    }

    std::vector<nodeData> JanetTree::myInsertInTreeStructure(const nodeData& triple) {
      CoCoA_ASSERT(IamInternalNode(LPP(triple->myGetPol())));
      long i(0);
      JanetIterator iter(*this);
      long n(NumIndets(myPolyRing));
      long maxPosDegree(-1);
      std::vector<nodeData> result;
      const std::vector<long> expv = exponents(LPP(*(triple->myGetPolPtr())));
      for (long j = n - 1; j != -1; --j) {
        if (expv[j] != 0) {
          maxPosDegree = j;
          break;
        }
      }
      while (i < n) {
        while ((iter.myCurrentDeg() < expv[i]) && (iter.myDisNextDeg())) {
          iter.myNextDeg();
        }
        // iter.myCurrentDeg() < exponent of expv[i] is not possible because of assertion above
        CoCoA_ASSERT(expv[i] <= iter.myCurrentDeg());
        if (expv[i] == iter.myCurrentDeg()) {
          if (i == maxPosDegree) {
    //      Every expv[j] = 0 for j in {i + 1, n - 1}
            if (iter.IamLeafNode()) {
              break;
            }
            result = iter.myAllNodesAboveMeIncludingMe();
            iter.myChangeToLeafNode(triple);
            break;
          }
    //      there is definitely a next var node because of assertion above
    // TODO     maybe it is exactly a leaf node???? (both triple or i)
          i += iter.myNextVar();
        } else { //expv[i] < iter.myCurrentDeg()
          CoCoA_ASSERT(expv[i] < iter.myCurrentDeg());
          result = iter.myAllNodesAboveMeIncludingMe();
          iter.myPrevDeg();
          CoCoA_ASSERT(expv[i] > iter.myCurrentDeg());
          iter.mySetNextDeg(expv[i] - iter.myCurrentDeg());
          iter.myNextDeg();
          iter.myChangeToLeafNode(triple);
          CoCoA_ASSERT(iter.myCurrentDeg() == expv[i]);
          break;
        }
      }
      return result;
    }

    std::vector<nodeData> JanetIterator::myAllNodesAboveMeIncludingMe() const {
      std::vector<nodeData> result;
      if (IamLeafNode()) {
        JanetLeafNodeImpl* node(dynamic_cast<JanetLeafNodeImpl*>(myCurIter->get()));
        CoCoA_ASSERT(node != 0);
        result.push_back(node->myGetTripleListPointer());
      } else {
        if (myDisNextVar() > 0) {
          JanetIterator varIter(*this);
          varIter.myNextVar();
          std::vector<nodeData> varRes(varIter.myAllNodesAboveMeIncludingMe());
          result.insert(result.end(), varRes.begin(), varRes.end());
        }

        if (myDisNextDeg() > 0) {
          JanetIterator degIter(*this);
          degIter.myNextDeg();
          std::vector<nodeData> degRes(degIter.myAllNodesAboveMeIncludingMe());
          result.insert(result.end(), degRes.begin(), degRes.end());
        }
      }
      return result;
    }

    void JanetContainer::myChangeBasis(std::list<JanetTriple>::const_iterator newListBegin,
                                       std::list<JanetTriple>::const_iterator newListEnd)
    {
      CoCoA_ASSERT(newListBegin != newListEnd);
      CoCoA_ASSERT(owner(newListBegin->myGetPol()) == myTree.myGetPolyRing());
      myList.clear();
      myTree.myDelete();

      myList.insert(myList.end(), newListBegin, newListEnd);
      myInitializeTreeFromList();
    }

    void JanetContainer::myInitializeTreeFromList()
    {
      if (!IsOne(myList.front().myGetPol()))
      {
        for (std::list<JanetTriple>::iterator iter(myList.begin()); iter != myList.end(); ++iter)
        {
          myTree.myInsert(iter);
        }
      }
    }
  } // end of namespace Involutive
} // end of namespace CoCoA
