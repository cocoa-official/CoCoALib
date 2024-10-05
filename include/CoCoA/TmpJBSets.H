#ifndef CoCoA_JBSets_H
#define CoCoA_JBSets_H

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
#include <set>
#include <bitset>
#include "CoCoA/ReductionCog.H"

namespace CoCoA {
  namespace Involutive {
    //typedef for easier writing
    typedef std::list<JanetTriple>::iterator Iterator;

    /*
     * TQSets stores the sets T and Q of our TQ-Degree algorithm. In our first version we store the
     * elements in two multisets. It turns out that the movement of elements from T to Q is very costly.
     * Therefore we create the class TQSets. This class has a basic Set which store all elements form T and Q.
     * The set T and Q contain only iterators of the basic Set. Now the movement of elements is very cheap.
     * In addition, the class contains some efficient reduction procedures (in comparison to the previous versions)
     */
    class TQSets {
    public:
      /*
       * definition of the inner class
       * for comparing elements in the set T and Q
       */
      class CompareIterator {
       public:
        /**
         * This function compares two JanetTriples: i1 < i2 iff (i1.pol < i2.pol) || (i1.pol == i2.pol) && (i1.anc < i2.anc)
         * @param  i1 is the Iterator which points to the first JanetTriple
         * @param  i2 is the Iterator which points to the second JanetTriple
         * @return returns true if i1 < i2 otherwise false
         */
        bool operator()(const Iterator& i1, const Iterator& i2) const {
          PPMonoid myPPM(PPM(owner(i1->myPolynomial)));
          //  return (LPP(i1->myPolynom) < LPP(i2->myPolynom)) || ((LPP(i1->myPolynom) == LPP(i2->myPolynom)) && ((i1->myAncestor < i2->myAncestor) || ((i1->myAncestor == i2->myAncestor) && (((AsPolyRing(owner(i1->myPolynom))->myNumTerms(raw(i1->myPolynom))) < AsPolyRing(owner(i1->myPolynom))->myNumTerms(raw(i2->myPolynom))) || ((AsPolyRing(owner(i1->myPolynom))->myNumTerms(raw(i1->myPolynom)) == AsPolyRing(owner(i1->myPolynom))->myNumTerms(raw(i2->myPolynom))) && (abs(LC(i1->myPolynom)) < abs(LC(i2->myPolynom))))))));
          //return (myRing->myLPP(raw(i1->myPolynom)) < myRing->myLPP(raw(i2->myPolynom))) ||((myRing->myLPP(raw(i1->myPolynom)) == myRing->myLPP(raw(i2->myPolynom))) && ((i1->myAncestor < i2->myAncestor) ));
          // int comp(myPPM->myCmp(raw(myRing->myLPP(raw(i1->myPolynom))), raw(myRing->myLPP(raw(i2->myPolynom)))));
          const int comp = myPPM->myCmp(raw(i1->myGetLppPol()), raw(i2->myGetLppPol()));
          return (comp < 0) || (comp == 0 && (myPPM->myCmp(raw(i1->myAncestor), raw(i2->myAncestor)) < 0));
          //  return (LPP(i1->myPolynom) < LPP(i2->myPolynom)) || ((LPP(i1->myPolynom) == LPP(i2->myPolynom)) && ((i1->myAncestor < i2->myAncestor) || ((i1->myAncestor == i2->myAncestor) && (((abs(LC(i1->myPolynom)) < abs(LC(i2->myPolynom))) || ((abs(LC(i1->myPolynom)) < abs(LC(i2->myPolynom))) &&((AsPolyRing(owner(i1->myPolynom))->myNumTerms(raw(i1->myPolynom))) < AsPolyRing(owner(i1->myPolynom))->myNumTerms(raw(i2->myPolynom)))))))));
        }
      };

      /**
       * typedef
       */
      typedef std::multiset<Iterator, CompareIterator>::iterator MultisetIterator;
      typedef std::multiset<Iterator, CompareIterator>::reverse_iterator MultisetReverseIterator;

      /*
       * constructor
       */
      TQSets(const SparsePolyRing& ring, std::bitset<3> criteria = std::bitset<3>(0))
          : myPolyRing(ring),
            myCriteria(criteria),
            myReductionCog(NewRedCogPolyField(myPolyRing)) {
      }

      virtual ~TQSets() {
      }

      /**
       * tests if the set mySetQ is empty
       *
       * @return true if mySetQ is empty
       */
      inline bool IamEmptySetQ() const {
        return mySetQ.empty();
      }

      /**
       * myHeadRecueSetQ performs an involutive headreduction on every element in mySetQ with respect to the JanetTree
       *
       * @param  JTree is the JanetTree
       *
       * @return It returns the size of mySetQ after the headreduction
       */
      long myHeadReduceSetQ(JanetTree& JTree);  //done

      /**
       * myMinimizeAndInsertSetT select in the first step the smallest element in the multiset mySetQ
       * (with respect to the number of terms and (if possible) the absolutue value of the leading coefficient)
       * Then the algorithm reduces the tail of the element and checks if the leading term of this element divides the leading
       * term of any element in mySetT. If this is the case it moves all elements which are greater then the selected element to mySetQ.
       * At the end it inserts the selected element in mySetT and returns the Iterator to this element.
       *
       * @param  JTree is the JanetTree
       *
       * @return       see above
       */
      Iterator myMinimizeAndInsertSetT(JanetTree& JTree);

      /**
       * myTailReduceSetT performs a involutive tailreduction on every element in mySetT which are bigger then TripleIter
       *
       * @param JTree      is the Janeet Tree
       * @param TripleIter an interator to a JanetTriple
       */
      void myTailReduceSetT(JanetTree& JTree, const std::list<JanetTriple>::iterator& TripleIter);

      /**
       * Returns an iterator to the beginning of mySetT
       *
       * @return a MultisetIterator
       */
      MultisetIterator myBeginSetT() const;  //done

      /**
       * Returns an iterator to the end of mySetT
       *
       * @return a MultisetIterator
       */
      MultisetIterator myEndSetT() const;  //done

      /**
       * Returns an iterator to the beginning of mySetQ
       *
       * @return a MultisetIterator
       */
      MultisetIterator myBeginSetQ() const;

      /**
       * Returns an iterator to the end of mySetQ
       *
       * @return a MultisetIterator
       */
      MultisetIterator myEndSetQ() const;

      /**
       * Virtual: Clear mySetT, mySetQ and myBasicSet
       */
      virtual void myClearSets();

      /**
       * myInsertSetQ insert triple in myBasicSet and the corresponding Iterator in mySetQ
       *
       * @param triple is a JanetTriple
       */
      void myInsertSetQ(JanetTriple& triple);

      /**
       * myInsertSetQ insert triple in myBasicSet and the corresponding Iterator in mySetQ
       *
       * @param triple is a JanetTriple
       */
      void myInsertSetT(JanetTriple& triple);

      /**
       *
       * @return returns size of mySetT
       */
      long mySizeSetT();

      /**
       * myJTailNormalForm computes the involutive tail normalform of iter w.r.t JTree
       *
       * @param JTree the janet tree
       * @param iter  the iterator which contains a janettriple
       */
      void myJTailNormalForm(JanetTree& JTree, Iterator iter);

      /**
       * insert the triple which is in ps in the JanetTree and prolong.
       *
       * @param t ref to ListIter to JanetTriple
       * @param JTree ref to JanetTree
       */
      void myJInsert(std::list<JanetTriple>::iterator& ps, JanetTree& JTree);

      /**
       * This algorithm prolong all triples which are behind iter with x_index and insert the prolongation into setQ.
       * If t is a JanetTriple behind iter, a prolongation with index is: t.pol() * x_index
       *
       * @param iter a JanetIterator which represents a tree containing everything 'behind'
       * this iterator
       * @param index is the index of the variable to prolong
       */
      void myJProlong(JanetIterator& iter, const long& index);

    protected:
      /*
       * our basic set
       */
      std::list<JanetTriple> myBasicSet;

      /*
       * the set T
       */
      std::multiset<Iterator, CompareIterator> mySetT;

      /*
       * the set Q
       */
      std::multiset<Iterator, CompareIterator> mySetQ;

      /*
       * our poly Ring
       */
      SparsePolyRing myPolyRing;

      /*
       * which involutive criteria the algorithm shall use?
       */
      std::bitset<3> myCriteria;

      /*
       * reduction cog to perform reductions
       */
      ReductionCog myReductionCog;

      /**
       * IAmCriteria checks if p satifies one of the selected criteria (in myCritieria) with respect to g
       *
       * @param  p is the first JanetTriple
       * @param  g is the second JanetTriple
       *
       * @return returns true if p satifies one of the selected criteria otherwiese false
       */
      bool IamCriteria(const JanetTriple& p, const JanetTriple& g) const;

      /**
       * IamCriteria3 checks if p (from IAmCriteria) satifies the third involutive Criteria with respect to g (from IAmCriteria)
       *
       * @param  AncP is the ancestor of p
       * @param  AncG is the ancestor of g
       *
       * @return returns true if AncP satifies the third criteria
       */
      bool IamCriteria3(const PPMonoidElem& AncP, const PPMonoidElem& AncG) const;


    private:
      /**
       * Choose the smallest element
       *
       * @return a MultisetIterator
       */
      MultisetIterator myChooseSmallestElement();
    };

    class TQPSets : public TQSets {
    public:
      enum Strategy {
        High,
        Low
      };

      /**
       * Constructor
       *
       * @param ring
       * @param criteria (default is nothing)
       * @param strategy (default is low)
       */
      TQPSets(const SparsePolyRing& ring, std::bitset<3> criteria = std::bitset<3>(0), Strategy strategy = Low)
          : TQSets(ring, criteria),
            myStrategy(strategy) {
      }

      /**
       * Inserts triple into mySetP
       *
       * @param triple
       */
      void myInsertSetP(JanetTriple& triple);

      /**
       * Virtual: Clear mySetT, mySetQ, mySetP and myBasicSet
       */
      void myClearSets();

      /*
       * Tests if the set mySetP is empty
       *
       * @return returns true if SetP is empty
       */
      inline bool IamEmptySetP() const {
        return mySetP.empty();
      }

      /**
       * myTailReduceSetTAll performs a involutive tailreduction on every element in mySetT
       *
       * @param JTree      is the Janeet Tree
       */
      void myTailReduceSetTAll(JanetTree& JTree);

      /**
       * Returns iterator to the beginning mySetP
       *
       * @return MultisetIterator
       */
      MultisetIterator myBeginSetP() const;

      /**
       * Returns iterator to the end mySetP
       *
       * @return MultisetIterator
       */
      MultisetIterator myEndSetP() const;

      /**
       * Moves elements from mySetQ to mySetQ. It move all elements which have the lowest StdDeg.
       * The function also performs a full involutive reduction at these elements.
       *
       * @param JTree is the JanetTree
       */
      void myMoveFromQtoP(JanetTree& JTree);

      /**
       * It performs a involutive Tailreduction on every element in mySetP
       *
       * @param JTree is the Janet Tree on which the reduction is based
       */
      void myTailReduceSetP(JanetTree& JTree);

      /**
       * Moves all elements from mySetP to mySetT
       */
      void myInsertSetPInSetT();

      /**
       * @brief IamJUpdate updates mySetP.
       * @details In IamUpdate an element f in mySetP (which is J−reduced modulo mySetT when the subalgorithm is invoked in the main algorithm)
       * is chosen at the initialization with the highest or lowest leading term, depends on the selection strategy used.
       * We put f in the set stack
       * After that, the other polynomial p in mySetP, if any, with the same leading monomial as that in f is reduced modulo f.
       * In the case of nonzero reduction (when monic p is different from monic f) the normal form obtained is added to the stack
       * to be involved in the further reductions.
       * Then, the processes of the selection and reduction of elements in mySetP is repeated until mySetP becomes empty.
       * At the end we copy the content of stack back to mySetP
       *
       * @return true if there is still a prolongation in mySetP
       */
      bool IamJUpdate();

      /**
       * Performs the process described for the function IamJUpdate() choosing always the biggest element in mySetP
       * after each reduction process we have to recompute the order in mySetP.
       *
       * @param stack is the set described in IamJUpdate
       */
      void myJUpdateHighPart(std::list<Iterator>& stack);

      /**
       * Performs the process described for the function IamJUpdate() choosing always the smalles element in mySetP
       * after each reduction process we have to recompute the order in mySetP.
       *
       * @param stack is the set described in IamJUpdate
       */
      void myJUpdateLowPart(std::list<Iterator>& stack);

      /**
       * This function refines the search for the smallest/highest element.
       * They also observed the number of terms and (if ring is an ordered domain) the leading coefficient
       *
       * @param ChooseIter (Call-by-ref) the iterator pointing to the choosen element
       */
      void myChoosingJUpdate(TQSets::MultisetIterator& ChooseIter);

      /**
       * Performs a reduction on PIter. Adepting the ancestor and nonmultiplicative variables
       *
       * @param PIter         iterator to be reduced
       * @param ReductionIter reducing iterator
       */
      void myReductionStepJUpdate(MultisetIterator& PIter, std::list<Iterator>::iterator ReductionIter);

      /**
       * The reduction step myReductionStepJUpdate destroies the ordering in mySetP.
       * This function recomputes the ordering.
       */
      void myRecomputeOrderingSetP();

      /**
       * myTailNormalForm computes the tail normal form of ToRed w.r.t. red (which is only one element!!!)
       *
       * @param ToRed the iterator which contains a janettriple which we want to reduce
       * @param red   the iterator which contains a janettriple which reduces the triple above
       */
      void myTailNormalForm(Iterator ToRed, Iterator red);

      /**
       * myMinimizeSetT moves all elements from mySetT to mySetQ, which have a smaller StdDeg as the beginning element of mySetP
       */
      void myMinimizeSetT();

    private:
      /*
       * the set P
       */
      std::multiset<Iterator, CompareIterator> mySetP;

      /*
       * which strategy should use blockTQ?  high strategy: true, low strategy: false
       */
      Strategy myStrategy;

      /**
       * Moves all elements from mySetP to mySetQ
       */
      void myMoveFromPToQ();
    };
  }
}
#endif
