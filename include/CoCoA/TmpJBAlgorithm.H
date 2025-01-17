#ifndef CoCoA_TmpJBAlgorithm_H_
#define CoCoA_TmpJBAlgorithm_H_

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

#include "TmpJBSets.H"

namespace CoCoA {
  namespace Involutive {

    /****************************************************************************************************************************/
    /*!\brief This class compares the LPP of the polynomials in Janet Triples
     */
    /*****************************************************************************************************************************/
    struct IsLowerLPP {
      /**
       * The "compare" method
       * @param t1 JanetTriple
       * @param t2 JanetTriple
       * The output is: LPP(t1.myPolynomial)<LPP(t2.myPolynomial)
       **/
      bool operator()(const JanetTriple& t1, const JanetTriple& t2) {
        return ((LPP(t1.myPolynomial) < LPP(t2.myPolynomial))
            || ((LPP(t1.myPolynomial) == LPP(t2.myPolynomial)) && (t1.myAncestor < t2.myAncestor)));
      }
    };

    /**
     * @brief Virtual algorithm class
     * @details Every algorithm for computing Janet basis should be a subclass of
     * JBAlgorithm. It provides interfaces for computing and returning Janet basis.
     */
    class JBAlgorithm {
    public:
      /*
        virtual constructor because of inheritance
      */
      virtual ~JBAlgorithm() {
      }

      /**
       * Purely virtual method to compute a Janet basis. Must be specified in subclasses.
       * @param beginInput a const iterator to a vector of RingElems. Should be the beginning of the vector.
       * @param endInput a const iterator to a vector of RingElems. Should be the end of the vector.
       */
      virtual void myComputer(std::vector<RingElem>::const_iterator beginInput,
                              std::vector<RingElem>::const_iterator endInput) = 0;

      /**
       * This purley virtual method returns the result. Must be specified in subclasses.
       * @return It returns a JanetContainer containing the computed basis.
       */
      virtual JanetContainer myOutputResult() = 0;

    protected:
      /**
       * Proteced Constructor. The user should not use this class directly.
       * The constructor initialize the class members myPolyRing and myPPMValue.
       * @param polyRing
       */
      JBAlgorithm(SparsePolyRing polyRing)
          : myPolyRing(polyRing),
            myPPMValue(PPM(myPolyRing)) {
      }

      const SparsePolyRing myPolyRing;

      const PPMonoid myPPMValue;
    };

    class TQAlgorithm : public JBAlgorithm {
    public:
      virtual ~TQAlgorithm() {
      }

      /**
       * See at base class
       */
      virtual void myComputer(std::vector<RingElem>::const_iterator beginInput,
                              std::vector<RingElem>::const_iterator endInput) = 0;

      /**
       * Returns the result.
       * It extracts the computed Janet Basis out of a TQSet as a list and inserts
       * the result in a JanetContainer.
       * @return JanetConatiner
       */
      JanetContainer myOutputResult();

    protected:
      /**
       * Constructor: It calls the constructor of the base class and initializes
       * a myJTree.
       * @param polyRing
       */
      TQAlgorithm(SparsePolyRing polyRing)
          : JBAlgorithm(polyRing),
            myJTree(myPolyRing, 0, 0) {
      }

      /**
       * A purely virtual method. Must be defined in a subclass.
       * It access TQSets (or a subclass of it) which is defined in subclasses
       * @return A reference to TQSets
       */
      virtual TQSets& myGetSets() = 0;

      /**
       * The method initializes a TQSet which was accessed by myGetSets(). It inserts
       * the elements of a vector<RingElem> between beginInput and endInput
       * @param beginInput
       * @param endInput
       */
      void myInitialization(std::vector<RingElem>::const_iterator beginInput,
                            std::vector<RingElem>::const_iterator endInput);

      /**
       * It modifies the instance such that it represents the trivial ideal:
       * The set T in TQSets contains only the JanetTriple of the RingElem 1.
       */
      void myTrivialIdeal();

      JanetTree myJTree;
    };

    /**
     * The class is a specialization of the TQAlgorithm. It uses the Degree TQ Algorithm
     * to compute a Janet Basis. This algorithm is described in
     * "Gerdt - Construction of Janet Bases II. Polynomial Bases"
     */
    class DegreeTQ : public TQAlgorithm {
    public:
      /**
       * Constructor: It calls the constructor of the base class (TQAlgorithm) and initializes
       * mySets with the option crit and the ring polyRing
       *
       * @param polyRing
       * @param crit a bitflag with three elements. If a bit is true mySets uses the corresponding involutive
       * criteria
       */
      DegreeTQ(SparsePolyRing polyRing, std::bitset<3> crit = std::bitset<3>(0))
          : TQAlgorithm(polyRing),
            mySets(myPolyRing, crit) {
      }

      /**
       * accessor for the base class to mySets
       * @return TQSets
       */
      TQSets& myGetSets() {
        return mySets;
      }

      /**
       * Implementation of the algorithm to compute the Janet Basis via the DegreeTQ strategy
       * @param beginInput
       * @param endInput
       */
      void myComputer(std::vector<RingElem>::const_iterator beginInput, std::vector<RingElem>::const_iterator endInput);

    private:
      TQSets mySets;
    };

    /**
     * The class is a specialization of the TQAlgorithm. It uses the Degree TQ Algorithm
     * to compute a Janet Basis. This algorithm is described in
     * "Gerdt - "On Computing Janet Bases for Degree Compatible Orderings"
     */
    class BlockTQ : public TQAlgorithm {
    public:

      /**
       * Constructor: It calls the constructor of the base class (TQAlgorithm) and initializes
       * mySets with the options crit, strategy and the ring polyRing
       *
       * @param polyRing
       * @param crit a bitflag with three elements. If a bit is true mySets uses the corresponding involutive
       * criteria
       * @param strategy specified the concrete strategy. One could choose between TQPSets::Low (default) or TQPSets::High
       */
      BlockTQ(SparsePolyRing polyRing, std::bitset<3> crit = std::bitset<3>(0), TQPSets::Strategy strategy =
                    TQPSets::Low)
          : TQAlgorithm(polyRing),
            mySets(myPolyRing, crit, strategy) {
      }

      /**
       * accessor for the base class to mySets
       * @return TQSets
       */
      TQSets& myGetSets() {
        return mySets;
      }

      /**
       * Implementation of the algorithm to compute the Janet Basis via the DegreeTQ strategy
       * @param beginInput
       * @param endInput
       */
      void myComputer(std::vector<RingElem>::const_iterator beginInput, std::vector<RingElem>::const_iterator endInput);

    private:
      TQPSets mySets;
    };

    /**
     * The class is a specialization of JBAlgorithm. It uses a polynomial completion to compute a Janet basis.
     * It first uses the method ComputeGBasis of CoCoALib to compute a Groebner basis for the input.
     * Afterwards it complete it to a Janet basis.
     */
    class CompletionGB : public JBAlgorithm {
    public:
      /**
       * Constructor: It calls the constructor of the base class (JBAlgorithm) and initializes myJTree
       *
       * @param polyRing
       */
      CompletionGB(SparsePolyRing polyRing)
          : JBAlgorithm(polyRing),
            myJTree(myPolyRing, 0, 0) {
      }

      /**
       * Implementation of the algorithm to compute the Janet Basis via the completion algorithm.
       * @param beginInput
       * @param endInput
       */
      void myComputer(std::vector<RingElem>::const_iterator beginInput, std::vector<RingElem>::const_iterator endInput);

      /**
       * Returns the result.
       * It initializes a JanetContainer with myPolyRing and myJBasis
       * @return JanetConatiner
       */
      JanetContainer myOutputResult();

    private:
      typedef std::pair<RingElem, PPMonoidElem> elemPair;
      typedef std::list<elemPair> elemPairList;
      typedef std::list<elemPair>::iterator elemPairIter;

      /**
       * This class is used to compare two pairs p1, p2 of <RingElem, PPMonoidElem>
       * p1 < p2 if LPP(p1.RingElem) < LPP(p2.RingElem)
       */
      struct LPPCompare {
        bool operator()(const elemPair& r1, const elemPair& r2) {
          return LPP(r1.first) < LPP(r2.first);
        }
      } myLPPComparison;

      /**
       * This class is used to compare two JanetTriples p1, p2
       * p1 < p2 if LPP(p1.pol) < LPP(p2.pol)
       */
      struct TripleLPPCompare {
        bool operator()(const JanetTriple& r1, const JanetTriple& r2) {
          return LPP(r1.myGetPol()) < LPP(r2.myGetPol());
        }
      } myLPPTripleComparison;

      /**
       * This method first computes the reduced Groebner Basis out out the Input beginInput .. endInput.
       * Then it inserts the result in the Janet Tree. Due to the fact, that the reduced Groebner Basis
       * is contained in the Janet basis we do not get any collisions.
       *
       * @param beginInput
       * @param endInput
       */
      void myBuildInitialTree(std::vector<RingElem>::const_iterator beginInput,
                              std::vector<RingElem>::const_iterator endInput);

      /**
       * This method checks if the PPMonoidElem m is multiplicative for index
       * To decide this it moves through the JanetTree
       *
       * @param m [description]
       * @param index [description]
       *
       * @return true is m is nonmultiplicative for index
       */
      bool myElemIsNonMultForIndex(const PPMonoidElem& m, long index);

      /**
       * It insert the pair (elem, anc) to list (at the end)
       *
       * @param list it is a reference to a list
       * @param elem
       * @param anc
       */
      void myPushToProlongList(elemPairList& list, const RingElem& elem, const PPMonoidElem& anc);

      /**
       * It insert the pair (elem, anc) to list (s.t. the list remains sorted after LPP(elem))
       *
       * @param list it is a reference to a list
       * @param elem
       * @param anc
       */
      void myInsertToProlongList(elemPairList& list, const RingElem& elem, const PPMonoidElem& anc);

      /**
       * This method takes the Janet basis which is represented by myJBasis and myJTree which is possibly
       * not minimal and sorts it by degree and minimizes it.
       */
      void mySortAndMinimizeJanetBasis();

      std::list<JanetTriple> myJBasis;

      JanetTree myJTree;

    };
  }
} /* namespace CoCoA */

#endif /* TMPJBALGORITHM_H_ */
