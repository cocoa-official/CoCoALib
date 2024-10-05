#ifndef CoCoA_TmpMorseGraph_H
#define CoCoA_TmpMorseGraph_H

//   Copyright (c) 2015  John Abbott,  Anna M. Bigatti
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

#include "CoCoA/TmpMorsePaths.H"

namespace CoCoA
{

  class matrix; // fwd decl -- defined in matrix.H

  namespace Involutive
  {
    /*
     * This is the basic class to deal with a MorseGraph. A MorseGraph consists
     * of points (so called morse elements (represented by a MorseElement)) and direct edges (called morse path) between them
     * It is natural that an edge is accessible by its starting point. But in our implementation
     * it turns out, that this would be a bad idea. In our implementation an edge is
     * accessible by its target. The reason for that is the reduction of our graph.
     * One can get nearer information about that during reading the documentation of the relevant methods.
     * All morse path which directed to one edge are collected in the class MorsePaths. We store the MorsePaths as a
     * value of a map, where the corresponding edge is the key.
     *
     * This class should not be used by the normal user. The normal user should use the
     * classes MorseBetti or MorseResolution which compute a graded betti diagram or a (minimal) free resolution
     * for an ideal.
     *
     * The ideal must be represented by a JBMill or PBMill.
     */
    class MorseGraph
    {
    public:
      /*
       * get myResolution as map<MorseElement, MorsePaths>
       */
      /**
       * Returns the internal representation of a resolution
       *
       * @return a map where the keys are the MorseElement and the values are the corresponding
       * MorsePaths
       */
      inline const std::map<MorseElement, MorsePaths>& myGetBasicResolution() const
      { return myResolution; }

      /**
       * prints the elements of the current morse Graph
       */
      std::string toStr() const;
      typedef std::map<std::pair<long, long>, long> BetaVector;
    protected:
      /*
       * Constructor: it uses the JBMill to set several members of MorseGraph
       * in addition to that it stores the janet basis as a vector of JBElems
       *
       * @param mill
       */
      MorseGraph(JBMill mill):
          myMill(mill),
          myRing(mill.myGetSparsePolyRing()),
          myMapRing(myRing)
      {
          myBasis = myTransform(mill);
          mySetBetaVector();
      }

      /**
       * Get begin and end vector of myBasis
       *
       * @return a pair of JBElemConstIters, which represents the begin and end of myBasis
       */
      inline std::pair<MorseElement::JBElemConstIter, MorseElement::JBElemConstIter> myGetBasisRange() const
      { return std::make_pair(myBasis.begin(), myBasis.end()); }

      /**
       * computes the general basis of the form: WedgeProduct |x| (1) (BasisElement)
       * It takes BasisElement b. Then it computes all ordered subsets of the corresponding
       * nonmultiplicative variables. We combine these subsets with the BasisElement to
       * WedgeProduct (the subset) |x| (1) (b)
       * It returns a vector of all generated elements of this form.
       * (For usage of this function take a look at the class documentation of MorseBetti
       * and MorseResolution)
       *
       * @return a vector of MorseElements
       */
      std::vector<MorseElement> myComputeGeneralBasis() const;



      /**
       * delete element at position iter of myResolution and returns iterator to
       * formerly previous element
       *
       * @param iter
       *
       * @return iterator to formerly previous element of iter
       */
      std::map<MorseElement, MorsePaths>::reverse_iterator myDeleteAndJumpToPrevious(std::map<MorseElement, MorsePaths>::reverse_iterator iter);

      /**
       * transform vector<pair<RingElem, vector<bool> > > to vector<pair<RingElem, DynamicBitset> >
       * This function transforms the JBMill to an list of JBElems. A MorseElement stores a reference to this list.
       *
       * @param mill const ref to JBMill
       *
       * @return vector of JBElems
       */
      std::vector<MorseElement::JBElem> myTransform(const JBMill& mill) const;

      /**
       * Morse Reduction the right part:
       * This function reduces our MorseGraph with so called reductions of type 2 or simply right reductions
       * It first checks if maximum and LongIter are different and LongIter is not an multiplicative
       * variable of the JBElem which is contained in origin.
       * If this is not the case there is nothing to do.
       * In the other case we compute the inv. standard representation (isr) of
       * (x_LongIter * m.rf * m.b)/x_maximum using the container (m.rf is the right product of origin and
       * m.b the janet basis element of origin).
       * The isr is of the form \sum(\sum(c_{ij} t_{ij}) b_i).
       * We iteratore over all c_{ij} t{ij} b_i to some prechecks (the author is too lazy to describe them)
       * and then generate a new MorseElement if the prechecks show that this could be a valid MorseElement
       * which not reduces directly to zero. After that it adepts the maps (represented by paths)
       *
       * @param paths
       * @param NewWedgeProduct
       * @param origin
       * @param maximum
       * @param LongIter
       * @param container StandardRepresentationConatainer to speed up the computation and avoid redundant normal form computations
       */
      void myRightMinimization(const PathMap& paths,
                               DynamicBitset NewWedgeProduct,
                               const MorseElement& origin,
                               long maximum,
                               long LongIter,
                               StandardRepresentationContainer& container);

      /**
       * add maps from origin to maps[i].first to the MorseGraph
       *
       * @param maps
       * @param origin
       */
      void myAddMapsToResolution(const std::vector<std::pair<MorseElement, RingElem> >& maps, const ConstResIter& origin);

      /**
       * This function performs a lexicographic comparison of the janet basis elems.
       * It compares the LPPs of e1 and e2 lexicographically
       *
       * @param e1
       * @param e2
       * @return returns true if e1 <_lex e2
       */
      static bool myLexCompareJBElems(std::vector< std::pair<RingElem, std::vector<bool> > >::const_iterator e1,
                                      std::vector< std::pair<RingElem, std::vector<bool> > >::const_iterator e2);

      /**
       * Recursive function. It generates all orderd subset of InputSet.
       * It is a call-by-ref function. The result is stored in result.
       *
       * @param result
       * @param CurrentResult
       * @param InputSet
       * @param length
       */
      void myVariationWithoutRepetition(std::vector<std::vector<long> >& result,
                                        const std::vector<long>& CurrentResult,
                                        const std::vector<long>& InputSet,
                                        long length)  const;

      /*
       * the basic datastructure
       * myRing and myPPM are basic datas about the ring itself
       * myMill represents the janet basis, but for faster access we also save the JB also
       * in myBasis which is a vector of JBElems (for more informations take a look at the class documentation of JBElems)
       * myResolution contains the resolution as a MorseGraph i.e. a list of nodes and directed edges between
       * this notes.
       * The edges have a value. For computing the betti numbers it is suitable to compute these maps
       * over the base field. For computing the free resolution we have to compute over the original polynomial ring (myRing)
       * myMapRing represents these two possibilities.
       */
      JBMill myMill;
      SparsePolyRing myRing;
      std::map<MorseElement, MorsePaths> myResolution;
      std::vector<MorseElement::JBElem> myBasis;
      ring myMapRing;
      BetaVector myBetaVector;
    private:
      void mySetBetaVector();
      /**
       * compute all possible wedge products of NonMultVars.
       *
       * @param NonMultVars
       * @param length
       *
       * @return a vector of DynamicBitsets which represents all possible wedges
       */
      std::vector<DynamicBitset> myPossibleWedges(const std::vector<long>& NonMultVars) const;

      /**
       * Creates a new basis element (in the next step we compute the standard rep of this basis elemenbt)
       * (x_IndexMult* m.rf * m.b)/x_IndexDiv
       *
       * @param m
       * @param IndexMult
       * @param IndexDiv
       * @return returns the new element as RingElem
       */
      RingElem myCreateNewBasisElement(const MorseElement& m, long IndexMult, long IndexDiv) const;
    };

    /**
     * MorsePaths to stream
     *
     * @param os
     * @param obj
     *
     * @return stream object
     */
    inline std::ostream& operator<<(std::ostream& os, const MorseGraph& obj)
    {
      os << obj.toStr();
      return os;
    }

    /**
     * Computes the graded Betti numbers of the ideal represented by mill.
     *
     * @param mill
     * @return matrix of BigInts which represents the graded Betti numbers
     */
    matrix BettiDiagram(JBMill mill);

    /**
     * Computes the graded Betti numbers of the ideal represented by mill which
     * are in the col-th column of the BettiDiagram
     *
     * @param mill
     * @param col
     *
     * @return matrix of RingZZ which consists of only one row
     */
    matrix BettiColumn(JBMill mill, long col);

    /**
     * Computes the graded Betti numbers of the ideal represented by mill which
     * are in the col-th column and in row greater or equal MinRow of the BettiDiagram
     *
     * @param mill
     * @param MinRow
     * @param col
     *
     * @return matrix of RingZZ which consists of only one row. The first entry in the matrix
     * is the BettiNumber at position (MinRow, col) in the associated BettiDiagram
     */
     matrix BettiPartialColumn(JBMill mill, long MinRow, long col);

    /**
     * Computes !!only!! the Betti number which is at position (row, col) in the BettiDiagram
     * of the ideal which is associated to mill
     *
     *
     * @param mill
     * @param row
     * @param col
     *
     * @return a RingElem of type RingZZ
     */
    RingElem BettiNumber(JBMill mill, long row, long col);

    /**
     * Computes the graded Betti numbers of the ideal represented by mill.
     *
     * @param mill
     * @return a pair of matrices of BigInts where the first entry represents the graded
     * Betti numbers and the second entry represents the graded pseudo-Betti numbers
     */
    std::pair<matrix, matrix> PseudoBettiDiagram(JBMill mill);

    /**
     * Computes the minimal free resolution of the ideal represented by mill.
     *
     * @param mill
     * @return a vector of matrices which represents the minimal free resolution.
     */
    std::vector<matrix> MinimalResolution(JBMill mill);

    /**
     * Computes a free resolution of the ideal represented by mill.
     *
     * @param mill
     * @return a vector of matrices which represents a free resolution.
     */
    std::vector<matrix> Resolution(JBMill mill);
  } // enf of namespace Involutive
} // end of namespace CoCoa
#endif

