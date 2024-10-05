#ifndef CoCoA_TmpMorseBetti_H
#define CoCoA_TmpMorseBetti_H

//   Copyright (c)  2015 John Abbott,  Anna M. Bigatti
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

#include "CoCoA/TmpMorseGraph.H"

namespace CoCoA
{
  namespace Involutive {

    class MorseBetti : public MorseGraph {
    public:
      MorseBetti(JBMill mill): MorseGraph(mill)
      {
        myMapRing = CoeffRing(myRing);
      }

      /**
       * computes the graded Betti numbers
       *
       * @return a matrix of BigInts which represents the Betti diagram
       */
      matrix myComputeBettiNumbers();

      /*
       * computes the pseudo Betti numbers and the 'normal' Betti numbers
       *
       * @return a pair of matrices of BigInts where the first entry represents the graded
       * Betti numbers and the second entry represents the graded pseudo-Betti numbers
       */
      std::pair<matrix, matrix> myComputePseudoBettiNumbers();
    protected:
      /**
       * This methods takes all elementary MorseElements (MorseElements of the form WedgeProduct \tensor (1) (BasisElement))
       * and compute all constant maps.
       *
       * @param elems
       * @param container StandardRepresentationConatainer to speed up the computation and avoid redundant normal form computations
       */
      void myComputeBasicConstantGraph(const std::vector<MorseElement>& elems, StandardRepresentationContainer& container);

      /**
       * Compute the constant Morse reduction of our resolution. It applies to every element which is not of the form
       * WedgeProduct \tensor (1) (BasisElement) the right minimization.
       *
       * @param container StandardRepresentationConatainer to speed up the computation and avoid redundant normal form computations
       */
      void myConstantDirectMorseReduction(StandardRepresentationContainer& container);

      /**
       * The list of MorseElements are ordered in a strict way. Now we extract all element with NumberWedges or
       * NumberWedges+1 Wedge Basis elements and divide them by degree
       * (degree = Number Wedge elements + degree of basis element) and number of wedges elements. Then we
       * determine the position of a specific MorseElement within all MorseElements with same degree and number
       * of wedge elements and save the position in the value of the map.
       * TODO: Do we really need a pair? Maybe only a MorseElement as key could be enough.
       *
       * @param NumberWedges
       *
       * @return a map with degree and MorseElement as key and position as described above as value.
       */
      std::map<std::pair<long, MorseElement>, long> myGradedPositionMorseElements(long NumberWedges) const;

      /**
       * This method extracts a submatrix out of F_PosInRes where 'image' and 'domain' have same degree
       * First we generate an empty matrix with rows rows and cols columns.
       * Then we iterate over all elements with PosInRes number of wedges and degree degree.
       * During the iteration we iterate over all elements which have a path to the current element and
       * have the same degree. The values of these paths define the submatrix.
       *
       * @param ResIter
       * @param rows
       * @param cols
       * @param degree
       * @param PosInRes
       * @param identifier
       *
       * @return a matrix of RingElems
       */
      matrix ConstructDegreeMatrix(std::map<MorseElement, MorsePaths>::const_iterator& ResIter,
                                   long rows,
                                   long cols,
                                   long degree,
                                   long PosInRes,
                                   const std::map<std::pair<long, MorseElement>, long>& identifier) const;
    private:

      /**
       * Computes the index of the last nonzero entry. If the first entry
       * is the only nonzero entry it returns 1. If all entries are zero it
       * returns 0.
       * @param vec
       *
       * @return the result as long
       */
      long myLastNonZeroIndex(const std::vector<long> vec) const;

      /**
       * Computes the wasteless ranks: This method computes the difference between
       * the Betti diagram and the pseudo Betti diagram using the Morse Graph.
       * For implementation details see inline documentation.
       *
       * @param ranks the pseudo Betti diagram
       * @return pseudo Betti diagram - Betti diagram
       */
      std::vector< std::vector<long> > myComputeWasteRanks(const std::vector< std::vector<long> >& ranks) const;

      /**
       * Computes the wasteless rank for the map from (PosInRes + 1)=F1 to PosInRes=F0.
       * In terms of Betti and pseudo Betti diagramms: We compute the difference
       * for the F1-th column.
       *
       * The map from F1 to F0 can be represented as matrix. These matrix can be divided into submatrices, s.t.
       * the 'domains' of these submatrices have the same degree per submatrice and the 'image' of these submatrices
       * have the same degree. If the degree of the 'domains' and 'images' are the same there can be some non-zero
       * constant entries in the submatrices. We compute the rank of these matrices which give us the difference
       * for the specific degree in the F1-th column between the Betti diagram and the pseudo Betti diagram.
       *
       *
       * @param WasteRanks at the end this will be the result. (Call by Ref)
       * @param ResIter
       * @param ranks
       * @param PosInRes
       */
      void myComputeWasteRanksPerMap(std::vector< std::vector<long> >& WasteRanks,
                                     std::map<MorseElement, MorsePaths>::const_iterator& ResIter,
                                     const std::vector< std::vector<long> >& ranks,
                                     long PosInRes) const;

      /**
       * Compute the ranks of the free resolution. It iterates over all
       * MorseElements and determine the rank of these elements and counts how
       * many elements of same degree with same number of wedge products are present.
       * @return a matrix represented by a vector of vector of longs.
       */
      std::vector< std::vector<long> > myComputeRanks() const;

      /**
       * Computes the number of rows in the Betti diagram.
       * This is length of rows of ranks minus length of cols of ranks
       *
       * @param ranks
       *
       * @return number of rows in Betti diagram
       */
      long myNumRowsBettiDiagram(const std::vector<std::vector<long> >& ranks) const;

      /**
       * Compute the number of cols in the Betti diagram.
       * This is the length of ranks minus number of zero cols at end of ranks
       *
       * @param ranks
       * @return number of cols in Betti diagarm
       */
      long myNumColsBettiDiagram(const std::vector<std::vector<long> >& ranks) const;

      /**
       * transform the rank matrix to a Betti matrix (so called 'Macaulay style').
       * This is only some crappy index shifting.
       *
       * @param ranks
       *
       * @return matrix of RingZZs
       */
      matrix myTransformRanksToBettis(const std::vector<std::vector<long> >& ranks) const;

      /**
       * m1 - m2 for these two pseudo matrices
       *
       * @param m1
       * @param m2
       * @return m1 - m2
       */
      std::vector<std::vector<long> > myMatrixMinus(std::vector<std::vector<long> > m1, const std::vector<std::vector<long> >& m2) const;

      /**
       * Computes wasteless ranks per degree for map from PosInRes + 1 to PosInRes position.
       * First it extracts the corresponding submatrices in F_PosInRes
       * then it computes the ranks of these submatrices
       *
       * @param ColRanks
       * @param PosInRes
       *
       * @return a vector of pairs. In a pair the first entry identifies the degree and the second entry
       * identifies the wasteless rank for this degree.
       */
      std::vector<std::pair<long, long> > myComputeWasteRanksPerDegree(std::map<MorseElement, MorsePaths>::const_iterator& ResIter, const std::vector<long>& RowRanks, const std::vector<long>& ColRanks, long PosInRes) const;


      /**
       * It computes the 'constant resolution': First it computes the basic constant graph
       * and then it applies the constant morse reduction.
       */
      void myComputeConstantResolution();

    };

  } // end of namespace Involutive
} // end of namespace CoCoa
#endif

