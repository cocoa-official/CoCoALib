#ifndef CoCoA_TmpPartialMorseBetti_H
#define CoCoA_TmpPartialMorseBetti_H

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

#include "CoCoA/TmpMorseBetti.H"
#include "CoCoA/matrix.H"
#include "CoCoA/MatrixOps.H"

namespace CoCoA
{
  namespace Involutive
  {
    class PartialMorseBetti : public MorseBetti
    {
    public:
      PartialMorseBetti(JBMill mill): MorseBetti(mill)
      {}

      /**
       * Computes !!!only!! the column c of the betti diagram
       *
       * @param c
       * @return Returns a matrix with 1 column, which represents the column c of
       * the betti diagram
       */
      matrix myComputeBettiColumn(long c)
      {
        return myComputeDownToBettiNumber(0, c);
      }

      /**
       * Computes the betti numbers at position (row + i, col) in the betti diagram
       * (i >= 0).B
       *
       * @param row
       * @param col
       *
       * @return returns a matrix with 1 column, which represents the column col of the
       * betti diagram. The first entry of the matrix is the betti number (row, col).
       */
      matrix myComputeDownToBettiNumber(long row, long col);

      /**
       * Computes the betti number at position (row, col) in the betti diagram.
       * This method uses the method myComputeDownToBettiNumber to extract the
       * betti number and everything with larger degree and simply extract the
       * specific betti number.
       *
       * @param row
       * @param col
       *
       * @return a RingElem from the ring RingZZ
       */
      long myComputeBettiNumber(long row, long col);

    // private:
      /**
       * Computes the general basis of the form: WedgeProduct |x| (1) (BasisElement)
       * It takes BasisElement b. Then it computes all ordered subsets of the corresponding
       * nonmultiplicative variables such that these orderd subsets have length lengthWedge or
       * (lengthWedge + 1). We combine these subsets with the BasisElement to
       * (WedgeProduct (the subset) |x| (1) (b)) such that the length of the WedgeProduct
       * plus the degree of b is not less than MinDeg. It returns a vector of all generated
       * elements of this form. (For usage of this function take a look at the class
       * documentation of MorseBetti and MorseResolution)
       *
       * @param MinDeg
       * @param lengthWedge
       *
       * @return a vector of MorseElements
       */
      std::vector<MorseElement> myComputeGeneralColumnBasis(long MinDeg, long lengthWedge) const;

      /**
       * computes the general basis of the form: WedgeProduct |x| (1) (BasisElement)
       * It takes BasisElement b. Then it computes all ordered subsets of the corresponding
       * nonmultiplicative variables such that these orderd subsets have length lengthWedge or
       * (lengthWedge + 1). We combine these subsets with the BasisElement to
       * (WedgeProduct (the subset) |x| (1) (b)) such that the length of the WedgeProduct
       * plus the degree of b is equal to degree. It returns a vector of all generated
       * elements of this form. (For usage of this function take a look at the class
       * documentation of MorseBetti and MorseResolution)
       *
       * @param degree
       * @param lengthWedge
       *
       * @return a vector of MorseElements
       */
      std::vector<MorseElement> myComputeGeneralSingleBasis(long degree, long lengthWedge) const;

      /**
       * It computes the relevant part of the 'constant resolution': First it computes the
       * relevant part of the basic constant graph  and then it applies the constant morse
       * reduction.
       *
       * @param row
       * @param col
       */
      void myComputeConstantSingleResolution(long deg, long numWedges);

      /**
       * It compute all possible elements of the form (WedgeProduct |x| (1) (it))
       * such that WedgeProduct has length lengthWedge. It appends the result to basis.
       * The indices of WedgeProduct have to be a subset of the nonmultiplicative
       * variables of it. The result is stored in basis via call-by-ref
       *
       * @param basis
       * @param it
       * @param lengthWedge
       */
      void myComputeBasisForElement(std::vector<MorseElement>& basis,
                                    MorseElement::JBElemConstIter it,
                                    long lengthWedge) const;

      /**
       * Compute the pseudo betti number \beta_{numWedges, numWedges+deg}.
       * For this it uses the formula to compute the pseudo betti numbers out of the
       * induced free resolution (e.g. see Book Involution)
       *
       * @param deg
       * @param numWedges
       *
       * @return long
       */
      long myComputeSinglePseudoBettiNumber(long deg, long numWedges) const;

      /**
       * compute all possible wedge products of NonMultVars of.
       * Length specified the length of the DynamicBitsets
       *
       * @param NonMultVars
       * @param length
       *
       * @return a vector of DynamicBitsets which represents all possible wedges
       */
      std::vector<DynamicBitset> myPossibleWedgesOfLength(const std::vector<long>& NonMultVars, long length) const;

      /**
       * It computes the relevant part of the 'constant resolution': First it computes the
       * relevant part of the basic constant graph  and then it applies the constant morse
       * reduction.
       *
       * @param row
       * @param col
       */
      void myComputeConstantColumnResolution(long deg, long numWedges);

      /**
       * It computes the pseudo betti number \beta_{numWedges, numWedges+deg} for all deg >= MinDeg
       * up to MaxDeg. MaxDeg is the maximal degree which can occur (in this setting the element
       * which has the highest degree in myResolution)
       * For this it uses the formula to compute the pseudo betti numbers out of the
       * induced free resolution (e.g. see Book Involution)
       *
       * @param MinDeg
       * @param numWedges
       *
       * @return a vector of longs
       */
      std::vector<long> myComputePartialPseudoBettiNumber(long MinDeg, long numWedges) const;

      /**
       * It computes the ranks of the submatrices which representing the concerning submatrices
       * of the differentials. The concerning matrices are all matrices which map from deg to deg
       * for MinDeg <= deg <= LenWasteRanks + MinDeg from the module F_{numWedges + 1} to F_{numWedges}
       * and F_{numWedges} to F_{numWedges - 1}. It adds the ranks of the matrices concerning to deg
       * and insert it at position deg - MinDeg in the resulting vector.
       * If numWedges = 0 the it only takes a look at the maps from F_{numWedges + 1} to F_{numWedges}.
       *
       * @param LenWasteRanks
       * @param MinDeg
       * @param numWedges
       * @return a vector of longs
       */
      std::vector<long> myComputePartialRanks(long LenWasteRanks, long MinDeg, long numWedges) const;

      /**
       * It computes the ranks of the submatrices which representing the concerning submatrices
       * of the differentials. The concerning matrices are the two matrices which map from deg to deg
       * from the module F_{numWedges + 1} to F_{numWedges} and F_{numWedges} to F_{numWedges - 1}.
       * It adds the ranks of the matrices concerning to deg and returns it.
       * If numWedges = 0 the it only takes a look at the map from F_{numWedges + 1} to F_{numWedges}.
       *
       * @param deg
       * @param numWedges
       *
       * @return long
       */
      long myComputeSingleRank(long deg, long numWedges) const;

      /**
       * It computes the ranks of the submatrices which representing the differentials. The
       * concerning matrices are the matrices which map from PosInRes + 1 to PosInRes
       * for degree deg, where deg >= MinDeg. The shape of the matrices is given in RowRanks[deg-MinDeg]
       * resp. ColRanks[deg-MinDeg]. It fills the matrices with values from myResolution. The start point
       * is given by ResIter. At the end it computes the rank of the matrix and insert (deg, rank) into
       * the result vector.
       *
       *
       * @param ResInIter
       * @param ColRanks
       * @param PosInRes
       * @param MinDeg
       * @return a vector of pair<long,long>
       */
      std::vector<std::pair<long, long> > myComputePartialRanksPerDegree(std::map<MorseElement, MorsePaths>::const_iterator& ResIter,
                                                                              const std::vector<long>& RowRanks,
                                                                              const std::vector<long>& ColRanks,
                                                                              long PosInRes,
                                                                              long MinDeg) const;

      /**
       * It computes the ranks of the submatrix which representing the concerning differential. The
       * concerning matrix is the matrix which map from PosInRes + 1 to PosInRes
       * for degree deg. The shape of the matrix is given in RowRank
       * resp. ColRank. It fills the matrix with values from myResolution. The start point
       * is given by ResIter. At the end it computes the rank of the matrix and return the value.
       *
       * @param RowRank
       * @param ColRank
       * @param PosInRes
       * @param deg
       * @return a long
       */
      long myComputeSingleRankPerDegree(std::map<MorseElement, MorsePaths>::const_iterator& ResIter,
                                             long RowRank,
                                             long ColRank,
                                             long PosInRes,
                                             long deg) const;

      /**
       * This method takes a vector of longs and transform it to a matrix with one column and
       * len(ranks) rows. The i-th row in the matrix corresponds to ranks[i].
       *
       * @param ranks
       * @return a matrix of type RingZZ
       */
      matrix myTransformPartialRanksToPartialBettis(const std::vector<long>& ranks) const;


    };
  } // end namespace Involutive
} // end namespace CoCoA

#endif
