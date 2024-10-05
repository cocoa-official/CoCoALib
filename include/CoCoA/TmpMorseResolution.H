#ifndef CoCoA_TmpMorseResolution_H
#define CoCoA_TmpMorseResolution_H

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

#include "CoCoA/TmpMorseGraph.H"

namespace CoCoA
{
  namespace Involutive
  {

    class MorseResolution : public MorseGraph {
    public:
      MorseResolution(JBMill mill): MorseGraph(mill)
      {
      }

      /**
       * Transforms the morse graph, represented by a map<MorseElement, MorsePaths>
       * to a vector of matrices
       *
       * @return a vector of matrices
       */
      std::vector<matrix> myMapsAsMatrices() const;

      /**
       * computes a free resolution but without minimalization
       */
      void myComputeResolution();

      /**
       * computes the minimal free resolution
       *
       * @return a vector of matrices, which represents the free resolution
       */
      std::vector<matrix> myComputeMinimalResolution();

      /*
       * returns the 0th matrix e.g the pommaret basis
       *
       * @return a matrix of RingElems
       */
      matrix myZerothMatrix() const;
    private:
      /**
       * Transforms the morse pathes which represent F_pos to a matrix
       *
       * @param pos
       * @return a matrix which represents F_pos
       */
      matrix myMapsAsMatrix(long pos) const;

      /**
       * The list of MorseElements are ordered in a strict way. Now we extract all
       * element with NumberWedges or NumberWedges+1 Wedge Basis elements and divide
       * them by number of wedges elements. Then we determine the position of a
       * specific MorseElement within all MorseElements with same number of wedge elements
       * and save this position in the value of the map.
       *
       * @param NumberWedges
       *
       * @return a map with MorseElement as key and position as described above as value.
       */
      std::map<MorseElement, long> myPositionMorseElements(long NumberWedges) const;

      /**
       * The resolution shall expressed as matrices. This method
       * determine the dimension of the matrix which corresponds
       * to the map F_pos
       *
       * @param pos
       * @return matrix with zeros, but correct number of rows and cols
       */
      matrix myInitialMapsAsMatrix(long pos) const;

      /**
       * Morse Reduction the left part:
       * This function reduces our MorseGraph with so called reductions of type 1 or simply left reductions.
       * It first computes a new right product, which is origin.rightFactor divided by x_maximum.
       * After that it does some sanity checks. If they fail we skip the left minimization for this element.
       * If they are valid we create a new MorseElement with NewWedgeProduct, the new rightFactor and the
       * basis element of origin.
       * Then we iter over all paths which ending at origin and adept them.
       *
       * @param paths
       * @param NewWedgeProduct
       * @param origin
       * @param maximum
       * @param LongIter
       */
      void myLeftMinimization(const PathMap& paths,
                              DynamicBitset NewWedgeProduct,
                              const MorseElement& origin,
                              long maximum,
                              long LongIter);

      /**
       * The method converts the MorsePaths to matrices, which represent the free resolution.
       *
       * @return a vector of matrices which represent the resolution.
       */
      std::vector<matrix> myGetResolution() const;

      /**
       * This methods takes all elementary MorseElements (MorseElements of the form WedgeProduct \tensor (1) (BasisElement))
       * and compute all maps.
       *
       * @param elems
       * @param container StandardRepresentationConatainer to speed up the computation and avoid redundant normal form computations
       */
      void myComputeBasicGraph(const std::vector<MorseElement>& elems, StandardRepresentationContainer& container);

      /**
       * Compute the Morse reduction of our resolution. It applies to every element which is not of the form
       * WedgeProduct \tensor (1) (BasisElement) the left and right minimization.
       *
       * @param container StandardRepresentationConatainer to speed up the computation and avoid redundant normal form computations
       */
      void myDirectMorseReduction(StandardRepresentationContainer& container);
    };

  } // end of namespace Involutive
} // end of namespace CoCoa
#endif

