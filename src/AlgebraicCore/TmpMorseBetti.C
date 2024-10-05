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

#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/matrix.H"

using std::make_pair;

namespace CoCoA
{
  namespace Involutive
  {

    std::vector< std::vector<long> > MorseBetti::myComputeWasteRanks(const std::vector< std::vector<long> >& ranks) const
    {
      // extract length of the resolution
      // maybe easier out of ranks???
      std::map<MorseElement, MorsePaths>::const_iterator iter(myResolution.end());
      --iter;
      const long size((iter->first).myCountWedgeBasis());
      // initialize the waste ranks matrix. It has the same dimension as ranks
      std::vector< std::vector<long> > WasteRanks (len(ranks), std::vector<long>(len(ranks[0]), 0));
      // compute the waste ranks for each map
      // ResIter moves forward during the loop via myComputeWasteRanksPerMap
      std::map<MorseElement, MorsePaths>::const_iterator ResIter(myResolution.begin());
      for (long i = 0; i < size; ++i)
      {
        myComputeWasteRanksPerMap(WasteRanks, ResIter, ranks, i);
      }
      return WasteRanks;
    }


    void MorseBetti::myComputeWasteRanksPerMap(std::vector< std::vector<long> >& WasteRanks, std::map<MorseElement, MorsePaths>::const_iterator& ResIter, const std::vector< std::vector<long> >& ranks, long PosInRes) const
    {
      // computes the ranks of the degree matrices. Move ResIter Forward
      std::vector<std::pair<long, long> > RanksOfDegreeMatrices(myComputeWasteRanksPerDegree(ResIter, ranks[PosInRes], ranks[PosInRes + 1], PosInRes));
      // adept the ranks
      for (std::vector<std::pair<long, long> >::iterator RMIter = RanksOfDegreeMatrices.begin(); RMIter != RanksOfDegreeMatrices.end(); ++RMIter)
      {
        WasteRanks[PosInRes][RMIter->first] = WasteRanks[PosInRes][RMIter->first] + RMIter->second;
        long TooMuch(0);
        if (WasteRanks[PosInRes][RMIter->first] > ranks[PosInRes][RMIter->first])
        {
          TooMuch = WasteRanks[PosInRes][RMIter->first] - ranks[PosInRes][RMIter->first];
          WasteRanks[PosInRes][RMIter->first] = ranks[PosInRes][RMIter->first];
        }
        WasteRanks[PosInRes + 1][RMIter->first] = RMIter->second - TooMuch;
      }
    }

    long MorseBetti::myLastNonZeroIndex(const std::vector<long> vec) const
    {
      long counter = len(vec);
      for (std::vector<long>::const_reverse_iterator i = vec.rbegin(); i != vec.rend(); ++i)
      {
        if (*i != 0)
        {
          return counter;
        }
        --counter;
      }
      return counter;
    }

    long MorseBetti::myNumRowsBettiDiagram(const std::vector<std::vector<long> >& ranks) const
    {
      CoCoA_ASSERT(!ranks.empty());
      long res(0); //skipping 0 degree
      long currentRow(1);

      for (std::vector<std::vector<long> >::const_iterator i = ranks.begin(); i != ranks.end(); ++i)
      {
        long index(myLastNonZeroIndex(*i));
        if (index - currentRow > res)
        {
          res = index - currentRow;
        }
        ++currentRow;
      }

      if (res <= 0)
      {
        return 1;
      }
      return res;
    }

  // we can now access an element via <NumberWedges, Position In NumberWedges>
  // orderd by degree
    std::map<std::pair<long, MorseElement>, long> MorseBetti::myGradedPositionMorseElements(long NumberWedges) const
    {
      std::map<std::pair<long, MorseElement>, long> res;
      long pos(0);
      long degree(0);
      long OldNumberWedges(0);
      for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
      {
        const long CurrentNumberWedges((it->first).myCountWedgeBasis());
        if (CurrentNumberWedges < NumberWedges)
        {
          continue;
        }
        if (CurrentNumberWedges > NumberWedges + 1)
        {
          break;
        }
        if (degree != (it->first).myDegree() || CurrentNumberWedges != OldNumberWedges)
        {
          pos = 0;
          degree = (it->first).myDegree();
          OldNumberWedges = CurrentNumberWedges;
        }
        res.insert(make_pair(make_pair(degree, it->first), pos));
        ++pos;
      }
      return res;
    }

    matrix MorseBetti::myTransformRanksToBettis(const std::vector<std::vector<long> >& ranks) const
    {
      const long RowsBetti(myNumRowsBettiDiagram(ranks));
      const long NumCols = (myNumColsBettiDiagram(ranks));
      if (NumCols == 1 && len(ranks[0]) == 1)
      {
        return NewDenseMat(RingZZ(), 1, 1);
      }
      const long LenRanksRow = len(ranks[0]);
      matrix bettis(NewDenseMat(RingZZ(), RowsBetti, NumCols + 1));
      for (long i = 0; i < NumCols; ++i)
      {
        for (long j = i + 1; j < RowsBetti + i + 1; ++j)
        {
          if (j < LenRanksRow)
          {
            SetEntry(bettis, j - i - 1, i + 1, ranks[i][j]);
          }
        }
      }
      SetEntry(bettis, 0, 0, 1);
      const long nrows = NumRows(bettis);
      for (long i = 1; i < nrows; ++i)
      {
        SetEntry(bettis, i, 0, 0);
      }
      return bettis;
    }

    long MorseBetti::myNumColsBettiDiagram(const std::vector<std::vector<long> >& ranks) const
    {
      long NumCols = len(ranks);
      for (std::vector<std::vector<long> >::const_reverse_iterator i = ranks.rbegin(); i != ranks.rend(); ++i)
      {
        if (!(i->empty()))
        {
          bool AllZero = true;
          for (std::vector<long>::const_iterator j = (*i).begin(); j != (*i).end(); ++j)
          {
            if ((*j) != 0)
            {
              AllZero = false;
              break;
            }
          }
          if (!AllZero)
          {
            break;
          }
        }
        --NumCols;
      }
      return NumCols;
    }

    std::vector<std::vector<long> > MorseBetti::myMatrixMinus(std::vector<std::vector<long> > m1, const std::vector<std::vector<long> >& m2) const
    {
      const long nrows = len(m1);
      const long ncols = len(m1[0]);
      for (long i = 0; i < nrows; ++i)
      {
        for (long j = 0; j < ncols; ++j)
        {
          m1[i][j] -= m2[i][j];
        }
      }
      return m1;
    }

    std::vector<std::pair<long, long> > MorseBetti::myComputeWasteRanksPerDegree(std::map<MorseElement, MorsePaths>::const_iterator& ResIter, const std::vector<long>& RowRanks, const std::vector<long>& ColRanks, long PosInRes) const
    {
      std::vector<std::pair<long, long> > RanksOfDegreeMatrices;

      std::map<std::pair<long, MorseElement>, long> identifier(myGradedPositionMorseElements(PosInRes));
      // find a real submatrix (number of rows and cols are not zero)
      for (long degree = 0; degree < len(RowRanks); ++degree)
      {
        if (RowRanks[degree] == 0 || ColRanks[degree] == 0)
        {
          if (RowRanks[degree] != 0)
          {
            while((ResIter->first).myDegree() == degree && (ResIter->first).myCountWedgeBasis() == PosInRes)
            {
              ++ResIter;
            }
          }
          continue;
        }
        // extract matrix
        matrix SubMat(ConstructDegreeMatrix(ResIter, RowRanks[degree], ColRanks[degree], degree, PosInRes, identifier));
        // compute rank and add to result
        RanksOfDegreeMatrices.push_back(make_pair(degree, rk(SubMat)));
      }
      return RanksOfDegreeMatrices;
    }

    matrix MorseBetti::ConstructDegreeMatrix(std::map<MorseElement, MorsePaths>::const_iterator& ResIter, long rows, long cols, long degree, long PosInRes, const std::map<std::pair<long, MorseElement>, long>& identifier) const
    {
      matrix SubMat(NewDenseMat(myMapRing, rows, cols));
      while((ResIter->first).myDegree() == degree && (ResIter->first).myCountWedgeBasis() == PosInRes)
      {
        std::map<std::pair<long, MorseElement>, long>::const_iterator RowIter(identifier.find(make_pair(degree, ResIter->first)));
        if (RowIter != identifier.end())
        {
          const long PositionRow(RowIter->second);
          PathMap maps((ResIter->second).myGetPaths());
          for (PathMap::const_iterator MapIter = maps.begin(); MapIter != maps.end(); ++MapIter)
          {
            std::map<std::pair<long, MorseElement>, long>::const_iterator ColIter(identifier.find(make_pair(degree, (MapIter->first)->first)));
            if (ColIter != identifier.end())
            {
              const long PositionCol(ColIter->second);
              SetEntry(SubMat, PositionRow, PositionCol, MapIter->second);
            }
          }
        }
        ++ResIter;
      }
      return SubMat;
    }

    void MorseBetti::myComputeConstantResolution()
    {
      myResolution.clear();
      StandardRepresentationContainer container(myMill);
      myComputeBasicConstantGraph(myComputeGeneralBasis(), container);

      myConstantDirectMorseReduction(container);
    }

    void MorseBetti::myConstantDirectMorseReduction(StandardRepresentationContainer& container)
    {
      std::map<MorseElement, MorsePaths>::reverse_iterator ResolutionIter(myResolution.rbegin());
      while (ResolutionIter != myResolution.rend())
      {
        if ((ResolutionIter->first).IAmBasisElement())
        {
          ++ResolutionIter;
          continue;
        }
        MorseElement morse(ResolutionIter->first);
        //assert !(ResolutionIter->second).IAmEmpty()
        // (ResolutionIter->second) should always be false
        // (because when it is empty this must be an BasisElement, but than we
        // does not reach this line)
        if ((ResolutionIter->second).IamEmpty())
        {
          ResolutionIter = myDeleteAndJumpToPrevious(ResolutionIter);
          continue;
        }
        const PathMap& paths((ResolutionIter->second).myGetPaths());
        std::vector<long> longs(morse.myGetWedgeProductAsLongs());
        // because of myMaxTypeOne i is not in longs
        const long maximum(morse.myMaxTypeOne());
        longs.push_back(maximum);
        for (std::vector<long>::iterator LongIter = longs.begin(); LongIter != longs.end(); ++LongIter)
        {
          DynamicBitset NewWedgeProduct(morse.myGetWedgeProduct());
          NewWedgeProduct.mySet(maximum, true);
          NewWedgeProduct.mySet(*LongIter, false);
          myRightMinimization(paths,
                              NewWedgeProduct,
                              morse,
                              maximum,
                              *LongIter,
                              container);
        }
        ResolutionIter = myDeleteAndJumpToPrevious(ResolutionIter);
      }
    }

    matrix MorseBetti::myComputeBettiNumbers()
    {
      myComputeConstantResolution();
      std::vector< std::vector<long> > ranks(myComputeRanks());
      std::vector< std::vector<long> > WasteRanks(myComputeWasteRanks(ranks));
      ranks = myMatrixMinus(ranks, WasteRanks);
      return myTransformRanksToBettis(ranks);
    }

    std::vector< std::vector<long> > MorseBetti::myComputeRanks() const
    {
      const long LengthRows = ((myResolution.rbegin())->first).myCountWedgeBasis() + 1;
      long LengthCols(0);
      std::vector< std::pair<long, long> > RanksVector;
      for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
      {
        RanksVector.push_back(make_pair((it->first).myCountWedgeBasis(), (it->first).myDegree()));
        if ((it->first).myDegree() > LengthCols)
        {
          LengthCols = (it->first).myDegree();
        }
      }
      std::vector< std::vector<long> > ranks (LengthRows, std::vector<long>(LengthCols + 1, 0));
      for (std::vector<std::pair<long, long> >::iterator it = RanksVector.begin(); it != RanksVector.end(); ++it)
      {
        ++ranks[it->first][it->second];
      }
      return ranks;
    }

    std::pair<matrix, matrix> MorseBetti::myComputePseudoBettiNumbers()
    {
      myComputeConstantResolution();
      std::vector< std::vector<long> > ranks(myComputeRanks());
      std::vector< std::vector<long> > WasteRanks(myComputeWasteRanks(ranks));
      return make_pair(myTransformRanksToBettis(ranks), myTransformRanksToBettis(myMatrixMinus(ranks,WasteRanks)));
    }

    void MorseBetti::myComputeBasicConstantGraph(const std::vector<MorseElement>& elems, StandardRepresentationContainer& container)
    {
      for (std::vector<MorseElement>::const_iterator it = elems.begin(); it != elems.end(); ++it)
      {
        ConstResIter IterToRes(myResolution.insert(make_pair(*it, MorsePaths())).first);
        std::vector<std::pair<MorseElement, RingElem> > maps(it->myComputeBasicConstantMaps(myGetBasisRange(), container));
        myAddMapsToResolution(maps, IterToRes);
      }
    }

    matrix BettiDiagram(JBMill mill)
    {
      if (!IsStdDegRevLex(ordering(mill.myGetPPMonoid())))
      {
        CoCoA_THROW_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
      }
      if (!mill.IamHomogenous())
      {
        CoCoA_THROW_ERROR("The ideal isn't homogenous", "JBMinimalResolution");
      }
      MorseBetti mg(mill);
      return mg.myComputeBettiNumbers();
    }

    std::pair<matrix, matrix> PseudoBettiDiagram(JBMill mill)
    {
      if (!IsStdDegRevLex(ordering(mill.myGetPPMonoid())))
      {
        CoCoA_THROW_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
      }
      if (!mill.IamHomogenous())
      {
        CoCoA_THROW_ERROR("The ideal isn't homogenous", "JBMinimalResolution");
      }
      MorseBetti mg(mill);
      return mg.myComputePseudoBettiNumbers();
    }

  } // end of namespace Involutive
} // end of namespace CoCoA
