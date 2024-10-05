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

#include "CoCoA/TmpPartialMorseBetti.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/DenseMatrix.H"

using std::make_pair;

namespace CoCoA
{
  namespace Involutive
  {
    std::vector<MorseElement> PartialMorseBetti::myComputeGeneralColumnBasis(long MinDeg, long lengthWedge) const
    {
      std::vector<MorseElement> basis;
      for (MorseElement::JBElemConstIter it = myBasis.begin(); it != myBasis.end(); ++it)
      {
        long deg(StdDeg(it->elem));
        if (deg < MinDeg - lengthWedge - 1)
        {
          continue;
        } else if (deg == MinDeg - lengthWedge - 1)
        {
          myComputeBasisForElement(basis, it, lengthWedge + 1);
        } else {
          myComputeBasisForElement(basis, it, lengthWedge);
          myComputeBasisForElement(basis, it, lengthWedge + 1);
        }
      }
      return basis;
    }


    std::vector<MorseElement> PartialMorseBetti::myComputeGeneralSingleBasis(long degree, long lengthWedge) const
    {
      std::vector<MorseElement> basis;
      for (MorseElement::JBElemConstIter it = myBasis.begin(); it != myBasis.end(); ++it)
      {
        long deg(StdDeg(it->elem));
        if (deg == degree - lengthWedge - 1)
        {
          myComputeBasisForElement(basis, it, lengthWedge + 1);
        } else if (deg == degree - lengthWedge)
        {
          myComputeBasisForElement(basis, it, lengthWedge);
        }
      }
      return basis;
    }


    void PartialMorseBetti::myComputeBasisForElement(std::vector<MorseElement>& basis,
                                                     MorseElement::JBElemConstIter it,
                                                     long lengthWedge) const
    {
      const std::vector<long> NonMultVars(myDynamicBitsetToLong(flip(it->multVars)));
      if (len(NonMultVars) >= lengthWedge)
        {
        std::vector<DynamicBitset> bitsets(myPossibleWedgesOfLength(NonMultVars, lengthWedge));
        for (std::vector<DynamicBitset>::const_iterator BitsetIter = bitsets.begin(); BitsetIter != bitsets.end(); ++BitsetIter)
        {
          // std::cout <<  "MorseElement(*BitsetIter, it) = " << MorseElement(*BitsetIter, it) << std::endl;
          basis.push_back(MorseElement(*BitsetIter, it));
        }
      }
    }


    std::vector<DynamicBitset> PartialMorseBetti::myPossibleWedgesOfLength(const std::vector<long>& NonMultVars, long length) const
    {
      std::vector<std::vector<long> > ResultAsVector;
      myVariationWithoutRepetition(ResultAsVector, std::vector<long>(), NonMultVars, length);
      std::vector<DynamicBitset> result;
      result.reserve(len(ResultAsVector));
      for (std::vector<std::vector<long> >::iterator it = ResultAsVector.begin(); it != ResultAsVector.end(); ++it)
      {
        result.push_back(myVectorLongToDynamicBitset(*it, NumIndets(myRing)));
      }
      return result;
    }


    void PartialMorseBetti::myComputeConstantColumnResolution(long deg, long numWedges)
    {
      myResolution.clear();
      StandardRepresentationContainer container(myMill);
      myComputeBasicConstantGraph(myComputeGeneralColumnBasis(deg, numWedges), container);

      myConstantDirectMorseReduction(container);
    }


    void PartialMorseBetti::myComputeConstantSingleResolution(long deg, long numWedges)
    {
      myResolution.clear();
      StandardRepresentationContainer container(myMill);
      myComputeBasicConstantGraph(myComputeGeneralSingleBasis(deg, numWedges), container);

      myConstantDirectMorseReduction(container);
    }


    long PartialMorseBetti::myComputeBettiNumber(long row, long col)
    {
      if (col == 0)
      {
        if (row == 0)
        {
          return 1;
        }
        return 0;
      }
      long degree = row + col;
      long numWedges = col - 1;
      myComputeConstantSingleResolution(degree, numWedges);
      if (myResolution.empty())
      {
        return 0;
      }
      long rank(myComputeSinglePseudoBettiNumber(degree, numWedges));
      long WasteRank(myComputeSingleRank(degree, numWedges));
      return (rank - WasteRank);
    }


    matrix PartialMorseBetti::myComputeDownToBettiNumber(long row, long col)
    {
      if (col == 0)
      {
        matrix bettis(NewDenseMat(RingZZ(), 1, 1));
        if (row == 0)
        {
          SetEntry(bettis, 0, 0, 1);
        } else {
          SetEntry(bettis, 0, 0, 0);
        }
        return bettis;
      }
      long degree = row + col;
      long numWedges = col - 1;
      // std::cout <<  "degree = " << degree << std::endl;
      // std::cout <<  "numWedges = " << numWedges << std::endl;
      myComputeConstantColumnResolution(degree, numWedges);
      // std::cout <<  "len(myResolution) = " << len(myResolution) << std::endl;
      if (myResolution.empty())
      {
        matrix bettis(NewDenseMat(RingZZ(), 1, 1));
        SetEntry(bettis, 0, 0, 0);
        return bettis;
      }

      std::vector<long> ranks(myComputePartialPseudoBettiNumber(degree, numWedges));
      // std::cout <<  "ranks = ";
      // for (std::vector<long>::iterator i = ranks.begin(); i != ranks.end(); ++i)
      // {
      //   std::cout << *i << "; ";
      // }
      // std::cout << std::endl;
      std::vector<long> WasteRanks(myComputePartialRanks(len(ranks), degree, numWedges));
      // std::cout <<  "WasteRanks = ";
      // for (std::vector<long>::iterator i = WasteRanks.begin(); i != WasteRanks.end(); ++i)
      // {
      //   std::cout << *i << "; ";
      // }
      // std::cout << std::endl;
      CoCoA_ASSERT(len(ranks) == len(WasteRanks));
      for (int i = 0; i < len(ranks); ++i)
      {
        ranks[i] -= WasteRanks[i];
        if (ranks[i] < 0)
        {
          ranks[i] = 0;
        }
        // std::cout <<  "ranks[i] = " << ranks[i] << std::endl;
      }
      return myTransformPartialRanksToPartialBettis(ranks);
    }


    matrix PartialMorseBetti::myTransformPartialRanksToPartialBettis(const std::vector<long>& ranks) const
    {
      long RowsBetti(len(ranks));
      long TooMuch(0);
      for (int i = (RowsBetti - 1); i >= 0; --i)
      {
        if (ranks[i] != 0)
        {
          break;
        }
        ++TooMuch;
      }
      RowsBetti -= TooMuch;
      if (RowsBetti == 0)
      {
        RowsBetti = 1;
      }
      matrix bettis(NewDenseMat(RingZZ(), RowsBetti, 1));
      for (long i = 0; i < RowsBetti; ++i)
      {
        SetEntry(bettis, i, 0, ranks[i]);
      }
      return bettis;
    }


    std::vector<long> PartialMorseBetti::myComputePartialRanks(long LenWasteRanks, long MinDeg, long numWedges) const
    {
      // initialize the rank matrix. It has the same dimension as PseudoBettiNumber
      std::vector<long> WasteRanks(LenWasteRanks, 0);
      // compute the waste ranks for each map
      // ResIter moves forward during the loop via myComputeWasteRanksPerMap
      std::map<MorseElement, MorsePaths>::const_iterator ResIter(myResolution.begin());
      // std::cout <<  "MinDeg = " << MinDeg << std::endl;
      // std::cout <<  "numWedges = " << numWedges << std::endl;
      std::vector<std::pair<long, long> > WasteRanksFirst;
      if (numWedges > 0)
      {
        WasteRanksFirst = myComputePartialRanksPerDegree(ResIter,
                                                         myComputePartialPseudoBettiNumber(MinDeg, numWedges - 1),
                                                         myComputePartialPseudoBettiNumber(MinDeg, numWedges),
                                                         numWedges - 1,
                                                         MinDeg);
      }
      std::vector<std::pair<long, long> > WasteRanksSecond = myComputePartialRanksPerDegree(ResIter,
                                                                                            myComputePartialPseudoBettiNumber(MinDeg, numWedges),
                                                                                            myComputePartialPseudoBettiNumber(MinDeg, numWedges + 1),
                                                                                            numWedges,
                                                                                            MinDeg);
      for (std::vector<std::pair<long, long> >::iterator wrf = WasteRanksFirst.begin(); wrf != WasteRanksFirst.end(); ++wrf)
      {
        WasteRanks[wrf->first - MinDeg] = wrf->second;
      }

      for (std::vector<std::pair<long, long> >::iterator wrs = WasteRanksSecond.begin(); wrs != WasteRanksSecond.end(); ++wrs)
      {
        WasteRanks[wrs->first - MinDeg] += wrs->second;
      }

      return WasteRanks;
    }


    long PartialMorseBetti::myComputeSingleRank(long deg, long numWedges) const
    {

      std::map<MorseElement, MorsePaths>::const_iterator ResIter(myResolution.begin());
      long res = 0;
      if (numWedges > 0)
      {
        res = myComputeSingleRankPerDegree(ResIter,
                                           myComputeSinglePseudoBettiNumber(deg, numWedges - 1),
                                           myComputeSinglePseudoBettiNumber(deg, numWedges),
                                           numWedges - 1,
                                           deg);
      }

      res += myComputeSingleRankPerDegree(ResIter,
                                          myComputeSinglePseudoBettiNumber(deg, numWedges),
                                          myComputeSinglePseudoBettiNumber(deg, numWedges + 1),
                                          numWedges,
                                          deg);
      // std::cout <<  "(" << deg << "," << numWedges << ") = " << res << std::endl;
      return res;
    }


    long PartialMorseBetti::myComputeSingleRankPerDegree(std::map<MorseElement, MorsePaths>::const_iterator& ResIter,
                                                         long RowRank,
                                                         long ColRank,
                                                         long PosInRes,
                                                         long deg) const
    {
      if (RowRank == 0 || ColRank == 0)
      {
        return 0;
      }
      std::map<std::pair<long, MorseElement>, long> identifier(myGradedPositionMorseElements(PosInRes));
      matrix SubMat(ConstructDegreeMatrix(ResIter, RowRank, ColRank, deg, PosInRes, identifier));
      return rk(SubMat);
    }


    long PartialMorseBetti::myComputeSinglePseudoBettiNumber(long deg, long numWedges) const
    {
      long res(0);
      long n(NumIndets(myRing));
      for (long k = 1; k <= (n - numWedges); ++k)
      {
        BetaVector::const_iterator BetaIter = myBetaVector.find(make_pair(deg - numWedges, k));
        if (BetaIter != myBetaVector.end())
        {
          res += ConvertTo<long>(binomial(n - k, numWedges) * BetaIter->second);
        }
      }
      return AsUnsignedLong(res);
    }


    std::vector<long> PartialMorseBetti::myComputePartialPseudoBettiNumber(long MinDeg, long numWedges) const
    {
      const long highDeg(((myResolution.rbegin())->first).myDegree());
      std::vector<long> ranks(highDeg - MinDeg + 1, 0);
      for (int i = MinDeg; i <= highDeg; ++i)
      {
        ranks[i - MinDeg] = myComputeSinglePseudoBettiNumber(i, numWedges);
      }

      return ranks;
    }


    std::vector<std::pair<long, long> > PartialMorseBetti::myComputePartialRanksPerDegree(std::map<MorseElement, MorsePaths>::const_iterator& ResIter,
                                                                                          const std::vector<long>& RowRanks,
                                                                                          const std::vector<long>& ColRanks,
                                                                                          long PosInRes,
                                                                                          long MinDeg) const
    {
      std::vector<std::pair<long, long> > RanksOfDegreeMatrices;

      std::map<std::pair<long, MorseElement>, long> identifier(myGradedPositionMorseElements(PosInRes));
      // find a real submatrix (number of rows and cols are not zero)
      for (long pos = 0; pos < len(RowRanks); ++pos)
      {
        if (RowRanks[pos] == 0 || ColRanks[pos] == 0)
        {
          if (RowRanks[pos] != 0)
          {
            while((ResIter != myResolution.end()) &&(ResIter->first).myDegree() == (pos + MinDeg) && (ResIter->first).myCountWedgeBasis() == PosInRes)
            {
              ++ResIter;
            }
          }
          continue;
        }
        // extract matrix
        // std::cout <<  "RowRanks[pos] = " << RowRanks[pos] << std::endl;
        // std::cout <<  "ColRanks[pos] = " << ColRanks[pos] << std::endl;
        // std::cout <<  "pos+MinDeg = " << pos+MinDeg << std::endl;
        // std::cout <<  "PosInRes = " << PosInRes << std::endl;
        matrix SubMat(ConstructDegreeMatrix(ResIter, RowRanks[pos], ColRanks[pos], pos + MinDeg, PosInRes, identifier));
        // compute rank and add to result
        RanksOfDegreeMatrices.push_back(make_pair(pos + MinDeg, rk(SubMat)));
      }
      return RanksOfDegreeMatrices;
    }


    matrix BettiColumn(JBMill mill, long col)
    {
      if (!IsStdDegRevLex(ordering(mill.myGetPPMonoid())))
      {
        CoCoA_THROW_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
      }
      if (!mill.IamHomogenous())
      {
        CoCoA_THROW_ERROR("The ideal isn't homogenous", __FUNCTION__);
      }
      PartialMorseBetti mg(mill);
      return mg.myComputeBettiColumn(col);
    }


    matrix BettiPartialColumn(JBMill mill, long MinRow, long col)
    {
      if (!IsStdDegRevLex(ordering(mill.myGetPPMonoid())))
      {
        CoCoA_THROW_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
      }
      if (!mill.IamHomogenous())
      {
        CoCoA_THROW_ERROR("The ideal isn't homogenous", __FUNCTION__);
      }
      PartialMorseBetti mg(mill);
      return mg.myComputeDownToBettiNumber(MinRow, col);
    }

    RingElem BettiNumber(JBMill mill, long row, long col)
    {
      if (!IsStdDegRevLex(ordering(mill.myGetPPMonoid())))
      {
        CoCoA_THROW_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
      }
      if (!mill.IamHomogenous())
      {
        CoCoA_THROW_ERROR("The ideal isn't homogenous", __FUNCTION__);
      }
      PartialMorseBetti mg(mill);
      return RingElem(RingZZ(), mg.myComputeBettiNumber(row, col));
    }

  } // end namespace Involutive
} // end namespace CoCoA
