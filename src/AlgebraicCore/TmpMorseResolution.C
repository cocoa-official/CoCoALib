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

#include "CoCoA/TmpMorseResolution.H"

#include "CoCoA/matrix.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/TmpResolutionMinimization.H"

using std::make_pair;

namespace CoCoA
{
  namespace Involutive
  {

    std::vector<matrix> MorseResolution::myMapsAsMatrices() const
    {
      std::vector<matrix> res;
      std::map<MorseElement, MorsePaths>::const_iterator iter(myResolution.end());
      --iter;
      long size((iter->first).myCountWedgeBasis());
      for (long i = 0; i < size; ++i)
      {
        res.push_back(myMapsAsMatrix(i));
      }
      return res;
    }


    matrix MorseResolution::myMapsAsMatrix(long pos) const
    {
      matrix res = myInitialMapsAsMatrix(pos);
      std::map<MorseElement, long> identifier(myPositionMorseElements(pos));
      for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
      {
        if ((it->first).myCountWedgeBasis() == pos)
        {
          const long PositionRow(identifier[it->first]);
          PathMap maps((it->second).myGetPaths());
          for (PathMap::iterator MapIter = maps.begin(); MapIter != maps.end(); ++ MapIter)
          {
            const long PositionCol(identifier[(MapIter->first)->first]);
            SetEntry(res, PositionRow, PositionCol, MapIter->second);
          }
        }
      }
      return res;
    }

  // we can now access an element via <NumberWedges, Position In NumberWedges>
    std::map<MorseElement, long> MorseResolution::myPositionMorseElements(long NumberWedges) const
    {
      std::map<MorseElement, long> res;
      long pos1(0);
      long pos2(0);
      for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
      {
        if ((it->first).myCountWedgeBasis() == NumberWedges)
        {
          res.insert(make_pair(it->first, pos1));
          ++pos1;
        }
        if ((it->first).myCountWedgeBasis() == NumberWedges + 1)
        {
          res.insert(make_pair(it->first, pos2));
          ++pos2;
        }
      }
      return res;
    }

    matrix MorseResolution::myInitialMapsAsMatrix(long pos) const
    {
      long row(0);
      long col(0);
      for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
      {
        const long NumWedgeBasis((it->first).myCountWedgeBasis());
        if (pos > NumWedgeBasis)
        {
          continue;
        }
        if (pos == NumWedgeBasis)
        {
          ++row;
        }
        if ((pos + 1) == NumWedgeBasis)
        {
          ++col;
        }
        if ((pos + 1) < NumWedgeBasis)
        {
          break;
        }
      }
      return NewDenseMat(myMapRing, row, col);
    }

    void MorseResolution::myLeftMinimization(const PathMap& paths, DynamicBitset NewWedgeProduct, const MorseElement& origin, long maximum, long LongIter)
    {
      // MorseElement NewMorseElement(origin);
      // NewMorseElement.myDivideRightProductWith(maximum);
      PPMonoidElem newRightProduct(origin.myGetRightFactor() / indet(owner(origin.myGetRightFactor()), maximum));
      if (!IsOne(newRightProduct) || !IsDisjoint(NewWedgeProduct, origin.myGetNCrit()))
      {
        if (IsSubset(NewWedgeProduct, origin.myGetNCrit()) && !(origin.myGetWedgeProduct()).IamAll0s())
        {
          return;
        }

        long max(origin.myMaxTypeOne(newRightProduct, NewWedgeProduct));
        if ( max == -1)
        {
          return;
        }
        if (origin.myMaxTypeTwo(newRightProduct, NewWedgeProduct, origin.myGetNCrit()) != max)
        {
          return;
        }
      }

      MorseElement NewMorseElement(origin);
      NewMorseElement.mySetRightFactor(newRightProduct);
      NewMorseElement.mySetWedgeProduct(NewWedgeProduct);

      for (PathMap::const_iterator PathIter = paths.begin(); PathIter != paths.end(); ++PathIter)
      {
        myResolution[NewMorseElement].myAddPath(PathIter->first, PathIter->second * origin.myEpsilon(maximum, maximum) * origin.myEpsilon(LongIter, maximum) * indet(myRing, LongIter));
      }
    }

    matrix MorseResolution::myZerothMatrix() const
    {
      std::vector<RingElem> vec;
      for (std::map<MorseElement, MorsePaths>::const_iterator it = myResolution.begin(); it != myResolution.end(); ++it)
      {
        if (0 != (it->first).myCountWedgeBasis())
        {
          break;
        }
        vec.push_back((it->first).myGetBasisElement());
      }
      CoCoA_ASSERT(!vec.empty());
      return NewDenseMat(RowMat(vec));
    }

    std::vector<matrix> MorseResolution::myGetResolution() const
    {
      std::vector<matrix> res(1, myZerothMatrix());
      std::vector<matrix> OtherMaps(myMapsAsMatrices());
      res.insert(res.end(), OtherMaps.begin(), OtherMaps.end());
      return res;
    }

    void MorseResolution::myComputeBasicGraph(const std::vector<MorseElement>& elems, StandardRepresentationContainer& container)
    {
      for (std::vector<MorseElement>::const_iterator it = elems.begin(); it != elems.end(); ++it)
      {
        ConstResIter IterToRes(myResolution.insert(make_pair(*it, MorsePaths())).first);
        myAddMapsToResolution(it->myComputeBasicMaps(myGetBasisRange(), container), IterToRes);
      }
    }

    void MorseResolution::myDirectMorseReduction(StandardRepresentationContainer& container)
    {
      // myPrintMorseElements();
      std::map<MorseElement, MorsePaths>::reverse_iterator ResolutionIter(myResolution.rbegin());
      while (ResolutionIter != myResolution.rend())
      {
        if ((ResolutionIter->first).IAmBasisElement())
        {
          ++ResolutionIter;
          continue;
        }
        if ((ResolutionIter->second).IamEmpty())
        {
          ResolutionIter = myDeleteAndJumpToPrevious(ResolutionIter);
          continue;
        }
        MorseElement morse(ResolutionIter->first);
        const long maximum(morse.myMaxTypeOne());
        const PathMap& paths((ResolutionIter->second).myGetPaths());
        std::vector<long> longs(morse.myGetWedgeProductAsLongs());
        // because of myMaxTypeOne i is not in longs
        longs.push_back(maximum); // for JB this is not an ordered list anymore!
        for (std::vector<long>::iterator LongIter = longs.begin(); LongIter != longs.end(); ++LongIter)
        {
          DynamicBitset NewWedgeProduct(morse.myGetWedgeProduct());
          NewWedgeProduct.mySet(maximum, true);
          NewWedgeProduct.mySet(*LongIter, false);
          myLeftMinimization(paths, NewWedgeProduct, morse, maximum, *LongIter);
          myRightMinimization(paths, NewWedgeProduct, morse, maximum, *LongIter, container);
        }
        ResolutionIter = myDeleteAndJumpToPrevious(ResolutionIter);
        // myPrintMorseElements();
      }
    }


    void MorseResolution::myComputeResolution()
    {
      myResolution.clear();
      StandardRepresentationContainer container(myMill);
      myComputeBasicGraph(myComputeGeneralBasis(), container);
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "unreduced graph" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << *this << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "unreduced graph end" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      myDirectMorseReduction(container);
    }

    std::vector<matrix> MorseResolution::myComputeMinimalResolution()
    {
      //initialization
      myResolution.clear();
      StandardRepresentationContainer container(myMill);
      //compute resolution
      myComputeBasicGraph(myComputeGeneralBasis(), container);
      myDirectMorseReduction(container);
      //minimalization
      std::vector<matrix> maps(1, myZerothMatrix());
      std::vector<matrix> OtherMaps(myMapsAsMatrices());
      maps.insert(maps.end(), OtherMaps.begin(), OtherMaps.end());
      ResolutionMinimization minim(myRing, maps);
      minim.myMinimization();
      return minim.myGetResolution();
    }

    /*
     * computes the minimal free Resolution
     */
    std::vector<matrix> MinimalResolution(JBMill mill)
    {
      if (!IsStdDegRevLex(ordering(mill.myGetPPMonoid())))
      {
        CoCoA_THROW_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
      }
      if (!mill.IamHomogenous())
      {
        CoCoA_THROW_ERROR("The ideal isn't homogenous", "JBMinimalResolution");
      }
      MorseResolution mg(mill);
      return mg.myComputeMinimalResolution();

    }


    /*
     * computes a free Resolution
     */
    std::vector<matrix> Resolution(JBMill mill)
    {
      if (!IsStdDegRevLex(ordering(mill.myGetPPMonoid())))
      {
        CoCoA_THROW_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
      }

      MorseResolution mg(mill);
      mg.myComputeResolution();
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "reduced graph" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << mg << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "reduced graph end" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      // std::cout << "==========================================" << std::endl;
      std::vector<matrix> maps(1, mg.myZerothMatrix());
      std::vector<matrix> OtherMaps(mg.myMapsAsMatrices());
      maps.insert(maps.end(), OtherMaps.begin(), OtherMaps.end());
      return maps;
    }

  } // end of namespace Involutive
} // end of namespace CoCoa
