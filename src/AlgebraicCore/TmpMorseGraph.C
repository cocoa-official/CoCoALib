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

#include "CoCoA/TmpMorseGraph.H"

#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/matrix.H"

using std::make_pair;

namespace CoCoA
{
  namespace Involutive
  {
    std::vector<MorseElement> MorseGraph::myComputeGeneralBasis() const
    {
      std::vector<MorseElement> basis;
      for (MorseElement::JBElemConstIter it = myBasis.begin(); it != myBasis.end(); ++it)
      {
        const std::vector<long> NonMultVars(myDynamicBitsetToLong(flip(it->multVars)));
        const std::vector<DynamicBitset> bitsets(myPossibleWedges(NonMultVars));
        for (std::vector<DynamicBitset>::const_iterator BitsetIter = bitsets.begin(); BitsetIter != bitsets.end(); ++BitsetIter)
        {
          basis.push_back(MorseElement(*BitsetIter, it));
        }
      }
      return basis;
    }


    std::vector<DynamicBitset> MorseGraph::myPossibleWedges(const std::vector<long>& NonMultVars) const
    {
      std::vector<std::vector<long> > ResultAsVector;
      const long NumNonMultVars = len(NonMultVars);
      for (long i = 0; i <= NumNonMultVars; ++i)
      {
        myVariationWithoutRepetition(ResultAsVector, std::vector<long>(), NonMultVars, i);
      }
      std::vector<DynamicBitset> result; result.reserve(len(ResultAsVector));
      for (std::vector<std::vector<long> >::iterator it = ResultAsVector.begin(); it != ResultAsVector.end(); ++it)
      {
        result.push_back(myVectorLongToDynamicBitset(*it, NumIndets(myRing)));
      }
      return result;
    }


    void MorseGraph::myVariationWithoutRepetition(std::vector<std::vector<long> >& result,
                                                  const std::vector<long>& CurrentResult,
                                                  const std::vector<long>& InputSet,
                                                  long length)  const
    {
      if (length == 0)
      {
        result.push_back(CurrentResult);
        return;
      }
      for (std::vector<long>::const_iterator it = InputSet.begin(); it != InputSet.end(); ++it)
      {
        std::vector<long> TmpCurrentResult(CurrentResult);
        TmpCurrentResult.push_back(*it);
        std::vector<long>::const_iterator iter(it);
        ++iter;
        std::vector<long> TmpInputSet(iter, InputSet.end());
        myVariationWithoutRepetition(result, TmpCurrentResult, TmpInputSet, length - 1);
      }
    }


   void MorseGraph::myAddMapsToResolution(const std::vector<std::pair<MorseElement, RingElem> >& maps, const ConstResIter& origin)
    {
      for (std::vector<std::pair<MorseElement, RingElem> >::const_iterator MapIter = maps.begin(); MapIter != maps.end(); ++MapIter)
      {
        myResolution[MapIter->first].myAddPath(origin, MapIter->second);
      }
    }


    std::map<MorseElement, MorsePaths>::reverse_iterator MorseGraph::myDeleteAndJumpToPrevious(std::map<MorseElement, MorsePaths>::reverse_iterator iter)
    {
      // We want to erase a_3
      // a_1 a_2 a_3 a_4
      //  -   -   i  i.b
      ++iter;
      // a_1 a_2 a_3 a_4
      //  -   i  i.b  -
      if (iter != myResolution.rend())
      {
        std::map<MorseElement, MorsePaths>::iterator TmpIter(iter.base());
        // a_1 a_2 a_3 a_4
        //  -   i  i.b  -
        //         Tmp
        ++iter;
        // a_1 a_2 a_3 a_4
        //  i  i.b Tmp  -
        myResolution.erase(TmpIter);
        // a_1 a_2 a_4
        //  i  i.b  -
        --iter;
        // a_1 a_2 a_4
        //  -   i  i.b
      } else {
        // iter == myResolution.rend()
        // --> we want to delete myResolution.begin()
        myResolution.erase(myResolution.begin());
        iter = myResolution.rend();
      }
      // base points to the next object in the sequence
      // e.g. to the original object before increment iter
      // something like this should work with C++11
      // iter = std::map<MorseElement, MorsePaths>::reverse_iterator(myResolution.erase(iter.base()));
      return iter;
    }

    RingElem MorseGraph::myCreateNewBasisElement(const MorseElement& m, long IndexMult, long IndexDiv) const
    {
      RingElem result(indet(myRing, IndexMult));
      result *= monomial(myRing, one(CoeffRing(myRing)), m.myGetRightFactor());
      result = result / indet(myRing, IndexDiv);
      result *= m.myGetBasisElement();
      return result;
    }


    void MorseGraph::myRightMinimization(const PathMap& paths,  DynamicBitset NewWedgeProduct, const MorseElement& origin, long maximum, long LongIter, StandardRepresentationContainer& container)
    {
      if (LongIter != maximum && (!origin.IamMultiplicativeIndex(LongIter)))
      {
        // std::cout << "Right minimize " << origin << "to ";
        // std::cout << "(Swap " << maximum << " with " << LongIter << "):" << std::endl;

        const bool InPolyRing(myMapRing == myRing);
        // NewBasis = (origin * x_LongIter)/x_maximum
        // RingElem NewBasis(myCreateNewBasisElement(origin, LongIter, maximum));
        const std::vector<RingElem>& NormalForm(container.myComputeStandardRepresentation(myCreateNewBasisElement(origin, LongIter, maximum)));
        MorseElement::JBElemConstIter BasisIter(myBasis.begin());
        for (std::vector<RingElem>::const_iterator NormalFormIter = NormalForm.begin(); NormalFormIter != NormalForm.end(); ++NormalFormIter)
        {
          CoCoA_ASSERT(BasisIter != myBasis.end());
          const SparsePolyIter& EndIter(myRing->myEndIter(raw(*NormalFormIter)));
          for (SparsePolyIter SPI = myRing->myBeginIter(raw(*NormalFormIter)); SPI != EndIter; ++SPI)
          {
            if (!IsOne(PP(SPI)) || !IsDisjoint(NewWedgeProduct, BasisIter->multVars))
            {
              if (IsSubset(NewWedgeProduct, BasisIter->multVars) && !(origin.myGetWedgeProduct()).IamAll0s())
              {
                continue;
              }
              long max(origin.myMaxTypeOne(PP(SPI), NewWedgeProduct));
              if ( max == -1 || origin.myMaxTypeTwo(PP(SPI), NewWedgeProduct, BasisIter->multVars) != max)
              {
                continue;
              }
              if (!InPolyRing && exponent(PP(SPI), NumIndets(myRing) - 1) != 0)
              {
                continue;
              }
            }
            MorseElement NewMorseElement(NewWedgeProduct, PP(SPI), BasisIter);
            // std::cout << "    " << NewMorseElement << std::endl;
            for (PathMap::const_iterator PathIter = paths.begin(); PathIter != paths.end(); ++PathIter)
            {
              RingElem path(PathIter->second);
              if (InPolyRing)
              {
                SparsePolyRingPtr(myMapRing)->myMulByCoeff(raw(path), raw(coeff(SPI)));
              } else {
                myMapRing->myMul(raw(path), raw(path), raw(coeff(SPI)));
              }
              if (origin.myEpsilon(maximum, maximum) * origin.myEpsilon(LongIter, maximum) == 1)
              {
                myMapRing->myNegate(raw(path), raw(path));
              }
              // std::cout << "    maximum = " << maximum << " LongIter = " << LongIter << std::endl;
              // std::cout << "    origin.myEpsilon(maximum, maximum) = " << origin.myEpsilon(maximum, maximum) << std::endl;
              // std::cout << "    origin.myEpsilon(LongIter, maximum) = " << origin.myEpsilon(LongIter, maximum) << std::endl;
              // std::cout << "    coeff(SPI) = " << coeff(SPI) << std::endl;
              // std::cout << "    PathIter->second = " << PathIter->second << std::endl;
              // std::cout << "    Add Path " << path << " to " << (PathIter->first)->first << std::endl;
              CoCoA_ASSERT(origin > NewMorseElement);
              myResolution[NewMorseElement].myAddPath(PathIter->first, path);
              // std::cout << "adding path (" << path << ") for " << NewMorseElement << std::endl;
              // myPrintMorseElements();
            }
          }
          ++BasisIter;
        }
        // std::cout << std::endl;
      }
    }


    std::string MorseGraph::toStr() const
    {
      std::ostringstream os;
      os << std::endl << std::endl << std::endl;
      os << "----------------------" << std::endl;
      os << "----------------------" << std::endl;
      os << "CURRENT MORSE ELEMENTS" << std::endl;
      for (std::map<MorseElement, MorsePaths>::const_iterator i = myResolution.begin(); i != myResolution.end(); ++i)
      {
        os << "---> g = " << i->first << " ----" << std::endl;
        PathMap maps((i->second).myGetPaths());
        for (PathMap::const_iterator map = maps.begin(); map != maps.end(); ++map)
        {
          os << "       " << (map->first)->first << " ---(" << map->second << ")--> g" << std::endl;
        }
        os << std::endl;
      }
      return os.str();
    }

    bool MorseGraph::myLexCompareJBElems(std::vector< std::pair<RingElem, std::vector<bool> > >::const_iterator e1,
                                         std::vector< std::pair<RingElem, std::vector<bool> > >::const_iterator e2)
    {
      PPMonoidElem e1LPP(LPP(e1->first));
      PPMonoidElem e2LPP(LPP(e2->first));
      const std::vector<long> e1Exps = exponents(e1LPP);
      const std::vector<long> e2Exps = exponents(e2LPP);

      return std::lexicographical_compare(e1Exps.begin(), e1Exps.end(), e2Exps.begin(), e2Exps.end());
    }

    std::vector<MorseElement::JBElem> MorseGraph::myTransform(const JBMill& mill) const
    {
      std::vector< std::pair<RingElem, std::vector<bool> > > basis(mill.myMultVarsWithRingElem());
      std::vector<std::vector< std::pair<RingElem, std::vector<bool> > >::const_iterator> lexSortBasis;
      for (std::vector< std::pair<RingElem, std::vector<bool> > >::const_iterator it = basis.begin(); it != basis.end(); ++it)
      {
        lexSortBasis.push_back(it);
      }





      std::sort(lexSortBasis.begin(), lexSortBasis.end(), MorseGraph::myLexCompareJBElems);
      // for (std::vector<std::vector< std::pair<RingElem, std::vector<bool> > >::const_iterator>::iterator i = lexSortBasis.begin(); i != lexSortBasis.end(); ++i)
      // {
      //   std::cout <<  "(*i)->first = " << (*i)->first << std::endl;
      // }
      std::vector<MorseElement::JBElem> res;
      for (std::vector< std::pair<RingElem, std::vector<bool> > >::const_iterator it = basis.begin(); it != basis.end(); ++it)
      {
        std::vector<std::vector< std::pair<RingElem, std::vector<bool> > >::const_iterator>::iterator match(std::lower_bound(lexSortBasis.begin(), lexSortBasis.end(), it, MorseGraph::myLexCompareJBElems));
        long disToSmallest(std::distance(lexSortBasis.begin(), match));
        // std::cout <<  "it->first = " << it->first << std::endl;
        // std::cout <<  "disToSmallest = " << disToSmallest << std::endl;
        // std::cout <<  "==================== " << std::endl;
        res.push_back(MorseElement::JBElem( it->first,
                              myVectorBoolToDynamicBitset(it->second),
                              disToSmallest));
      }
      return res;
    }

    void MorseGraph::mySetBetaVector()
    {
      // first entry degree
      // second entry number of multVars
      for (std::vector<MorseElement::JBElem>::iterator i = myBasis.begin(); i != myBasis.end(); ++i)
      {
        myBetaVector[make_pair(deg(i->elem), count(i->multVars))] += 1;
      }

      // for (std::map<std::pair<long,long>, long>::iterator i = myBetaVector.begin(); i != myBetaVector.end(); ++i)
      // {
      //   std::cout << "(" << (i->first).first <<", " << (i->first).second << ") = " << i->second << std::endl;
      // }
    }
  } // end of namespace Involutive
} // end of namespace CoCoA
