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

#include "CoCoA/TmpMorseElement.H"

#include "CoCoA/SparsePolyIter.H"

#include <sstream>

using std::make_pair;

namespace CoCoA
{
  namespace Involutive
  {

    long myMaxTrueOfBitset(DynamicBitset bitset)
    {
      const long length(bitset.myLen());
      for (long i = length - 1; i > -1; --i)
      {
        if (bitset.Iam1At(i))
        {
          return i;
        }
      }
      return -1;
    }


    DynamicBitset myVectorLongToDynamicBitset(const std::vector<long>& longs, const long& length)
    {
      DynamicBitset result(length);
      for (std::vector<long>::const_iterator it = longs.begin(); it != longs.end(); ++it)
      {
        result.mySet(*it);
      }
      return result;
    }

    DynamicBitset myVectorBoolToDynamicBitset(const std::vector<bool>& BoolVector)
    {
      const long n = len(BoolVector);
      DynamicBitset result(n);
      for (long i=0; i < n; ++i)
        result.mySet(i, BoolVector[i]);
      return result;
    }

    std::vector<long> myDynamicBitsetToLong(const DynamicBitset& DynBitset)
    {
      std::vector<long> result;
      long length(len(DynBitset));
      for (long i = 0; i < length; ++i)
      {
        if (DynBitset.Iam1At(i))
        {
          result.push_back(i);
        }
      }
      return result;
    }

    const std::vector<RingElem>& StandardRepresentationContainer::myComputeStandardRepresentation(ConstRefRingElem r)
    {
      ++myOriginNormalForms; // statistic
      ContainerIter res(myContainer.end());
      // extracts range of container with equal LPP
      std::pair<ContainerIter, ContainerIter> iters(myContainer.equal_range(LPP(r)));

      // checks if r is in this range
      for (ContainerIter it = iters.first; it != iters.second; ++it)
      {
        if ((it->second).first == r)
        {
          res = it;
          break;
        }
      }
      // r is completly new -> compute inv Standard Representation
      if (res == myContainer.end())
      {
        ++myReallyNormalForms; // statistic
        res = myContainer.insert(make_pair(LPP(r), make_pair(r, myMill.myStandardRepresentationWithoutRestShort(r))));
      }
      // returns inv standard rep
      return (res->second).second;
    }


    //------------------------------------------------------------------
    // MorseElement

    // 2 ctors
    MorseElement::MorseElement(const DynamicBitset& WedgeProduct, const JBElemConstIter basis)
        : myWedgeProduct(WedgeProduct)
        , myRightFactor(LPP(one(owner(basis->elem))))
        , myBasis(basis)
        , myRightProduct(LPP(myBasis->elem))
        , myProduct(myRightProduct * NewPP(owner(myRightFactor), myWedgeProduct))
    {
      myWedgeProductAsLongs = myDynamicBitsetToLong(myWedgeProduct);
    }

    MorseElement::MorseElement(const DynamicBitset& WedgeProduct, const PPMonoidElem& RightFactor, const JBElemConstIter basis)
        : myWedgeProduct(WedgeProduct)
        , myRightFactor(RightFactor)
        , myBasis(basis)
        , myRightProduct(LPP(myBasis->elem) * myRightFactor)
        , myProduct(myRightProduct * NewPP(owner(myRightFactor), myWedgeProduct))
    {
      myWedgeProductAsLongs = myDynamicBitsetToLong(myWedgeProduct);
    }


    void MorseElement::mySetWedgeProduct(const DynamicBitset& elem)
    {
      myProduct = (myProduct / NewPP(owner(myRightFactor), myWedgeProduct)) * NewPP(owner(myRightFactor), elem);
      myWedgeProduct = elem;
      myWedgeProductAsLongs = myDynamicBitsetToLong(myWedgeProduct);
      CoCoA_ASSERT(myProduct == NewPP(owner(myRightFactor), myWedgeProduct) * myRightFactor * LPP(myBasis->elem));
    }


    void MorseElement::mySetRightFactor(const PPMonoidElem& elem)
    {
      myProduct = (myProduct / myRightFactor) * elem;
      myRightProduct = (myRightProduct / myRightFactor) * elem;
      myRightFactor = elem;
      CoCoA_ASSERT(myRightProduct == LPP(myBasis->elem) * myRightFactor);
      CoCoA_ASSERT(myProduct == NewPP(owner(myRightFactor), myWedgeProduct) * myRightFactor * LPP(myBasis->elem));
    }


    void MorseElement::myDivideRightProductWith(long i)
    {
      const PPMonoidElem t = indet(owner(myProduct), i);
      myProduct = myProduct / t;
      myRightProduct = myRightProduct / t;
      myRightFactor = myRightFactor / t;
      CoCoA_ASSERT(myRightProduct == LPP(myBasis->elem) * myRightFactor);
      CoCoA_ASSERT(myProduct == NewPP(owner(myRightFactor), myWedgeProduct) * myRightFactor * LPP(myBasis->elem));
    }


    bool operator <(const MorseElement& m1, const MorseElement& m2)
    {
      const long l1 = len(m1.myWedgeProductAsLongs);
      const long l2 = len(m2.myWedgeProductAsLongs);
      if (l1 != l2) return (l1 < l2);
      // Same length...
      const int CmpProduct(owner(m1.myProduct)->myCmp(raw(m1.myProduct), raw(m2.myProduct)));
      if (CmpProduct != 0) return (CmpProduct < 0);

      if (m1.myBasis->lexPos != m2.myBasis->lexPos) return m1.myBasis->lexPos > m2.myBasis->lexPos;

      const int CmpRightProduct(owner(m1.myProduct)->myCmp(raw(m1.myRightProduct), raw(m2.myRightProduct)));
      return (CmpRightProduct > 0);
    }

    bool operator <=(const MorseElement& m1, const MorseElement& m2) { return !(m2 < m1); }
    bool operator > (const MorseElement& m1, const MorseElement& m2) { return  (m2 < m1); }
    bool operator >=(const MorseElement& m1, const MorseElement& m2) { return !(m1 < m2); }


    long MorseElement::myEpsilon(long test, long add) const
    {
      CoCoA_ASSERT(test < len(myWedgeProduct));
      CoCoA_ASSERT(add < len(myWedgeProduct));

      bool even = true;
      if ((add < test) && !myWedgeProduct.Iam1At(add))
        even = !even;
      for (long i = 0; i < test; ++i)
      {
        even ^= myWedgeProduct.Iam1At(i);
      }
      if (even) return 1;
      else return -1;
    }


    long MorseElement::myMaxTypeOne() const
    {
      return myMaxTypeOne(myRightFactor, myWedgeProduct);
    }

    long MorseElement::myMaxTypeOne(const PPMonoidElem& supp, const DynamicBitset& wedgeProduct) const
    {
      const DynamicBitset support(supp);
      return myMaxTrueOfBitset(support - wedgeProduct);
    }


    long MorseElement::myMaxTypeTwo() const
    {
      return myMaxTypeTwo(myRightFactor, myWedgeProduct, myBasis->multVars);
    }

    long MorseElement::myMaxTypeTwo(const PPMonoidElem& supp, const DynamicBitset& wedgeProduct, const DynamicBitset& nCritVars) const
    {
      const DynamicBitset support(supp);
      return myMaxTrueOfBitset((support | wedgeProduct) & nCritVars);
    }


    std::vector<std::pair<MorseElement, RingElem> > MorseElement::myComputeBasicMaps(const std::pair<JBElemConstIter, JBElemConstIter>& BasisIters, StandardRepresentationContainer& container) const
    {
      std::vector< std::pair<MorseElement, RingElem> > result;
      const long length(len(myWedgeProduct));
      // remove one index from wedge product
      for (std::vector<long>::const_iterator it = myWedgeProductAsLongs.begin(); it != myWedgeProductAsLongs.end(); ++it)
      {
        // initialize new wedge product and compute maps
        DynamicBitset NewWedge(myWedgeWithOneRemoved(myWedgeProductAsLongs, length, it));
        myComputeLeftMap(result, *it, NewWedge);
        myComputeRightMaps(result, *it, NewWedge, BasisIters, container, myGetPolyRing());
      }
      return result;
    }


    std::vector<std::pair<MorseElement, RingElem> > MorseElement::myComputeBasicConstantMaps(const std::pair<JBElemConstIter, JBElemConstIter>& BasisIters, StandardRepresentationContainer& container) const
    {
      std::vector<std::pair<MorseElement, RingElem> > result;
      const long length(len(myWedgeProduct));
      // remove one index from wedge product
      for (std::vector<long>::const_iterator it = myWedgeProductAsLongs.begin(); it != myWedgeProductAsLongs.end(); ++it)
      {
        // initialize new wedge product and compute maps
        DynamicBitset NewWedge(myWedgeWithOneRemoved(myWedgeProductAsLongs, length, it));
        myComputeRightMaps(result, *it, NewWedge, BasisIters, container, CoeffRing(myGetPolyRing()));
      }
      return result;
    }


    DynamicBitset MorseElement::myWedgeWithOneRemoved(const std::vector<long>& WedgeAsLong, long LengthWedge, std::vector<long>::const_iterator it) const
    {
      DynamicBitset ans = myVectorLongToDynamicBitset(WedgeAsLong, LengthWedge);
      ans.mySet(*it, false);
      return ans;
    }


    void MorseElement::myComputeLeftMap(std::vector<std::pair<MorseElement, RingElem> >& maps, long i, DynamicBitset NewWedge) const
    {
      MorseElement m(NewWedge, myBasis);
      const SparsePolyRing P = owner(myBasis->elem);
      maps.push_back(make_pair(m, myEpsilon(i,i)*indet(P,i)));
    }


    void MorseElement::myComputeRightMaps(std::vector<std::pair<MorseElement, RingElem> >& maps,
                                          long i,
                                          DynamicBitset NewWedge,
                                          const std::pair<JBElemConstIter, JBElemConstIter>& BasisIters,
                                          StandardRepresentationContainer& container,
                                          const ring& MapRing)  const
    {
      const SparsePolyRing P = owner(myBasis->elem);
      bool constantMap(MapRing != P);

      // The order of the JB is here important... Lex Order is wrong... we need DegRevLexOrder!!!!!!
      // compute standard representation
      const std::vector<RingElem>& SecondPart(container.myComputeStandardRepresentation(myBasis->elem * indet(P,i)));
      JBElemConstIter BasisIter(BasisIters.first);

      // iterate over inv standard rep.
      for (std::vector<RingElem>::const_iterator RepIter = SecondPart.begin(); RepIter != SecondPart.end(); ++RepIter)
      {
        // split factor of JBElem in terms
        for (SparsePolyIter SPI = BeginIter(*RepIter); !IsEnded(SPI); ++SPI)
        {
          // if below is true this will reduces to zero -> we can omit this
          if (!IsOne(PP(SPI)) || !IsDisjoint(NewWedge, BasisIter->multVars))
          {
            if (IsSubset(NewWedge, BasisIter->multVars) && !NewWedge.IamAll0s())
            {
              continue;
            }
            if (constantMap && exponent(PP(SPI), NumIndets(P) - 1) != 0)
            {
              continue;
            }
            long max(myMaxTypeOne(PP(SPI), NewWedge));
            if ( max == -1)
            {
              continue;
            }
            if (myMaxTypeTwo(PP(SPI), NewWedge, BasisIter->multVars) != max)
            {
              continue;
            }
          }
          // create new MorseElement
          MorseElement m(NewWedge, PP(SPI), BasisIter);
          // compute new Map
          RingElem map(MapRing, (-myEpsilon(i, i)) * coeff(SPI));
          maps.push_back(make_pair(m, map));
        }
        ++BasisIter;
      }
      CoCoA_ASSERT(BasisIter == BasisIters.second);
    }

    std::string MorseElement::toStr() const
    {
      std::ostringstream os;
      for (std::vector<long>::const_iterator i = myWedgeProductAsLongs.begin(); i != myWedgeProductAsLongs.end(); ++i)
      {
        if(i != myWedgeProductAsLongs.begin())
        {
          os << "\u2227";
        }
        os << *i;
      }
      os << "\u2297" << "(" << myGetRightFactor() << ")(" << myGetBasisElement() << ")";
      os << " LexPos(" << myBasis->lexPos << ")";
      return os.str();
    }

  } // end of namespace Involutive
} // end of namespace CoCoA
