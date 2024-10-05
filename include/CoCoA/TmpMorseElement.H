#ifndef CoCoA_TmpMorseElement_H
#define CoCoA_TmpMorseElement_H

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

#include "CoCoA/DynamicBitset.H"
#include "CoCoA/TmpJBMill.H"
#include "CoCoA/utils.H"

namespace CoCoA
{
  namespace Involutive
  {
    // index of last "true" bit (or -1)
    long myMaxTrueOfBitset(DynamicBitset bitset);

    // conversions between DynamicBitset, vector<long> & vector<bool>
    DynamicBitset myVectorLongToDynamicBitset(const std::vector<long>& longs, const long& length);
    DynamicBitset myVectorBoolToDynamicBitset(const std::vector<bool>& BoolVector);
    std::vector<long> myDynamicBitsetToLong(const DynamicBitset& DynBitset);

    class MorseGraph;


    /* During the computation of the resolution/betti diagramm we have to compute many
     * involutive standard representations. It turns out, that we often compute invStrReps
     * for the same Elements. To avoid redundant computations we store an already computed
     * invStrRep in an Container. There we can first check it we already know this StrRep and
     * ommit an redundant computation.
     */
    class StandardRepresentationContainer
    {
    public:
      StandardRepresentationContainer(const JBMill& mill):
          myOriginNormalForms(0),
          myReallyNormalForms(0),
          myMill(mill)
      {}


      /*
       * First this function searches if we have already computed a StandardRepresentation for r.
       * If yes, it returns the matching StdRep; if not, it computes it, stores it in the container and
       * then returns it.
       */
      const std::vector<RingElem>& myComputeStandardRepresentation(ConstRefRingElem r);

      /*
       * Statistics
       */
      long myOriginNormalForms;
      long myReallyNormalForms;
    private: // data members
      typedef std::pair<RingElem, std::vector<RingElem> > StdRepr;
      typedef std::multimap<PPMonoidElem, StdRepr>::iterator ContainerIter;
      /*
       * stores all StandardRepresentations sorted after PPMonoidElem. But it is possible that
       * there are severals different RingElems with same LPP, therefore this is only a multimap
       * and we have to perform an additional search for the correct RingElem
       */
      std::multimap<PPMonoidElem, std::pair<RingElem, std::vector<RingElem> > > myContainer;

      /*
       * corresponding mill to compute StandardRepresenations
       */
      JBMill myMill;
    };

    /* GeneralDescription
     * The Morse Graph consists of many so called morse elements. A MorseElement consists
     * of an WedgeProduct, a
     * Basis Element of the current Janet Basis and a Factor. Beside the representation of
     * MorseElements the class
     * --> provides several functions to compute the corresponding paths between MorseElements
     * --> some helper functions to represent the data members of a MorseElement in several different structures
     * --> some member functions to compute signs
     *
     * The wedge product is an element of the exterior algebra /\V_i. It is represented as a Dynamic Bitset of length
     * n (number of variables) where i entries are true. An alternative representation of the
     * wedge product is a list of integers (the i indices which are true in the dynamic bitset)
     * The Factor is an simple power product and is represented by an PPMonoidElem.
     * The Factor must only consists of multiplicative vars w.r.t. JBElem.
     * The Basis Element of the current Janet Basis is represented as a JBElem (see below)
     *
     * Performance:
     * Instances of a MorseElements are used very often during a computation of a Resolution/BettiDiagramm.
     * Especially the comparison of MorseElements is used very often. To avoid redundant computations we
     * precomputed many datas to make the comparison of MorseElements more effective.
     * For more information read the comments above the data members of MorseElement.
     */
    class MorseElement
    {
    public:
      /*
       * The JanetBasis Element does not only consist of the basis element itself. It also
       * needs information of his multiplicative variables and his lexicographic position
       * in the list of all janet basis elements (the second fact is required for the ordering
       * of Morse Elements)
       */
      struct JBElem {
        /**
         * Constructor
         *
         * @param newElem
         * @param newMultVars
         * @param newLexPos
         */
        JBElem(RingElem newElem, DynamicBitset newMultVars, long newLexPos)
          : elem(newElem), multVars(newMultVars), lexPos(newLexPos)
        {
        }

        // The Basis Element
        RingElem elem;
        // The multiplicative variables of the basis element
        DynamicBitset multVars;
        // the lex Position in the list of janet basis elements:
        // b_1 <_lex b_2 => lexPos_b_1 < lexPos_b2
        long lexPos;
      };

      // The MorseElement does not save the JBElem itself. To avoid copying
      // we only store a reference to the position in the list of the Janet Basis
      typedef std::vector<JBElem>::const_iterator JBElemConstIter;

      /**
       * Constructor: Generates WedgeProduct \tensor (1) (basis)
       *
       * @param WedgeProduct
       * @param basis
       */
      MorseElement(const DynamicBitset& WedgeProduct, const JBElemConstIter basis);
      /**
       * Constructor: Generates WedgeProduct \tensor (RightFactor) (basis)
       *
       * @param WedgeProduct
       * @param RightFactor
       * @param basis
       */
      MorseElement(const DynamicBitset& WedgeProduct, const PPMonoidElem& RightFactor, const JBElemConstIter basis);

      /**
       * Returns myWedgeProduct
       *
       * @return myWedgeProduct as DynamicBitset
       */
      inline const DynamicBitset& myGetWedgeProduct() const
      { return myWedgeProduct; }

      /**
       * Returns myWedgeProduct as a vector of longs. If the i-th entry of myWedgeProduct is true
       * i is in myWedgeProductAsLongs
       *
       * @return myWedgeProduct as vector of longs, sorted in ascending order
       */
      inline const std::vector<long>& myGetWedgeProductAsLongs() const
      { return myWedgeProductAsLongs; }

      /**
       * Returns myRightFactor
       *
       * @return a PPMonoidElem
       */
      inline const PPMonoidElem& myGetRightFactor() const
      { return myRightFactor; }

      /**
       * Returns the noncritical variables of the MorseElement. These
       * are exactly the multiplicative variables of the basis element.
       *
       * @return noncritical variables as DynamicBitset
       */
      inline const DynamicBitset& myGetNCrit() const
      { return myBasis->multVars; }

      /**
       * Returns myBasis->elem
       *
       * @return RingElem
       */
      inline const RingElem& myGetBasisElement() const
      { return myBasis->elem; }

      /**
       * Sets the WedgeProduct. Adapting related data members (myProduct, myWedgeProductAsLongs)
       *
       * @param elem WedgeProduct as DynamicBitset
       */
      void mySetWedgeProduct(const DynamicBitset& elem);

      /**
       * Sets the RightFactor. Adapting related data members (myProduct)
       *
       * @param elem
       */
      void mySetRightFactor(const PPMonoidElem& elem);


      /**
       * test if MorseElement is an ordinary BasisElement
       * This is the case if the RightFactor is One &&
       * myWedgeProduct is a subset of flip(myNCrit)
       * last condition is equivalent to myWedgeProduct and
       * myNCrit are Disjoint
       *
       * @return true if MorseElement is an ordinary basis element
       */
      inline bool IAmBasisElement() const
      {
        return IsOne(myRightFactor) && IsDisjoint(myWedgeProduct, myBasis->multVars);
      }

      /**
       * max((support(myRightFactor) without WedgeProductAsLongs)
       *
       * @return result as long, -1 if there is no maximum
       */
      long myMaxTypeOne() const;

      /**
       * max((support(myRightFactor) without AsLongs(wedgeProduct))
       *
       * @param supp
       * @param wedgeProduct
       *
       * @return result as long, -1 if there is no maximum
       */
      long myMaxTypeOne(const PPMonoidElem& supp, const DynamicBitset& wedgeProduct) const;

      /**
       * max((support(myRightFactor) \/ myWedgeProduct) /\ myBasis->multVars)
       *
       * @return result as long, -1 if there is no maximum
       */
      long myMaxTypeTwo() const;

      /**
       * max((support(supp) \/ wedgeProduct) /\ nCritVars)
       *
       * @param supp
       * @param wedgeProduct
       * @param nCritVars
       *
       * @return result as long, -1 if there is no maximum
       */
      long myMaxTypeTwo(const PPMonoidElem& supp, const DynamicBitset& wedgeProduct, const DynamicBitset& nCritVars) const;

      /**
       * -1^max{i| (i in myWedgeProduct + add) and i < test}
       * assume test,add < len(myWedgeProduct)
       * assume add not in set (not necessary)
       *
       * @param test
       * @param add
       *
       * @return result as long
       */
      long myEpsilon(long test, long add) const;

      /**
       * Divides indeterminate x[i] from RightProduct
       * no check if this is possible!
       *
       * @param i
       */
      void myDivideRightProductWith(long i);


      /**
       * Computes the basic maps for a basic Morse Graph.
       * The maps are divided in left and right maps.
       * A map of an MorseElement points to an MorseElement again. But in this function
       * we only compute the wedge product of the new MorseElement. This product is simply
       * the same than the current, but with one 'true' entry converted to false.
       * We iterate over all possible new wedge products of this structure and computing the corresponding
       * right and left maps.
       * For computing the right and left maps we need the JanetBasis (BasisIters) and the StandardRepresentationContainer (container)
       * It returns the list of computed maps. (A computed map consists of the target and the map)
       *
       * @param BasisIters
       * @param container StandardRepresentationConatainer to speed up the computation and avoid redundant normal form computations
       *
       * @return a vector of pairs where the first element in a MorseElement (origin of the path) and the second one a RingElem (the value of the path)
       */
      std::vector<std::pair<MorseElement, RingElem> > myComputeBasicMaps(const std::pair<JBElemConstIter, JBElemConstIter>& BasisIters, StandardRepresentationContainer& container) const;

      /*
       * This function is nearly identical to myComputeBasicMaps, but omit to compute the left maps.
       * e.g. we only compute the constant part of the maps
       */
            /**
       * Computes the constant basic maps for a basic Morse Graph.
       * The maps are divided in left and right maps.
       * A map of an MorseElement points to an MorseElement again. But in this function
       * we only compute the wedge product of the new MorseElement. This product is simply
       * the same than the current, but with one 'true' entry converted to false.
       * We iterate over all possible new wedge products of this structure and computing the corresponding
       * maps. In contrast to myComputeBasicMaps we only need to compute the right maps, because the left maps only produce
       * non constant elements.
       * For computing the right maps we need the JanetBasis (BasisIters) and the StandardRepresentationContainer (container)
       * It returns the list of computed maps. (A computed map consists of the target and the map)
       *
       * @param BasisIters
       * @param container StandardRepresentationConatainer to speed up the computation and avoid redundant normal form computations
       *
       * @return a vector of pairs where the first element in a MorseElement (origin of the path) and the second one a RingElem (the value of the path)
       */
      std::vector<std::pair<MorseElement, RingElem> > myComputeBasicConstantMaps(const std::pair<JBElemConstIter, JBElemConstIter>& BasisIters, StandardRepresentationContainer& container) const;


      /**
       * returns the number of wedges in WedgeProduct
       *
       * @return number of wedges in myWedgeProduct
       */
      inline long myCountWedgeBasis() const
        { return len(myWedgeProductAsLongs); }

      /**
       * Returns the total degree of a MorseElement.
       * myCoundWedgeBasis + deg(myRightFactor) + deg(myBasis->elem)
       *
       * @return the result as long
       */
      inline long myDegree() const
      {
        return myCountWedgeBasis() + deg(myBasis->elem) + deg(myRightFactor);
      }

      /**
       * Checks if i is multiplicative for the basis element.
       *
       * @return true if i is multiplicative.
       */
      inline bool IamMultiplicativeIndex(long i) const
      { return (myBasis->multVars).Iam1At(i); }

      /**
       * Returns the MorseElement as string.
       *
       * @return a string
       */
       std::string toStr() const;
    private: // data members

      /**
       * Gets current SparsePolyRing
       * It is the owner of myBasisElement
       *
       * @return a poly ring.
       */
      inline SparsePolyRing myGetPolyRing() const
      { return owner(myBasis->elem); }

      /**
       * Computes the RightMaps of the Basic Graph
       * The map shall point to the new MorseElement which consists of NewWedge
       * the current JBElem and the right factor x_i. But x_i is in general nonmultiplicative.
       * Therefore this would not be a valid MorseElement. We compute the inv. StandardRep
       * of right factor times JBElem and split these to the form right_factor * JBElem.
       * Then we construct new MorseElements out of this with wedge product NewWedge.
       * The map to this new Element is simply the coefficient of right_factor in the StRep.
       *
       * If we only compute the betti numbers we know that all maps which might appear during the computation
       * are always in the base field. In this case we always save the maps as an element of the base field.
       *
       * In several cases we already know that a specific map will reduces later to zero. In this case
       * we omit this map.
       * The function returns nothing, but because of call by reference we add to the list of maps
       * the new maps.
       *
       * @param maps
       * @param i
       * @param NewWedge
       * @param BasisIters
       * @param container
       * @param MapRing
       */
      void myComputeRightMaps(std::vector<std::pair<MorseElement, RingElem> >& maps,
                              long i,
                              DynamicBitset NewWedge,
                              const std::pair<JBElemConstIter, JBElemConstIter>& BasisIters,
                              StandardRepresentationContainer& container,
                              const ring& MapRing) const;

      /**
       * Computes the LeftMap of the Basic Graph
       * The Map points to a MorseElement which consists of NewWedge
       * and the same JBElem than the current MorseElement
       * The map to the new Element is x_i multiplied by a sign.
       * The function returns nothing, but because of call by reference we add to the list of maps
       * the new map.
       *
       * @param maps
       * @param i
       * @param NewWedge
       */
      void myComputeLeftMap(std::vector<std::pair<MorseElement, RingElem> >& maps,
                             long i,
                             DynamicBitset NewWedge) const;

      /**
       * Constructs a new WedgeProduct with it in old WedgeAsLong removed
       *
       * @param WedgeAsLong
       * @param LengthWedge
       * @param it
       *
       * @return a DynamicBitset
       */
      DynamicBitset myWedgeWithOneRemoved(const std::vector<long>& WedgeAsLong, long LengthWedge, std::vector<long>::const_iterator it) const;

      friend bool operator <(const MorseElement& m1, const MorseElement& m2);
      friend bool operator ==(const MorseElement& m1, const MorseElement& m2);

      /*
       * A MorseElement consists of a WedgeProduct, a RightFactor, and a BasisElement
       */
      DynamicBitset myWedgeProduct;
      PPMonoidElem myRightFactor;
      JBElemConstIter myBasis;

      /*
       * for faster access
       * == myRightFactor * LPP(myBasis->elem)
       */
      PPMonoidElem myRightProduct;

      /*
       * for faster access
       * == AsPPMonoidElem(myWedgeProduct) * myRightFactor * LPP(myBasis->elem)
       */
      PPMonoidElem myProduct;

      /*
       * for faster access
       * myWedge as longs
       */
      std::vector<long> myWedgeProductAsLongs;
    }; // end of class MorseElement


    /**
     * Equal operator: Two MorseElements are equal if
     * myWedgeProduct, myRightFactor and myBasis are equal.
     *
     * @param m1
     * @param m2
     *
     * @return true if m1 and m2 are equal.
     */
    inline bool operator ==(const MorseElement& m1, const MorseElement& m2)
    {
      return (m1.myWedgeProduct == m2.myWedgeProduct) &&
        (m1.myRightFactor == m2.myRightFactor) &&
        (m1.myBasis->elem == m2.myBasis->elem); // Maybe Replacing me by LPP(m1.myBasis->elem)??
    }



    /* operator contains the following rules (old)
     * 1. NumTrues(m1) < NumTrues(m2)
     * 2. AsPPM(m1.myWedgeProduct) * m1.myRightFactor * LPP(m1.myBasisElement)
     *    < AsPPM(m2.myWedgeProduct) * m2.myRightFactor * LPP(m2.myBasisElement)
     * 3. m1.myNCrit contains m2.myNCrit
     * 4. The i-th index of myNCrit let be the i-th true entry.
     *    let i be the minimal index of m1.myNCrit and m2.myNcrit for which
     *    m1.myNCrit[i] != m2.myNCrit[i]. Then
     *    m1.myNCrit[i] < m2.myNCrit[i].
     * 5. m1.myRightFactor * LPP(m1.myBasisElement) < m2.myRightFactor * LPP(m2.myBasisElement)
     */
    /**
     * Operator <: operator contains the following rules (new)
     * 1. NumTrues(m1) < NumTrues(m2)
     * 2. m1.myProduct < m2.myProduct
     * 3. m1.myBasis >_lex m2.myBasis
     * 4. m1.myRightProduct > m2.myRightProduct
     *
     * @param m1
     * @param m2
     *
     * @return true if m1 < m2
     */
    bool operator <(const MorseElement& m1, const MorseElement& m2);

    /**
     * Operator <=
     *
     * @param m1
     * @param m2
     *
     * @return return true if m1 <= m2
     */
    bool operator <=(const MorseElement& m1, const MorseElement& m2);

    /**
     * Operator >
     *
     * @param m1
     * @param m2
     *
     * @return return true if m1 <= m2
     */
    bool operator >(const MorseElement& m1, const MorseElement& m2);

    /**
     * Operator >=
     *
     * @param m1
     * @param m2
     *
     * @return return true if m1 <= m2
     */
    bool operator >=(const MorseElement& m1, const MorseElement& m2);

    /**
     * MorseElement to stream
     *
     * @param os
     * @param obj
     *
     * @return stream object
     */
    inline std::ostream& operator<<(std::ostream& os, const MorseElement& obj)
    {
      os << obj.toStr();
      return os;
    }

  } // end of namespace Involutive
} // end of namespace CoCoA


#endif

