#ifndef CoCoA_JBMill_H
#define CoCoA_JBMill_H

//   Copyright (c)  2015-2016  John Abott,  Anna M. Bigatti
//   Original author: 2015-2016  Mario Albert

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

#include "CoCoA/TmpJBAlgorithm.H"
#include "CoCoA/FreeModule.H"
#include <map>

namespace CoCoA
{
  namespace Involutive
  {
    // different strategys for computing Janet Basis
    enum StrategyFlag {TQDegree, TQBlockHigh, TQBlockLow, GBCompletion};

    /**
     * JBMill represents an interface to work with Janet bases represented by a JanetContainer,
     * JBMill offers several functions to analyse a given Janet basis. In addition to that
     * JBMill is the base class for the class PBMill which represents an interface to work
     * with Pommaret bases (Janet bases in delta-regular coordinate). It must be noted that
     * JBMill is not designed to compute a Janet basis. For this the user should use
     * the algorithms which are derived from JBAlgorithm or the internal class JBMill::Builder.
     */
    class JBMill
    {
    public:
      class Builder;

      /**
       * Constructors.
       */
      /**
       * Constructor which constructs a JBMill out of
       * a JanetContainer
       * @param basis a JanetContainer
       */
      JBMill(const JanetContainer& basis)
        : myBasis(basis),
          myPolyRing(basis.myGetPolyRing()),
          myPPMValue(PPM(basis.myGetPolyRing()))
      {}


      /**
       * Constructor which constructs a JBMill out of
       * a JanetContainer and a SparsePolyRing
       * @param basis a JanetContainer
       * @param polyRing a SparsePolyRing
       */
      JBMill(const JanetContainer& basis, SparsePolyRing polyRing)
        : myBasis(basis),
          myPolyRing(polyRing),
          myPPMValue(PPM(myPolyRing))
      {}

      /**
       * Copy-Constructor
       *
       * @param mill
       */
      JBMill(const JBMill& mill) : myBasis(mill.myGetJanetContainer()),
                                   myPolyRing(mill.myGetSparsePolyRing()),
                                   myPPMValue(PPM(myPolyRing))
      {}

      /**
       * Constructs a JBMill out of a JBMill::Builder object.
       *
       * @param builder
       */
      JBMill(const JBMill::Builder& builder);

      /**
       * Returns the minimal Janet Basis
       * @return The method returns a vector of RingElem's
       */
      std::vector<RingElem> myReturnJB() const;

      /**
       * Returns the reduced Groebner basis
       * @return The method returns a vector of RingElem's
       */
      std::vector<RingElem> myReturnGB() const;

      /**
       * This method prints the multiplicative variables of every element
       * in the given Janet basis
       */
      void myPrintMultVar() const;

      /**
       * This method prints the nonmultiplicative variables of every element
       * in the given Janet basis
       */
      void myPrintNonMultVar() const;

      /**
       * This method returns the maximal degree of an element in the Janet Basis
       */
      long myDeg() const;

      /**
       * This method computes the multiplicative variables of the given Janet
       * basis
       *
       * @return The method returns a key-value-map, where the key is a PPMonoidElem
       * (the LPP of an element in the Janet basis) and the value is a vector of bools, which
       * represents the multiplicative variables of the corresponding basis element.
       */
      std::map<PPMonoidElem, std::vector<bool> > myMultVars() const;

      /**
       * This method computes the nonmultiplicative variables of the given Janet
       * basis
       * @return The method returns a key-value-map, where the key is a PPMonoidElem
       * (the LPP of an element in the Janet basis) and the value is a vector of bools, which
       * represents the nonmultiplicative variables of the corresponding basis element.
       */
      std::map<PPMonoidElem, std::vector<bool> > myNonMultVars() const;

      /**
       * This method computes the nonmultiplicative variables of the given Janet
       * basis
       *
       * @return The method returns a vector of pairs. The first entry of a pair is the
       * corresponding polynomial of the Janet Basis. The second entry is a vector of bools, which
       * represents the nonmultiplicative variables of the corresponding basis element.
       */
      std::vector< std::pair<RingElem, std::vector<bool> > > myNonMultVarsWithRingElem() const;

      /**
       * This method computes the multiplicative variables of the given Janet
       * basis
       *
       * @return The method returns a vector of pairs. The first entry of a pair is the
       * corresponding polynomial of the Janet Basis. The second entry is a vector of bools, which
       * represents the multiplicative variables of the corresponding basis element.
       */
      std::vector< std::pair<RingElem, std::vector<bool> > > myMultVarsWithRingElem() const;

      /**
       * This method computes the mnonultiplicative variables of a given basis elem.
       * If the input is not contained in the basis we assume that every variable is
       * nonmultiplicative
       *
       * @param elem a RingElem which should be contained in the Janet basis
       *
       * @return a vector of bools, where true indicates that a variable is nonmultiplicative
       */
      std::vector<bool> myNonMultVarsOf(ConstRefRingElem elem) const;


      /**
       * This method checks if the Janet basis is also a Pommaret basis
       * @return It returns true if the Janet basis is also a Pommaret basis
       */
      bool IamPommaretBasis() const;

      /*
       * This method checks if the ideal (resp. Janet basis) is homogeneous
       * @return It returns true if the ideal (resp. Janet basis) is homogeneous
       */
      bool IamHomogenous() const;

      /*
       * This method checks if the Janet basis is monomial
       * @return It returns true if the Janet basis is monomial
       */
      bool IamMonomialIdeal() const;


      /**
       * This method computes an involutive standard representation of the element r
       * with respect to our given Janet-basis
       *
       * @param r a RingElem for that we shall compute the standard representation
       *
       * @return It returns a pair. The second entry is a RingElem which represents the remainder.
       * The first entry is a key-value-map, where the key is a PPMonoidElem
       * (the LPP of an element in the Janet basis) and the value is a RingElem such that
       * this element times the corresponding basis element appears in the standard representation.
       */
      std::pair<std::map<PPMonoidElem, RingElem>, RingElem> myStandardRepresentation(ConstRefRingElem r) const;

      /**
       * This method computes an involutive standard representation of the element r
       * with respect to our given Janet-basis. The method do not returns the remainder.
       *
       * @param r a RingElem for that we shall compute the standard representation
       *
       * @return  The method returns a vector of pairs. The first entry of a pair is the
       * corresponding polynomial of the Janet basis. The second entry  is a RingElem such that
       * this element times the corresponding basis element appears in the standard representation.
       */
      std::vector<std::pair<RingElem, RingElem> > myStandardRepresentationWithoutRest(ConstRefRingElem r) const;

      /**
       * This method computes an involutive standard representation of the element r
       * with respect to our given Janet-basis. The method do not returns the remainder.
       *
       * @param r a RingElem for that we shall compute the standard representation
       *
       * @return  The method returns a vector of RingElem. The i-th entry corresponds to the
       * i-th basis element of the Janet Basis (sorted by the given term order), that is
       * the i-th basis element times the i-th entry in the vector appears as summand in
       * the standard representation.
       */
      std::vector<RingElem> myStandardRepresentationWithoutRestShort(ConstRefRingElem r) const;


      /**
       * This method computes an involutive standard representation of the element r
       * with respect to our given Janet-basis. The result is printed to cout
       *
       * @param r a RingElem for that we the method shall compute the standard representation
       */
      void myOutputStandardRepresentation(ConstRefRingElem r) const;

      /**
       *Computes the Hilbert polynomial hp_R/I.
       *
       * @param s a RingElem which must be a single variable.
       * The Hilbert polynomial is expressed in this variable.
       *
       * @result a RingElem (a polynomial) in variable s
       */
      RingElem myHilbertPol(ConstRefRingElem s) const;

      /**
       * Computes the Hilbert function h_R/I(s)       *
       * @param s argument for the Hilbert function as BigInt
       * @return a BigInt
       */
      BigInt myHilbertFunc(const BigInt& s) const;

      /**
       * Computes the Hilbert Function. Output is the functional expression
       * The output is printed to cout.
       */
      void myHilbertFunc() const;


      /**
       * Computes the rational function of the Hilbert Series
       *
       * @param s a RingElem which must be a single variable.
       * The Ring of this element must be a fraction field.
       *
       * @return a RingElem (a rational function) in variable s.
       */
      RingElem myHilbertSeries(ConstRefRingElem s) const;

      /**
       * Computes the first involutive syzygy.
       *
       * @return a finitely generated module which represents the
       * first syzygy. The generators of this module are a Janet basis
       * of the syzygy.
       */
      FGModule mySyzygy() const;

      /*
       * Dimension of the P/I where I is generated by the current Janet-Basis
       *
       * @return the dimension as long.
       */
      long myDim() const;

      /**
       * Computes the class of an PPMonoidElem in the setting.
       * The class starts counting at 0 up to n - 1.
       * The cls of '1' is '-1'
       * cls(x^\mu) = max{i | \deg_i{x^\mu} != 0}
       *
       * @param elem a PPMonoidElem
       * @return returns the cls of elem as long
       */
      long myCls(ConstRefPPMonoidElem elem) const;

      /**
       * The minimal class of all LPP's of the Janet basis
       *
       * @return returns the minimal class as long
       */
      long myMinCls() const;

      /**
       * The maximal class of all LPP's of the Janet basis
       *
       * @return returns the maximal class as long
       */
      long myMaxCls() const;

      /**
       * All elements of the Janet basis where the class of the LPP
       * is InputCls
       *
       * @param InputCls
       *
       * @return returns all elements with the same class as vector of RingElems
       */
      std::vector<RingElem> myElementsWithClass(long InputCls) const;

      /**
       * Computes the complementary decomposition of I.
       * In this case I must be a monomial ideal.
       *
       * @return a vector of pairs. The first entry of a pair represents
       * a PPMonoidElem and the second one represents it corresponding multiplicative variables
       */
      std::vector< std::pair<PPMonoidElem, std::vector<bool> > > myComplementaryDecomposition() const;

      /**
       * Computes the complementary decomposition of lt(I).
       *
       * @return a vector of pairs. The first entry of a pair represents
       * a PPMonoidElem and the second one represents it corresponding multiplicative variables
       */
      std::vector< std::pair<PPMonoidElem, std::vector<bool> > > myComplementaryDecompositionPolynomial() const;

      /*
       * Computes the standard pairs of a monomial ideal
       *
       * @return a vector of pairs, which represents the standard pairs.
       */
      std::vector< std::pair<PPMonoidElem, std::vector<bool> > > myStandardPairs() const;

      /**
       * Computes the involutive normal form of a polynomial
       *
       * @param elem a RingElem.
       *
       * @return a RingElem, which is the involutive normal form of the input
       */
      RingElem myJNormalForm(const RingElem& elem) const;

      /**
       * Computes the j-divisor of a leading term of a polynomial.
       * It is only a wrapper to the 'deep' method
       *
       * @param elem a RingElem
       *
       * @result a RingElem which leading term is an involutive divisor of the leading
       * term of the input
       */
      inline RingElem myJDivisor(const RingElem& elem) const
      {
        if (myBasis.IamConstant())
        {
          if (IsZero(elem))
          {
            return zero(myPolyRing);
          }
          return one(myPolyRing);
        }
        JanetTriple* DivPtr(myBasis.myGetTree().myJDivisor(LPP(elem)));
        if(DivPtr == nullptr)
        {
          return zero(myPolyRing);
        }
        return DivPtr->myGetPol();
      }

      /*
      * noether normalization (Robertz: Noether normalization guided by monomial cone decompositions, page 8/9)
      * !!!!!NOT TESTED!!!!!
      * Only works in a few cases...
      */
      // std::pair<RingHom, std::vector<bool> > myNoetherNormalization() const;

      /**
       * Returns the SparsePolyRing of the JBMill
       *
       * @return a SparsePolyRing
       */
      SparsePolyRing myGetSparsePolyRing() const
      {
        return myPolyRing;
      }

      /**
       * Returns the PPMonoid of the JBMill
       *
       * @return a PPMonoid
       */
      PPMonoid myGetPPMonoid() const
      {
        return myPPMValue;
      }

      /**
       * Returns a const ref to the JanetContainer containing the Janet basis
       *
       * @return a const ref to the internal JanetContainer
       */
      const JanetContainer& myGetJanetContainer() const
      {
        return myBasis;
      }

    protected:
      /**
       * This method count the true values in a given
       * boolean vector
       *
       * @param vec a vector of bools
       * @return number of true's in the vector
       */
      long myCountTrues(const std::vector<bool>& vec) const;

      /**
       * Pommaret multiplicative variables for pp: cls(pp)..(n-1)
       * @param pp
       *
       * @return vector of bools which represents the multiplicative variables
       */
      std::vector<bool> myPommaretMultVar(ConstRefPPMonoidElem pp) const;

      JanetContainer myBasis;
      SparsePolyRing myPolyRing;
      PPMonoid myPPMValue;

    private:

      //For C++11 there is a more elegant way to avoid using this operator!!!
      // JBMill& operator=(const JBMill& mill); //never used

      /**
       * disjunction of two boolean vectors of same length:
       * res[i] = v1[i] or v2[i]
       *
       * @param v1 boolean vector
       * @param v2 boolean vector
       *
       * @return boolean vector who is disjunction of input
       */
      // std::vector<bool> myUnionBoolVectors(std::vector<bool> v1, std::vector<bool> v2) const;

      /**
       * Extract the maximal index of a true value wrt the term order.
       * Returns the indet of the polynomial ring with this index.
       *
       * @param v a vector of bools which represents the multiplicative variables
       * @return a RingElem which is an indet
       */
      // RingElem myGreatestMultVar(std::vector<bool> v) const;

      /*
      * check if mon a multiple of the multiplicative vars
      */

      // bool IamMultipleOfMultVars(PPMonoidElem ppm, std::vector<bool> v) const;

      /*
       * we need this method for PrintMultVar and PrintNonMultVar. It returns the (non)multiplicative variables
       */


      /**
       * Printing the multiplicative/nonmultiplicative variables
       *
       * @param multVars
       * @param outputMultVar -> multiplicative or nonmultiplicative?
       */
      void myOutputVar(const std::map<PPMonoidElem,std::vector<bool> >& multVars, bool outputMultVar) const;


      /**
       * Iterates recursively through the Janet Tree to determine the
       * multiplicative variables.
       *
       * @param multVars key value map containing the result (call by ref!)
       * @param iter current position in the Janet Tree
       * @param currentNonMultVars indicies of current multiplicative vars wrt iter
       * @param curVar the current variable
       */
      void myRekComputeMultVar(std::map<PPMonoidElem,std::vector<bool> >& multVars, JanetIterator iter, std::vector<int> currentNonMultVars, int curVar) const;


      /**
       * like the binomial coeffient but with rings
       * warning: maybe infinite loop!!!!
       * ((s n)) = (s * (s - 1) * ... * (s - n + 1)) / n!
       *
       * @param ring the polynomial ring
       * @param polAbove
       * @param intBelow
       * @return the binomial coefficient
       */
      RingElem myBinLike(PolyRing ring, RingElem polAbove, long intBelow) const;

      /**
       * Recursive part of the complementary decomposition.
       *
       * @param output
       * @param JIter
       */
      void myComplementaryDecompositionRecPart(std::vector< std::pair<PPMonoidElem, std::vector<bool> > >& output, JanetIterator JIter) const;

      /**
       * Computes a complementary decomposition of the leading ideal
       *
       * @return the complementary decomposition
       */
      std::vector< std::pair<PPMonoidElem, std::vector<bool> > > myComplementaryDecompositionLeadingIdeal() const;

      /**
       * Computes the multiplicative variables of a JanetIterator
       *
       * @param iter
       *
       * @return vector of bool where true indicates multiplicative
       */
      std::vector<bool> myMultVar(JanetIterator iter) const;

      /**
       * returns a reversed bool vector (true to false and vica versa)
       *
       * @param vec
       *
       * @return !vec
       */
      std::vector<bool> myReverseBoolVec(const std::vector<bool>& vec) const;

      /**
       * Computes the nonmultiplicative variables of a JanetIterator
       *
       * @param iter
       *
       * @return vector of bool where true indicates nonmultiplicative
       */
      std::vector<bool> myNonMultVar(JanetIterator iter) const;

      /**
       * Computes a free module for the first syzygy
       *
       * @return a free module
       */
      FreeModule myMakeNewFreeModuleForSyz() const;
    };

    /**
     * This class is designed to construct a Janet basis. The goal of this
     * class is to seperate the construction of the JBMill from its representation.
     * The 'Gang of Four' (Gamma, Helm, Johnson, Vlissides - Design Patterns) served
     * as template for the construction. The corresponding Pattern is called Building Pattern
     */
    class JBMill::Builder{
    public:
      friend class JBMill;

      // default value for the Strategy
      static const StrategyFlag defaultStrategy = TQBlockLow;
      // default value for the involutive criteria
      static const std::bitset<3> defaultInvCriteria;

      /**
       * Constructor: Sets the default Strategy and 2 involutive criteria as default
       */
      Builder() : myStrategy(defaultStrategy),
                  myInvolutiveCriteria(defaultInvCriteria)
      {
      }

      /**
       * Set the input
       *
       * @param input
       *
       * @return a reference to the current Builder object
       */
      Builder& setInput(const std::vector<RingElem>& input)
      {
        myInput = input;
        myRemoveZeroInput();
        return *this;
      }

      /**
       * Set the input
       *
       * @param inputBegin
       * @param inputEnd
       *
       * @return a reference to the current Builder object
       */
      Builder& setInput(std::vector<RingElem>::const_iterator inputBegin, std::vector<RingElem>::const_iterator inputEnd)
      {
        myInput.clear();
        myInput.insert(myInput.begin(), inputBegin, inputEnd);
        myRemoveZeroInput();
        return *this;
      }

      /**
       * Set the strategy to compute a Janet Basis
       *
       * @param strategy
       *
       * @return a reference to the current Builder object
       */
      Builder& setStrategy(StrategyFlag strategy)
      {
        myStrategy = strategy;
        return *this;
      }

      /**
       * Set the involutive criteria to compute a Janet Basis.
       * The involutive criteria are only applied in the TQ-Algos.
       *
       * @param invCriteria
       *
       * @return a reference to the current Builder object
       */
      Builder& setInvolutiveCriteria(const std::bitset<3> invCriteria)
      {
        myInvolutiveCriteria = invCriteria;
        return *this;
      }

      friend JBMill::JBMill(const JBMill::Builder& builder);

    private:
      std::vector<RingElem> myInput;
      StrategyFlag myStrategy;
      std::bitset<3> myInvolutiveCriteria;

      /**
       * Computes the a Janet basis for the given input with the given parameters.
       * Before computing the basis it performs some sanity checks.
       *
       * @return a JBMill
       */
      JBMill build()
      {
        myCheckEmptyInput();
        myCheckMixedRingElems();

        // return JBMill(myInput, myPolyRing(), myChooseStrategy());
        return JBMill(myUseStrategy());
      }

      /**
       * The task of this method is to return the choosen strategy as unique_ptr.
       * @return an unique_ptr to JBAlgorithm
       */
      std::unique_ptr<JBAlgorithm> myChooseStrategy() const;

      /**
       * Uses myStrategy to compute a JanetContainer (the internal representation)
       * of a Janet basis
       * @return JanetContainer
       */
      JanetContainer myUseStrategy() const;

      /**
       * Determine the SparsePolyRing of the input
       *
       * @return the sparse poly ring
       */
      SparsePolyRing myPolyRing() const
      {
        return owner(myInput.front());
      }

      /**
       * Checks if there is an empty input. If this is the case it raises an error
       */
      void myCheckEmptyInput() const;

      /**
       * Checks if there is an zero input. If this is the case it removes the zero
       * element from the input.
       */
      void myRemoveZeroInput();

      /**
       * Checks if the input consists only of elements of the same ring. If this
       * is the case it raises an error
       */
      void myCheckMixedRingElems() const;

    };
  }
}
#endif
