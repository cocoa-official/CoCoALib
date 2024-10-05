#ifndef CoCoA_UniversalInvolutiveBasisContainer_H
#define CoCoA_UniversalInvolutiveBasisContainer_H

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

#include "CoCoA/TmpPBMill.H"
#include "CoCoA/SmartPtrIRC.H"

#include <vector>
#include <string>
#include <map>
#include <memory>

namespace CoCoA
{
  namespace Involutive
  {

    /**
    * The UniverstalInvolutiveBasisContainer (In the following UIBC)
    * Goal of UIBC
    * The UIBC acts as an interface between JBMill/PBMill and an ideal. It handles
    * every request of an ideal to a JBMill/PBMill and computes if necessary a
    * Janet or Pommaret basis. In addition to that the UIBC stores every computed
    * value which is related to a Janet or Pommaret basis. Therefore it needs to
    * compute values like the dimension only once, because for additional calls it
    * already knows these values.
    *
    * Implementation of UIBC
    * The UIBC is designed as a subclass of IntrusiveReferenceCount. The reason for that
    * is that we use this class as a SmartPtrIRC. Nearly every datamember of UIBC is mutable one.
    * The only non-mutable one is gens, which is a vector of RingElems. gens is constant e.g.
    * after the initialization we can't change this value anymore. As all other data members of
    * UIBC are directly dependent from this generating set we make them mutable. In general
    * every mutable data member is 'empty'. Only if the user requests one of this values we check
    * if this value is non-empty and return it. If is empty we compute this value and return it.
    *
    * The class ``SparsePolyRing::IdealImpl`` contains a SmartPtrIRC<Involutive::UniversalInvolutiveBasisContainer>
    * as data member. It implements several methods to compute value via the UIBC. In addition to that
    * we have some non-member methods in the namespace Involutive, which access these methods.
    */
    class UniversalInvolutiveBasisContainer: protected IntrusiveReferenceCount
    {
      friend class SmartPtrIRC<const UniversalInvolutiveBasisContainer>;
      friend class SmartPtrIRC<UniversalInvolutiveBasisContainer>;

    public:
      typedef std::vector< std::pair<PPMonoidElem, std::vector<bool> > >  ComplementaryDecomposition;

      /**
       * Constructor, it sets the generators and inititalizes the numerical
       * values with -1.
       *
       * @param gens
       */
      UniversalInvolutiveBasisContainer(const std::vector<RingElem>& gens):
        IntrusiveReferenceCount(),
        generators(gens),
        strategy(TQBlockLow),
        degree(-1),
        dimension(-1),
        depth(-1),
        projDim(-1),
        regularity(-1),
        satiety(-1)
      {}

      /**
       * If the jbMillPtr is not present it builds the JBMill.
       * If jbMillPtr is present it does nothing.
       */
      void myInitializeJBMill() const;

      /**
       * First it initializes the JBMill and checks if the ideal is delta regular.
       * If this is the case and the pbMillPtr is not present it builds the PBMill.
       * If pbMillPtr is present it does nothing.
       * If the ideal is not delta regular it raises an error with fcnName
       *
       * @param fcnName
       */
      void myInitializeAndValidatePBMill(std::string fcnName) const;

      /**
       * It checks if the given ideal is homogeneous. If this is not the
       * case it raises an error with fcnName
       *
       * @param fcnName
       */
      void myValidateHomogeneous(std::string fcnName) const;

      /**
       * It checks if the given polynomial ring has DegRevLexOrdering. If this is not the
       * case it raises an error with fcnName
       *
       * @param fcnName
       */
      void myValidateDegRevLexOrdering(std::string fcnName) const;

      /**
       * It checks if the given ideal is monomial. If this is not the
       * case it raises an error with fcnName
       *
       * @param fcnName
       */
      void myValidateMonomialIdeal(std::string fcnName) const;

      /**
       * It checks if the given ideal is Cohen-Macualay. If this is not the
       * case it raises an error with fcnName
       *
       * @param fcnName
       */
      void myValidateCohenMacaulay(std::string fcnName) const;

      /**
       * It checks if the given ideal is delta-regular.
       * @return true if the ideal is delta-regular
       */
      bool IamDeltaRegular() const;

      /**
       * It checks if the given ideal is monomial.
       * @return true if the ideal is monomial
       */
      bool IamMonomial() const;

      /**
       * It checks if the given ideal is homogeneous.
       * @return true if the ideal is homogeneous
       */
      bool IamHomogeneous() const;

      /**
       * It compute the multiplicative variables of the Janet basis.
       * If multVars is not present it computes the multiplicative and
       * nonmultiplicative variables at once.
       * @return a map where the key is a PPMonoidElem which corresponds to the leading term
       * of a generator of the Janet basis, and the value is a vector of bools which represents
       * the multiplicative variables.
       */
      std::map<PPMonoidElem, std::vector<bool> > myMultVars() const;

      /**
       * It compute the nonmultiplicative variables of the Janet basis.
       * If multVars is not present it computes the nonmultiplicative and
       * nonmultiplicative variables at once.
       * @return a map where the key is a PPMonoidElem which corresponds to the leading term
       * of a generator of the Janet basis, and the value is a vector of bools which represents
       * the nonmultiplicative variables.
       */
      std::map<PPMonoidElem, std::vector<bool> > myNonMultVars() const;

      /**
       *Computes the Hilbert polynomial hp_R/I.
       *
       * @param s, a RingElem which must be a single variable.
       * The Hilbert polynomial is expressed in this variable.
       *
       * @result a RingElem (a polynomial) in variable s
       */
      RingElem myHilbertPol(ConstRefRingElem s) const;

      /**
       * Computes the rational function of the Hilbert Series
       *
       * @param s, a RingElem which must be a single variable.
       * The Ring of this element must be a fraction field.
       *
       * @return a RingElem (a rational function) in variable s.
       */
      RingElem myHilbertSeries(ConstRefRingElem s) const;

      /**
       * Computes the Hilbert Function h_R/I(s)       *
       * @param s argument for the hilbert function as BigInt
       * @return a BigInt
       */
      BigInt myHilbertFunc(const BigInt& s) const;

      /**
       * Computes the first involutive syzygy.
       *
       * @return a finitely generated module which represents the
       * first syzygy. The generators of this module are a Janet basis
       * of the syzygy.
       */
      FGModule myFirstSyzygy() const;

      /*
       * Dimension of the P/I where I is generated by the current Janet-Basis
       *
       * @return the dimension as long.
       */
      long myDimension() const;

      /**
       * Computes the complementary decomposition of lt(I).
       *
       * @return a vector of pairs. The first entry of a pair represents
       * a PPMonoidElem and the second one represents it corresponding multiplicative variables
       */
      ComplementaryDecomposition myComplementaryDecomposition() const;

      /**
       * Computes the depth of I
       *
       * @return the depth as long
       */
      long myDepth() const;

      /**
       * Computes the projective dimension of I
       *
       * @return the projective dimension as long
       */
      long myProjDim() const;

      /**
       * Computes the socle of I
       *
       * @return the socle as vector of RingElems
       */
      std::vector<RingElem> mySocle() const;

      /**
       * Computes the extremal betti numbers of I
       *
       * @return a map where the key is the position of the extremal betti number
       * as pair of longs and the value is the value of the betti diagram at this
       * position
       */
      std::map<std::pair<long, long>, long> myExtremalBettiNumbers() const;

      /**
       * Computes the maximal regular sequence of P/I
       *
       * @return the regular sequence as vector of RingElems
       */
      std::vector<RingElem> myRegularSequence() const;

      /**
       * Computes the unique maximal strongly independent of P/I
       *
       * @return the strongly independent set as vector of RingElems
       */
      std::vector<RingElem> myMaximalStronglyIndependentSet() const;

      /**
       * Computes the regularity of I
       *
       * @return the regularity as long
       */
      long myRegularity() const;

      /**
       * Computes the satiety of I
       *
       * @return the satiety as long
       */
      long mySatiety() const;

      /**
       * Computes the saturation of I
       *
       * @return the saturation as vector of RingEles
       */
      std::vector<RingElem> mySaturation() const;

      /**
       * Tests if I is Cohen Macualay
       *
       * @return true if I is Cohen Macaulay
       */
      bool IamCohenMacaulay() const;

      long myDegree() const;

      /**
       * Returns the Janet basis. The Janet basis is stored
       * in a additional vector of RingElems, to avoid the (more or less)
       * costly extraction of the Janet basis form the JBMill
       * @return a reference to a vector of RingElems
       */
      const std::vector<RingElem>& myJanetBasis() const;

      void setStrategy(Involutive::StrategyFlag strat) const
      {
        strategy = strat;
      }

    private:
      // gens
      const std::vector<RingElem> generators;

      mutable std::vector<RingElem> janetBasis;

      // mills
      mutable std::unique_ptr<JBMill> jbMillPtr;
      mutable std::unique_ptr<PBMill> pbMillPtr;
      mutable Involutive::StrategyFlag strategy;

      // values
      // computable with JB
      mutable std::map<PPMonoidElem, std::vector<bool> > multVars;
      mutable std::map<PPMonoidElem, std::vector<bool> > nonMultVars;
      mutable RingElem hilbertPol;
      mutable RingElem hilbertSeries;
      mutable std::unique_ptr<FGModule> firstSyzygyPtr;
      mutable long degree;
      // non-static data member initializers only available with -std=c++11 or -std=gnu++11
      // mutable long dimension = -1;
      mutable long dimension;
      // complementary decomposition of the leading ideal
      mutable std::unique_ptr<ComplementaryDecomposition> complementaryDecompositionPtr;
      mutable bool3 isHomogeneous;
      mutable bool3 isMonomial;
      mutable bool3 isDeltaRegular;

      // computable only with PB
      mutable long depth;
      mutable long projDim;
      mutable std::unique_ptr<std::vector<RingElem> > soclePtr;
      mutable std::map<std::pair<long, long>, long> extremalBettiNumbers;
      mutable std::vector<RingElem> regularSequence;
      mutable std::vector<RingElem> maximalStronglyIndependentSet;
      mutable long regularity;
      mutable long satiety;
      mutable std::vector<RingElem> saturation;
      mutable bool3 isCohenMacaulayRing;

      /**
       * This method computes the multiplicative and nonmultiplicative variables
       * at once.
       */
      void myMultAndNonMultVars() const;
    };

  } // end namespace Involutive
} // end namespace CoCoA

#endif
