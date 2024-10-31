#ifndef CoCoA_PBMill_H
#define CoCoA_PBMill_H

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

#include "CoCoA/TmpJBMill.H"
#include "CoCoA/TmpStabilityAlgorithm.H"

namespace CoCoA
{
  namespace Involutive
  {

    /**
     * PBMill is a subclass of JBMill. It models the mathematical relationship between
     * Janet and Pommaret bases (A Pommaret basis is always a Janet basis). It implements
     * some methods to compute several invariants of an ideal from a Pommaret basis.
     *
     * In addition the class contains some nested classes which convert/transform a given
     * ideal to a delta-regular ideal which is represented by the Pommaret basis.
     *
     */
    class PBMill : public JBMill
    {
    public:
      class Builder;
      class Converter;
      class Transformator;
      class DeltaRegularTransformator;
      class StableLTITransformator;
      class StronglyStableLTITransformator;

      /**
       * Constructor. It creates a PBMill and checks if the
       * JanetContainer contains a Pommaret basis. If this is not
       * the case it raises an error.
       *
       * @param basis is a JanetContainer
       */
      PBMill(const JanetContainer& basis) : JBMill(basis)
      {
        if (!IamPommaretBasis())
        {
          CoCoA_THROW_ERROR("This is not a Pommaret Basis", "PBMill");
        }
      }

      /**
       * Constructs a PBMill out of a PBMill::Builder object
       *
       * @param builder
       */
      PBMill(const PBMill::Builder& builder);

      /**
       * Returns the minimal Pommaret basis.
       * This method is only an alias for myReturnJB().
       *
       * @return a vector of RingElems
       */

      inline std::vector<RingElem> myReturnPB() const
      {
       return myReturnJB();
      }

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
      std::vector<RingElem> myRegSeq() const;

      /**
       * Computes the unique maximal strongly independent set of P/I
       *
       * @return the strongly independent set as vector of RingElems
       */
      std::vector<RingElem> myMaxStronglyIndependentSet() const;

      /**
       * Tests if I is Cohen Macualay
       *
       * @return true if I is Cohen Macaulay
       */
      bool IamCohenMacaulay() const;

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
       * Computes the saturation of I. The output is not minimized.
       *
       * @return the saturation as vector of RingElems
       */
      std::vector<RingElem> mySaturation() const;

    private:
      /**
       * Counts every element in vec which has degree degree
       *
       * @param vec const ref to vector of RingElems
       * @param degree long
       *
       * @return Number of elements in vec with degree degree as long
       */
      long myCountElementsWithDegree(const std::vector<RingElem>& vec, long degree) const;

      /**
       * This method takes the Pommaret basis and returns it divided by classes
       *
       * @return a multimap where the key is the class and the value is a RingElem
       * which corresponds to the class
       */
      std::multimap<long, RingElem> myBasisSeperatedByNumberOfMultVar() const;

      /*
       * the maximal degree of all elements with
       * class index in the pommaret basis
       */

      /**
       * This method returns the maximal degree of an element of the PB
       * with index multiplicative variables
       *
       * @param index
       * @return the degree as long
       */
      long myDegForNumMultVar(long index) const;
      /*
       * a map with all classes in the pommaret basis +
       * associated degree (see above)
       */
      /**
       * This method computes the maximal degree of elements in the PB
       * divided by number of mult vars.
       *
       * @return a map where the first entry represents the number of mult vars and the second
       * the maximal degree which occurs with this number of mult vars.
       */
      std::map<long, long> myDegSeperatedByNumMultVar() const;
    };

    /**
     * This nested class is the base class for all builder-classes which
     * build a Pommaret basis. This class is purely virtual.
     * The 'Gang of Four' (Gamma, Helm, Johnson, Vlissides - Design Patterns) served
     * as template for the construction. The corresponding Pattern is called Building Pattern
     */
    class PBMill::Builder{
      friend PBMill::PBMill(const PBMill::Builder& builder);
    protected:
      /**
       * This is a purely virtual method, which should build a JanetContainer which is in
       * correct position
       *
       * @return JanetConatiner
       */
      virtual JanetContainer build() const = 0;



      /**
       * Protected Constructor
       */
      Builder()
      {}

    };

    /**
     * This is the converter-builder. It takes a JBMill
     * and converts it to a PBMill.
     *
     */
    class PBMill::Converter : public PBMill::Builder{
    public:
      Converter()
      {}

      /**
       * This method sets the JBMill which should be converted.
       * It raises an error if the JBMill is not delta-regular
       * (i.e. the Janet basis is not a Pommaret basis)
       *
       * @param mill a JBMill
       * @return a reference to the current Builder object
       */
      Builder& setJBMill(const JBMill& mill)
      {
        if (!mill.IamPommaretBasis())
        {
          CoCoA_THROW_ERROR("Janet basis is not a Pommaret basis", "PBMill::Converter::setJBMill");
        }
        myContainerPtr = std::unique_ptr<JanetContainer> (new JanetContainer(mill.myGetJanetContainer()));
        return *this;
      }

    private:
      /**
       * Builds a JanetContainer
       *
       * @return JanetContainer
       */
      JanetContainer build() const  // wasteful copy?
      {
        if (myContainerPtr.get() == nullptr)
        {
          CoCoA_THROW_ERROR("There must be a Janet basis as Input", "PBMill::Converter::build");
        }
        return *myContainerPtr;
      }

      std::unique_ptr<JanetContainer> myContainerPtr;

    };

    /**
     * This class is a subclass of PBMill::Builder. It is designed as a virtual class,
     * e.g. the constructor is protected and the user can't call this class directly.
     * It acts as base class for all Builders which perform coordinate-transformations
     * to get a Pommaret basis.
     */
    class PBMill::Transformator : public PBMill::Builder{
    public:

      /**
       * Builds a PBMill and returns all coordinate transformations which were done
       *
       * @return a pair: The first entry is the PBMill, the second entry is again a pair
       * where the first entry is a vector of JBMills, which occur during the coordinate transformations
       * and the second entry is a vector of vector of RingElems which represent the coordinate transformations.
       */
      std::pair<PBMill, StabilityAlgorithm::Statistics> BuildWithStatistics()
      {
        std::pair<JanetContainer, StabilityAlgorithm::Statistics> result(myUseStrategyWithStatistics());
        return std::make_pair(PBMill(result.first), result.second);
      }

      /**
       * Set the input to the given range. The method raises
       * an error if the input is empty.
       *
       * @param inputBegin
       * @param inputEnd
       *
       * @return a reference to the current Builder object
       */

      Transformator& setInput(std::vector<RingElem>::const_iterator inputBegin, std::vector<RingElem>::const_iterator inputEnd)
      {
        myInput.clear();
        myInput.insert(myInput.begin(), inputBegin, inputEnd);
        if (myInput.empty())
        {
          CoCoA_THROW_ERROR("empty input", "setInput");
        }
        return *this;
      }

      Transformator& setJanetStrategy(StrategyFlag JanetStrat)
      {
        myJanetStrategy = JanetStrat;
        return *this;
      }

      /**
       * Set the input to input. The method raises an error
       * if the input is empty.
       *
       * @param input a vector of RingElems
       *
       * @return a reference to the current Builder object
       */
      Transformator& setInput(const std::vector<RingElem>& input)
      {
        return setInput(input.begin(), input.end());
      }

      Transformator& setStatisticLevel(StabilityAlgorithm::StatisticLevel NewStatLevel)
      {
        StatLevel = NewStatLevel;
        return *this;
      }

    protected:
      /**
       * Builds a JanetConatiner
       *
       * @return JanetConatiner
       */
      JanetContainer build() const // cannot return ref since myUseStrategy returns temp obj
      {
        if (myInput.empty())
        {
          CoCoA_THROW_ERROR("empty input", "build");
        }
        return myUseStrategy();
      }

      std::vector<RingElem> myInput;

      StabilityAlgorithm::StatisticLevel StatLevel;

      StrategyFlag myJanetStrategy;

      /**
       * Protected constructor
       */
      Transformator() : PBMill::Builder(),
                        StatLevel(StabilityAlgorithm::None),
                        myJanetStrategy(TQBlockLow)
      {}

      /**
       * Extracts the current poly ring
       *
       * @return SparsePolyRing
       */
      SparsePolyRing myPolyRing() const
      {
        if (myInput.empty())
        {
          CoCoA_THROW_ERROR(ERR::ReqNonEmpty, __FUNCTION__);
        }
        return owner(myInput.front());
      }

      /**
       * This method uses the specified strategy to compute a JanetContainer.
       * The build method uses this JanetContainer to build the PBMill.
       *
       * @return a JanetContainer
       */
      JanetContainer myUseStrategy() const;

      /**
       * This method uses the specified strategy to compute a JanetContainer.
       * The build method uses this JanetContainer to build the PBMill.
       * In addition to the method myUseStrategy it returns some statistics.
       *
       * @return a pair: The first entry is the JanetContainer, the second entry is again a pair
       * where the first entry is a vector of JBMills, which occur during the coordinate transformations
       * and the second entry is a vector of vector of RingElems which represent the coordinate transformations.
       */

      std::pair<JanetContainer, StabilityAlgorithm::Statistics> myUseStrategyWithStatistics();

      /**
       * Purely virtual method. The task of this method is to return the chosen strategy as unique_ptr.
       * @return an unique_ptr to StabilityAlgorithm
       */
      virtual std::unique_ptr<StabilityAlgorithm> myChooseStrategy() const = 0;
    };

    /**
     * This class transforms a given input to delta-regular coordinates. It is a subclass
     * of PBMill::Transformator
     */
    class PBMill::DeltaRegularTransformator : public PBMill::Transformator{
    public:
      /**
       * This enum represents the different strategies:
       * SingleWithPermutation: Uses class DeltaRegular with flag DeltaRegular::UsePermutations
       * SingleWithoutPermutation: Uses class DeltaRegular with flag DeltaRegular::NotUsePermutations
       * AllWithPermutation: Uses class DeltaRegularAll with flag DeltaRegular::UsePermutations
       * AllWithoutPermutation: Uses class DeltaRegularAll with flag DeltaRegular::NotUsePermutations
       */
      enum StrategyFlag
      {
        SingleWithPermutation,
        SingleWithoutPermutation,
        AllWithPermutation,
        AllWithoutPermutation
      };

      static const StrategyFlag DEFAULT_STRATEGY = SingleWithoutPermutation;

      /**
       * Constructor. Initialises class with default strategy
       */
      DeltaRegularTransformator() : PBMill::Transformator(), myStrategy(DEFAULT_STRATEGY)
      {}

      /**
       * Set the strategy to compute a delta-regular coordinates
       *
       * @param strategy
       *
       * @return a reference to the current Builder object
       */
      DeltaRegularTransformator& setStrategy(StrategyFlag strategy)
      {
        myStrategy = strategy;
        return *this;
      }

    private:
      StrategyFlag myStrategy;

      /**
       * The task of this method is to return the chosen strategy as unique_ptr.
       * @return an unique_ptr to StabilityAlgorithm
       */
      std::unique_ptr<StabilityAlgorithm> myChooseStrategy() const;
    };

    /**
     * This class transforms a given input to a stable position, e.g. the leading ideal
     * will be in stable position. It is a subclass of PBMill::Transformator
     */
    class PBMill::StableLTITransformator : public PBMill::Transformator{
    public:
      /**
       * This enum represents the different strategies:
       * Single: Uses class StableLTI
       * All: Uses class StableLTIAll
       */
      enum StrategyFlag
      {
        Single,
        All
      };

      static const StrategyFlag DEFAULT_STRATEGY = Single;

      /**
       * Constructor. Initialises class with default strategy
       */
      StableLTITransformator() : PBMill::Transformator(), myStrategy(DEFAULT_STRATEGY)
      {}

      /**
       * Set the strategy to compute a stable position.
       *
       * @param strategy
       *
       * @return a reference to the current Builder object
       */
      StableLTITransformator& setStrategy(StrategyFlag strategy)
      {
        myStrategy = strategy;
        return *this;
      }


    private:
      StrategyFlag myStrategy;

      /**
       * The task of this method is to return the chosen strategy as unique_ptr.
       * @return an unique_ptr to StabilityAlgorithm
       */
      std::unique_ptr<StabilityAlgorithm> myChooseStrategy() const;
    };

    /**
     * This class transforms a given input to a strongly stable position, e.g. the leading ideal
     * will be in stable position. It is a subclass of PBMill::Transformator
     */
    class PBMill::StronglyStableLTITransformator : public PBMill::Transformator{
    public:
      /**
       * This enum represents the different strategies:
       * Single: Uses class StronglyStableLTI
       * All: Uses class StronglyStableLTIAll
       */
      enum StrategyFlag
      {
        Single,
        All
      };

      static const StrategyFlag DEFAULT_STRATEGY = Single;

      /**
       * Constructor. Initialises class with default strategy
       */
      StronglyStableLTITransformator() : PBMill::Transformator(), myStrategy(DEFAULT_STRATEGY)
      {}

      /**
       * Set the strategy to compute a strongly stable position.
       *
       * @param strategy
       *
       * @return a reference to the current Builder object
       */
      StronglyStableLTITransformator& setStrategy(StrategyFlag strategy)
      {
        myStrategy = strategy;
        return *this;
      }


    private:
      StrategyFlag myStrategy;

      /**
       * The task of this method is to return the chosen strategy as unique_ptr.
       * @return an unique_ptr to StabilityAlgorithm
       */
      std::unique_ptr<StabilityAlgorithm> myChooseStrategy() const;
    };
  }
}
#endif

