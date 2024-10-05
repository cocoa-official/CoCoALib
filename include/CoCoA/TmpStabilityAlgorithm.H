#ifndef CoCoA_StabilityAlgorithm_H
#define CoCoA_StabilityAlgorithm_H

//   Copyright (c)  2014-2015  John Abbot,  Anna M. Bigatti
//   Original authors:  Pierre Pytlik, Mario Albert

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

#include "TmpJBMill.H"
#include <algorithm>

namespace CoCoA {
  namespace Involutive {
    class JBMill;


    /**
     * This class serves as base class for all classes which transform ideals
     * to various 'stability' positions.
     */
    class StabilityAlgorithm {
    public:
      enum StatisticLevel {None, Tracking, Logging};
      /**
       * This nested class defines some variables to track statistics.
       * At the moment we can tracking JBMills and the transformations represented
       * by a vector of RingElems
       */
      struct Statistics
      {
        std::vector<JBMill> mills;
        std::vector<std::vector<RingElem> > transformations;
      };

      virtual ~StabilityAlgorithm()
      {}

      /**
       * Purely virtual. This method shall compute the desired position for the given Input.
       * The desired position is given by the class where the method is reimplemented.
       * @param beginInput
       * @param endInput
       */
      virtual void myComputer(std::vector<RingElem>::const_iterator beginInput,
                              std::vector<RingElem>::const_iterator endInput) = 0;

      /**
       * Purely virtual. This method shall return the JanetContainer which contains the JB
       * with desired position
       *
       * @return JanetConatiner
       */
      virtual JanetContainer myOutputResult() = 0;

      /**
       * This method sets the statistic Level
       *
       * @param NewStatLevel
       */
      void SetStatisticLevel(StatisticLevel NewStatLevel)
      {
        StatLevel = NewStatLevel;
      }

      /**
       * This method sets the strategy which should be used to
       * compute a Janet basis
       *
       * @param strat
       */
      void SetJanetStrategy(StrategyFlag strat)
      {
        myJanetStrategy = strat;
      }

      /**
       * This method returns the collected statistics during the coordinate transformations.
       * stats is nonempty if the StatisticLevel was Tracking or Logging
       *
       * @return a pair where the first entry is a vector of JBMills and the second entry is a vector
       * of vectors of RingElems which represents the coordinate transformations.
       */
      const Statistics& myOutputStatistics()
      {
        return stats;
      }

    protected:
      /**
       * Constructor: initializes myPolyRing and myContainer with PolyRing
       *
       * @param PolyRing
       */
      StabilityAlgorithm(const SparsePolyRing& PolyRing)
        : myPolyRing(PolyRing),
          myContainer(PolyRing),
          StatLevel(None),
          myJanetStrategy(TQBlockLow),
          stats()
      {}

      /**
       * This method applies the ring homomorphism hom to the values between
       * beg and end
       * NOT EXCEPTION SAFE!
       *
       * @param beg
       * @param end
       * @param hom
       */
      virtual void myApplyHomToVector(std::vector<RingElem>::iterator beg, std::vector<RingElem>::iterator end, RingHom& hom) const;

      /**
       * Checks if StatLevel is not None. If this is the case
       * it pushes the given mill to stats.mills and if StatLevel = Logging it prints
       * the result to cout
       *
       * @param mill a const ref JBMill
       */
      inline void myTrackJBMill(const JBMill& mill)
      {
        if (StatLevel != None)
        {
// #ifdef MARIO11
          stats.mills.push_back(mill);
// #endif
          myPrintJBMill(mill);
        }
      }

      /**
       * Checks if StatLevel is not None. If this is the case
       * it pushes the given transformation to stats.transformations and
       * if StatLevel = Logging it prints the result to cout
       *
       * @param transformation a const ref vector of RingElems
       */
      inline void myTrackTransformation(const std::vector<RingElem>& transformation)
      {
        if (StatLevel != None)
        {
          stats.transformations.push_back(transformation);
          myPrintTransformation(transformation);
        }
      }

      const SparsePolyRing myPolyRing;
      JanetContainer myContainer;
      StatisticLevel StatLevel;
      StrategyFlag myJanetStrategy;

      Statistics stats;

    private:
      /**
       * Checks if StatLevel is Logging. If this is the case it prints
       * the Janet basis of mill to cout
       *
       * @param mill
       */
      void myPrintJBMill(const JBMill& mill) const;

      /**
       * Checks if StatLevel is Logging. If this is the case it prints
       * the transformation of mill to cout
       *
       * @param transformation
       */
      void myPrintTransformation(const std::vector<RingElem>& transformation) const;
    };
      // struct StabilityAlgorithm::Statistics
      // {
      //   std::vector<JBMill> mills;
      //   std::vector<std::vector<RingElem> > transformations;
      // };


    /**
     * This class and all its subclasses transform a given input to delta-regular coordinates.
     * This class always performs variable transformations of the form x_i + x_j
     */
    class DeltaRegular : public StabilityAlgorithm
    {
    public:
      /**
       * This enum represents the different strategies for the usage of variable permutations to
       * compute delta-reg. coordinates
       * UsePermuations
       * NotUsePermuations
       */
      enum UsagePermutations
      {
        UsePermutations, NotUsePermutations
      };

      /**
       * Constructor: Initializes the class with the default Permutations strategy
       *
       * @param PolyRing
       * @param use
       */
      DeltaRegular(const SparsePolyRing& PolyRing, UsagePermutations use = NotUsePermutations)
        : StabilityAlgorithm(PolyRing),
          myUsagePermutations(use)
      {}

      /**
       * This computes the delta-regular position for the given Input.
       * @param beginInput
       * @param endInput
       */
      void myComputer(std::vector<RingElem>::const_iterator beginInput,
                      std::vector<RingElem>::const_iterator endInput);

      /**
       * This method returns the JanetContainer which contains the JB
       * with delta-regular position
       *
       * @return JanetContainer
       */
      JanetContainer myOutputResult()
      {
        return myContainer;
      }

    protected:
      UsagePermutations myUsagePermutations;

      /**
       * Virtual: Computes the images which used for coordinate transformations.
       * This method only produces images of the form x_i + x_j
       *
       * @param mill
       * @return a vector of RingElems which represents the image
       */
      virtual std::vector<RingElem> DoComputeImage(const JBMill& mill) const;

    private:
      /**
       * Checks if r is of type x_i^r.
       * @param r
       * @return true if r is of type x_i^r
       */
      static bool myIsNotIndetPower(ConstRefRingElem r)
      {
        return UnivariateIndetIndex(r) == -1;
      }

      /**
       * r1, r2 must be of type x_i^r. It sorts first by degree and then by index
       *
       * @param r1
       * @param r2
       *
       * @return returns true if deg(r1) < deg(r2) or (deg(r1) = deg(r2) and index(r1) < index(r2))
       */
      static bool mySortByPowerAndIndex(ConstRefRingElem r1, ConstRefRingElem r2);

      /**
       * Permutes the variables in the following manner:
       * At first there are the univariate indets (sorted by function mySortByPowerAndIndex)
       * then it appends the rest of the variables
       *
       * Example:
       * generators = {x0 + x3 +x5, x3^5, x2, x1^2, x4^2}
       * result = {x2, x1, x4, x3, x0, x5}
       *
       * @param generators
       * @return a vector of RingElems, which contains the permuted variables.
       */
      std::vector<RingElem> ReturnPermutation(const std::vector<RingElem>& generators) const;

      /**
       * This method applies the homomorphism which is represented by images
       * on the range specified by beginGens and endGens.
       *
       * @param images
       * @param beginGens
       * @param endGens
       */
      void myPerformTransformation(const std::vector<RingElem>& images,
                                   std::vector<RingElem>::iterator beginGens,
                                   std::vector<RingElem>::iterator endGens) const;
    };

    /**
     * This class transforms a given input to delta-regular coordinates.
     * Its base class is DeltaRegular. It uses nearly all methods from the class DeltaRegular.
     * The only difference occurs in the method DoComputeImage. The class DeltaRegular only computes images
     * of the form x_i -> x_i + x_j. The implementation of DoComputeImage supplies images of the form
     * x_i -> x_i + x_j_1 + x_j_2 + ...
     */
    class DeltaRegularAll : public DeltaRegular
    {
    public:
      /**
       * Constructor: Initializes the class with the default Permutations strategy (see class DeltaRegular)
       *
       * @param PolyRing
       * @param use
       */
      DeltaRegularAll(const SparsePolyRing& PolyRing,UsagePermutations use = UsePermutations)
        : DeltaRegular(PolyRing, use)
      {}

    private:
      /**
       * Computes the images which used for coordinate transformations.
       * This method produces images of the form x_i -> x_i + x_j_1 + x_j_2 + ...
       *
       * @param mill
       * @return a vector of RingElems which represents the image
       */
      std::vector<RingElem> DoComputeImage(const JBMill& mill) const;
    };

    /**
     * This class and all its subclasses transform a given input at least to a stable position (resp. the
     * leading ideal is stable)
     * This class always performs variable transformations of the form x_i + x_j
     */
    class StableLTI : public StabilityAlgorithm {
    public:
      /**
       * Constructor: Initializes the class
       *
       * @param PolyRing
       */
      StableLTI(const SparsePolyRing& PolyRing) : StabilityAlgorithm(PolyRing)
      {}

      /**
       * This computes the stable position for the given Input.
       * @param beginInput
       * @param endInput
       */
      void myComputer(std::vector<RingElem>::const_iterator beginInput,
                      std::vector<RingElem>::const_iterator endInput);

      /**
       * This method returns the JanetContainer which contains the JB
       * with stable position
       *
       * @return JanetConatiner
       */
      JanetContainer myOutputResult()
      {
        return myContainer;
      }
    protected:

      /**
       * Computes the class of LPP(r). This is
       * max{i | deg_i(LPP(r)) != 0}
       *
       * @param r
       * @return the class as long
       */
      long myCls(ConstRefRingElem r) const;

      /**
       * Virtual: Computes the images which used for coordinate transformations.
       * This method only produces images of the form x_i + x_j
       *
       * @param mill
       * @return a vector of RingElems which represents the image
       */
      virtual std::vector<RingElem> DoComputeImage(const std::vector<RingElem>& GensLTIdeal, const ideal& LTIdeal) const;

      /**
       * Virtual: Check if the coordinates of generators have the correct position.
       *
       * @param generators
       *
       * @return returns true, if the input is in correct position.
       */
      virtual bool IamInPosition(const std::vector<RingElem>& generators) const
      {
        return myIsLTIStable(generators);
      }

      /**
       * Extracts the leading terms of all elements in v as RingElems.
       *
       * @param v
       * @return a vector of RingElems, which are the LPP's of v
       */
      std::vector<RingElem> myLeadingMonomials(const std::vector<RingElem>& v) const;
    private:

      /**
       * Checks if generators is stable that is:
       * Let v = LPP(generators).
       * generators is stable iff for all x^\mu \in v
       * (x^\mu / x_{\cls(x^\mu)}) * x_i \in <v>. (i < \cls(x^\mu))
       *
       * @param generators
       * @return returns true if generators is in stable position.
       */
      bool myIsLTIStable(const std::vector<RingElem>& generators) const;
    };

    /**
     * This class transforms a given input to a stable position (resp. the
     * leading ideal is stable). Its base class is StableLTI. It uses nearly
     * all methods from the class StableLTI. The only difference occurs in the
     * method DoComputeImage. The class StableLTI only computes images
     * of the form x_i -> x_i + x_j. The implementation of DoComputeImage
     * supplies images of the form x_i -> x_i + x_j_1 + x_j_2 + ...
     */
    class StableLTIAll : public StableLTI {
    public:
      /**
       * Constructor: Initializes the class
       *
       * @param PolyRing
       */
      StableLTIAll(const SparsePolyRing& PolyRing) : StableLTI(PolyRing)
      {}

    private:
      /**
       * Computes the images which used for coordinate transformations.
       * This method produces images of the form x_i -> x_i + x_j_1 + x_j_2 + ...
       *
       * @param mill
       * @return a vector of RingElems which represents the image
       */
      std::vector<RingElem> DoComputeImage(const std::vector<RingElem>& GensLTIdeal, const ideal& LTIdeal) const;
    };

    /**
     * This class and all its subclasses transform a given input at least to a
     * strongly stable position (resp. the leading ideal is strongly stable)
     * This class always performs variable transformations of the form x_i + x_j.
     * It is a subclass of StableLTI, because the positions strongly stable and stable
     * and the computation of this positions are similar.
     */
    class StronglyStableLTI : public StableLTI {
    public:
      /**
       * Constructor: Initializes the class
       *
       * @param PolyRing
       */
      StronglyStableLTI(const SparsePolyRing& PolyRing) : StableLTI(PolyRing)
      {}

    protected:
      /**
       * Virtual: Computes the images which used for coordinate transformations.
       * This method only produces images of the form x_i + x_j
       *
       * @param mill
       * @return a vector of RingElems which represents the image
       */
      virtual std::vector<RingElem> DoComputeImage(const std::vector<RingElem>& GensLTIdeal, const ideal& LTIdeal) const;

    private:
      /**
       * Virtual: Check if the coordinates of generators have the correct position.
       *
       * @param generators
       *
       * @return returns true, if the input is in correct position.
       */
      bool IamInPosition(const std::vector<RingElem>& generators) const
      {
        return myIsLTIStronglyStable(generators);
      }

      /**
       * Checks if generators is strongly stable that is:
       * Let v = LPP(generators).
       * generators is strongly stable iff for all x^\mu \in v
       * (x^\mu / x_j}) * x_i \in <v>.
       * (i < j and j <= \cls(x^\mu) and x_j divides x^\mu)
       *
       * @param generators
       * @return returns true if generators is in stable position.
       */
      bool myIsLTIStronglyStable(const std::vector<RingElem>& generators) const;
    };

    /**
     * This class transforms a given input to a strongly stable position (resp. the
     * leading ideal is strongly stable). Its base class is StronglyStableLTI. It uses nearly
     * all methods from the class StableLTI and StronglyStableLTI. The only difference occurs in the
     * method DoComputeImage. The class StronglyStableLTI only computes images
     * of the form x_i -> x_i + x_j. The implementation of DoComputeImage
     * supplies images of the form x_i -> x_i + x_j_1 + x_j_2 + ...
     */
    class StronglyStableLTIAll : public StronglyStableLTI {
    public:
      /**
       * Constructor: Initializes the class
       *
       * @param PolyRing
       */
      StronglyStableLTIAll(const SparsePolyRing& PolyRing) : StronglyStableLTI(PolyRing)
      {}

    private:
      /**
       * Computes the images which used for coordinate transformations.
       * This method produces images of the form x_i -> x_i + x_j_1 + x_j_2 + ...
       *
       * @param mill
       * @return a vector of RingElems which represents the image
       */
      std::vector<RingElem> DoComputeImage(const std::vector<RingElem>& GensLTIdeal, const ideal& LTIdeal) const;

    };
  } // end namespace Involutive
} // end namespace CoCoA

#endif
