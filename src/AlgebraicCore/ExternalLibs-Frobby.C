//   Copyright (c)  2008  John Abbott,  Anna M. Bigatti
//   Original author: 2008  Bjarke Hammersholt Roune

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


// Source code for Frobby integration

#include "CoCoA/ExternalLibs-Frobby.H"

#ifdef CoCoA_WITH_FROBBY

#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/convert.H"
#include "CoCoA/ideal.H"
#include "CoCoA/ring.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H" // for len

#include "gmp.h"
#include "frobby.h"

#include <string>
using std::string;
#include <vector>
using std::vector;


namespace CoCoA
{

  const char* FrbVersion()
  {
    return "UNKNOWN (0.9.4?)";
//    return constants::version;
  }

  namespace // anonymous
  {
	PPMonoidElem ExponentsToMonomial(mpz_ptr* exponentVector,
									 const PPMonoid& monoid)
	{
	  // ??? use this in the monomial ideal consumer.
	  const std::size_t numIndets = NumIndets(monoid);
	  PPMonoidElem m = one(monoid);
	  for (size_t indet = 0; indet < numIndets; ++indet)
            m *= IndetPower(monoid, indet, BigIntFromMPZ(exponentVector[indet]));
	  return m;
	}

    // Stores the consumed monomials into the passed-in ideal.
    class FrobbyMonomialIdealConsumer : public ::Frobby::IdealConsumer
    {
    public:
      FrobbyMonomialIdealConsumer(const SparsePolyRing& ring): myRing(ring) {}

      void idealBegin(size_t /*varCount*/) override
      {
        myIdeals.push_back(ideal(myRing, vector<RingElem>()));
      }

      void consume(mpz_ptr* exponentVector) override
      {
        const std::size_t count = NumIndets(myRing);

        const std::vector<RingElem>& x = indets(myRing);
        RingElem m(myRing);
        if (count > 0)
          m = power(x[0], 0);
        for (size_t indet = 0; indet < count; ++indet)
          m = m * power(x[indet], BigIntFromMPZ(exponentVector[indet]));

        myIdeals.back() += ideal(m);
      }

      std::vector<ideal>& getIdealsRef() { return myIdeals; }

    private:
      SparsePolyRing myRing;
      std::vector<ideal> myIdeals;
    };

	class FrobbyPolynomialConsumer : public ::Frobby::PolynomialConsumer
	{
	public:
	  FrobbyPolynomialConsumer(const SparsePolyRing& ring):
		myRing(ring),
		myPoly(ring) {}

	  void polynomialBegin(size_t /*varCount*/) override
	  {
		myPoly = zero(myRing);
	  }

	  void consume(const mpz_t coefficient, mpz_ptr* exponentVector) override
	  {
		const PPMonoid& ppm = PPM(myRing);
		myPoly += monomial(myRing,
                                   BigIntFromMPZ(coefficient),
                                   ExponentsToMonomial(exponentVector, ppm));
	  }

	  RingElem& getPolyRef() {return myPoly;}

	private:
	  SparsePolyRing myRing;
	  RingElem myPoly;
	};

	class FrobbyEvalPolynomialConsumer : public ::Frobby::PolynomialConsumer
	{
	public:
	  FrobbyEvalPolynomialConsumer(const RingElem& base):
		myBase(base),
		myPoly(owner(base)) {}

	  void polynomialBegin(size_t /*varCount*/) override
	  {
		myPoly = zero(owner(myBase));
	  }

	  void consume(const mpz_t coefficient, mpz_ptr* exponentVector) override
	  {
            myPoly += BigIntFromMPZ(coefficient) *
              power(myBase, BigIntFromMPZ(exponentVector[0]));
	  }

	  RingElem& getPolyRef() {return myPoly;}

	private:
	  RingElem myBase;
	  RingElem myPoly;
	};

    // The purpose of this class is to encapsulate initialization and
    // deallocation of mpz_t, while guaranteeing that the only member
    // variable is an mpz_t, so that an array of WrappedMpzT can be
    // safely cast to an array of mpz_t. mpz_class does not guarantee
    // this. In the absense of this class, it becomes very
    // inconvenient to allocate and initialize an array of mpz_t's in
    // an exception-safe way.
    class WrappedMpzT
    {
    public:
      WrappedMpzT() { mpz_init(_mpz); }

      ~WrappedMpzT() { mpz_clear(_mpz); }

      mpz_t& getMpzT() { return _mpz; }

    private:
      mpz_t _mpz;
    };

    void ToFrobbyIdeal(Frobby::Ideal& frobbyIdeal, const ideal& I)
    {
      const std::size_t nvars = NumIndets(RingOf(I));
      for (const RingElem& g: gens(I))
      {
        if (IsZero(g)) continue;
        for (size_t indet = 0; indet < nvars; ++indet)
            frobbyIdeal.addExponent(mpzref(BigExponent(LPP(g), indet)));
      }
    }

    void MustHaveMonomialGens(const ideal& I,
                              const string& operation) {
      if (!AreGensMonomial(I))
      {
        string msg = operation + " (Frobby interface)";
        for (long i=0; i<len(gens(I)); ++i) // zeroes are acceptable
          if ((!IsZero(gens(I)[i]) && !IsMonomial(gens(I)[i])))
            CoCoA_THROW_ERROR(ERR::ReqMonomialGens, msg);
      }
    }
  } // end of namespace anonymous


  ideal FrbAlexanderDual(const ideal& I, ConstRefPPMonoidElem pp)
  {
    MustHaveMonomialGens(I, "FrbAlexanderDual");
    const SparsePolyRing polyRing = RingOf(I);
    if (PPM(polyRing) != owner(pp))
      CoCoA_THROW_ERROR(ERR::MixedPPMs, "FrbAlexanderDual");

    const std::size_t count = NumIndets(polyRing);

    Frobby::Ideal frobbyIdeal(count);
    ToFrobbyIdeal(frobbyIdeal, I);

    FrobbyMonomialIdealConsumer consumer(polyRing);

    // Set up the point to dualize on. 
    WrappedMpzT* exponentVector = new WrappedMpzT[count];
    try
    {
      for (size_t indep = 0; indep < count; ++indep)
        mpz_set(exponentVector[indep].getMpzT(), mpzref(BigExponent(pp, indep)));

      // Compute Alexander dual using Frobby.
      Frobby::alexanderDual(frobbyIdeal, (mpz_t*)exponentVector, consumer);
    }
    catch (...)
    {
      delete[] exponentVector;
      throw;
    }

    delete[] exponentVector;

    return consumer.getIdealsRef().back();
  }

  ideal FrbAlexanderDual(const ideal& I)
  {
	MustHaveMonomialGens(I, "AlexanderDual");

    const SparsePolyRing polyRing = RingOf(I);
    Frobby::Ideal frobbyIdeal(NumIndets(polyRing));
    ToFrobbyIdeal(frobbyIdeal, I);

    FrobbyMonomialIdealConsumer consumer(polyRing);
    Frobby::alexanderDual(frobbyIdeal, (mpz_t*)nullptr, consumer);
    return consumer.getIdealsRef().back();
  }

  void FrbIrreducibleDecomposition(std::vector<ideal>& components, const ideal& I)
  {
	MustHaveMonomialGens(I, "IrreducibleDecomposition");

    const SparsePolyRing polyRing = RingOf(I);
    Frobby::Ideal frobbyIdeal(NumIndets(polyRing));
    ToFrobbyIdeal(frobbyIdeal, I);

    FrobbyMonomialIdealConsumer consumer(polyRing);
    Frobby::irreducibleDecompositionAsIdeals(frobbyIdeal, consumer);
    consumer.getIdealsRef().swap(components);
  }

  void FrbPrimaryDecomposition(std::vector<ideal>& components, const ideal& I)
  {
	MustHaveMonomialGens(I, "PrimaryDecomposition");

    const SparsePolyRing polyRing = RingOf(I);
    Frobby::Ideal frobbyIdeal(NumIndets(polyRing));
    ToFrobbyIdeal(frobbyIdeal, I);

    FrobbyMonomialIdealConsumer consumer(polyRing);
    Frobby::primaryDecomposition(frobbyIdeal, consumer);
    consumer.getIdealsRef().swap(components);
  }

  ideal FrbMaximalStandardMonomials(const ideal& I)
  {
	MustHaveMonomialGens(I, "MaximalStandardMonomials");

    const SparsePolyRing polyRing = RingOf(I);
    Frobby::Ideal frobbyIdeal(NumIndets(polyRing));
    ToFrobbyIdeal(frobbyIdeal, I);

    FrobbyMonomialIdealConsumer consumer(polyRing);
    Frobby::maximalStandardMonomials(frobbyIdeal, consumer);
    return consumer.getIdealsRef().back();
  }

  long FrbDimension(const ideal& I) {
	MustHaveMonomialGens(I, "Dimension");
	
    const SparsePolyRing polyRing = RingOf(I);
    Frobby::Ideal frobbyIdeal(NumIndets(polyRing));
    ToFrobbyIdeal(frobbyIdeal, I);

    BigInt dim;
    Frobby::dimension(frobbyIdeal, mpzref(dim));
    return ConvertTo<long>(dim);
  }

  void FrbAssociatedPrimes(std::vector<ideal>& primes, const ideal& I)
  {
	MustHaveMonomialGens(I, "AssociatedPrimes");

    const SparsePolyRing polyRing = RingOf(I);
    Frobby::Ideal frobbyIdeal(NumIndets(polyRing));
    ToFrobbyIdeal(frobbyIdeal, I);

    FrobbyMonomialIdealConsumer consumer(polyRing);
    Frobby::associatedPrimes(frobbyIdeal, consumer);
    consumer.getIdealsRef().swap(primes);
  }

  RingElem FrbMultigradedHilbertPoincareNumerator(const ideal& I) {
	MustHaveMonomialGens(I, "AssociatedPrimes");

    const SparsePolyRing polyRing = RingOf(I);
    Frobby::Ideal frobbyIdeal(NumIndets(polyRing));
    ToFrobbyIdeal(frobbyIdeal, I);

    FrobbyPolynomialConsumer consumer(polyRing);
    Frobby::multigradedHilbertPoincareSeries(frobbyIdeal, consumer);
    return consumer.getPolyRef();
  }

  RingElem FrbTotalDegreeHilbertPoincareNumerator(const ideal& I) {
	SparsePolyRing ring = NewPolyRing(RingZZ(), symbols("t"));
    return FrbTotalDegreeHilbertPoincareNumerator(I, indet(ring, 0));
  }

  RingElem FrbTotalDegreeHilbertPoincareNumerator(const ideal& I,
											   const RingElem& base) {
	MustHaveMonomialGens(I, "AssociatedPrimes");

    const SparsePolyRing polyRing = RingOf(I);
    Frobby::Ideal frobbyIdeal(NumIndets(polyRing));
    ToFrobbyIdeal(frobbyIdeal, I);

    FrobbyEvalPolynomialConsumer consumer(base);
	Frobby::univariateHilbertPoincareSeries(frobbyIdeal, consumer);
    return consumer.getPolyRef();	
  }
}

#endif
