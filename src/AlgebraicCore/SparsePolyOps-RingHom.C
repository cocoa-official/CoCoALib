//   Copyright (c)  2018  John Abbott and Anna M. Bigatti
//   Authors:  2005-2018  John Abbott and Anna M. Bigatti

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


// Source code for RingHom on SparsePolyRing

#include "CoCoA/SparsePolyRing.H"


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/DenseMatrix.H" // for MultiplicationMat/myDiv
#include "CoCoA/FGModule.H"  // for myGcd
#include "CoCoA/MatrixOps.H" // for LinSolve
#include "CoCoA/MatrixView.H" // for ZeroMat
#include "CoCoA/OpenMath.H"
#include "CoCoA/PPMonoidHom.H"
#include "CoCoA/QuotientRing.H" // for IsQuotientRing
#include "CoCoA/ReductionCog.H"
#include "CoCoA/RingDistrMPolyClean.H" // for NewPolyRing_DMP
#include "CoCoA/RingDistrMPolyInlFpPP.H" // for NewPolyRing_DMPII
#include "CoCoA/RingDistrMPolyInlPP.H" // for NewPolyRing_DMPI
#include "CoCoA/RingFp.H" // for IsRingFp
#include "CoCoA/RingQQ.H" // for IsQQ
#include "CoCoA/RingTwinFloat.H" // for IsRingTwinFloat
#include "CoCoA/RingZZ.H" // for IsZZ
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/TmpGOperations.H"  // for myIntersect, my Elim..
#include "CoCoA/TmpUniversalInvolutiveBasisContainer.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"  // for myGcd
#include "CoCoA/geobucket.H" // for myMul
#include "CoCoA/ideal.H"     // for myGcd
#include "CoCoA/matrix.H" // for OrdMat, myDivMod
#include "CoCoA/module.H"    // for myGcd
#include "CoCoA/random.H" // for RandomLongStream
#include "CoCoA/submodule.H"  // for myGcd
#include "CoCoA/symbol.H"
#include "CoCoA/verbose.H"

#include <algorithm>
using std::max;     // for MaxExponent, StdDeg
using std::remove;  // for myColon
using std::sort;    // for AreGoodIndetNames, QuotientBasisSorted
#include <functional>
using std::not1;    // for AreLPPSqFree
using std::ptr_fun; // for AreLPPSqFree
#include <iostream>
// using std::ostream in SparsePolyRingBase::myOutput
#include <iterator>
using std::back_inserter;
#include <list>
#include <map>
using std::map;
// using std::list;
#include <utility>
using std::make_pair;
using std::pair;
//#include <vector>
using std::vector;

namespace CoCoA
{


  //-- HomImpl ----------------------------------------

  SparsePolyRingBase::HomImpl::HomImpl(const SparsePolyRing& domain, const ring& codomain, const RingHom& CoeffHom, const vector<RingElem>& IndetImages):
      RingHomBase(domain, codomain),
      myCoeffHom(CoeffHom),
      myIndetImages(IndetImages)
  {
    // No need to check anything: checks already made when CoeffHom was built.
  }

namespace
{
  // ??? appropriate use of inheritance here?  this is getting pretty ugly

  // SPECIAL CASE: when all PPs map to just PPs **AND** coeffs map to coeffs!!
  void ApplySPRCodomain_(RingElem& image, ConstRefRingElem arg, const RingHom CoeffHom, const vector<RingElem>& IndetImages)
  {
    // NOTE: CoeffHom maps into S (and **not** into CoeffRing(S))
    const SparsePolyRing S = owner(image);
    const PPMonoid PPMS = PPM(S);
    geobucket gbk(S);

    const long NumInd = len(IndetImages);
    vector<PPMonoidElem> IndetImage; IndetImage.reserve(NumInd); for (int i=0; i < NumInd; ++i) IndetImage.push_back(LPP(IndetImages[i]));
    vector<long> expv(NumInd); // workspace for exponents in loop below
    for (SparsePolyIter it=BeginIter(arg); !IsEnded(it); ++it)
    {
      const RingElem CoeffImageAsPoly = CoeffHom(coeff(it));
      if (IsZero(CoeffImageAsPoly))  continue; // efficiency hack????
      const RingElem CoeffImage = LC(CoeffImageAsPoly);
      PPMonoidElem SummandPP = one(PPMS);
      exponents(expv, PP(it));
      for (long ind=0; ind < NumInd; ++ind)
      {
        if (expv[ind] == 0) continue;
        if (/*IsPPMonoidOv(myCodomain) && */IsMatrixOrdering(ordering(PPM(S))))
          if (expv[ind] >= 32749)  // BUG BUG: nasty hack to avoid exp overflow!!!
            CoCoA_THROW_ERROR1(ERR::ExpTooBig);
        SummandPP *= power(IndetImage[ind], expv[ind]);
      }
      RingElem ImageTerm = monomial(S, CoeffImage, SummandPP);
      gbk.myAddClear(ImageTerm, 1);
    }
    AddClear(image, gbk);
  }

  // assume image==0
  void ApplySPRCodomain(RingElem& image, ConstRefRingElem arg, const RingHom CoeffHom, const vector<RingElem>& IndetImages)
  {
    const ring& R = CoeffRing(owner(arg));
    if (IsZZ(R) || IsQQ(R) || IsFiniteField(R)) // WORKAROUND/BUG: should really test whether coeffs are mapped to coeffs
    {
      // Check for case where all indet images are just PPs
      bool SimpleCase = true;
      for (const auto& f: IndetImages)
        if (!IsMonomial(f) || !IsOne(LC(f)))  { SimpleCase = false; break; }
      if (SimpleCase)
        return ApplySPRCodomain_(image, arg, CoeffHom, IndetImages);
    }
    const SparsePolyRing S = owner(image);
    geobucket gbk(S);

    const long NumInd = len(IndetImages);
    vector<long> expv(NumInd); // workspace for exponents in loop below
    for (SparsePolyIter i=BeginIter(arg); !IsEnded(i); ++i)
    {
      RingElem SummandImage = CoeffHom(coeff(i));
      CoCoA_ASSERT(owner(SummandImage) == S);
      if (IsZero(SummandImage)) continue; // efficiency hack????
      exponents(expv, PP(i));
      for (long ind=0; ind < NumInd; ++ind)
      {
        if (expv[ind] == 0) continue;
        if (/*IsPPMonoidOv(myCodomain) && */IsMatrixOrdering(ordering(PPM(S))))
          if (expv[ind] >= 32749)  // BUG BUG: nasty hack to avoid exp overflow!!!
            CoCoA_THROW_ERROR1(ERR::ExpTooBig);
        SummandImage *= power(IndetImages[ind], expv[ind]);
      }
      gbk.myAddClear(SummandImage, NumTerms(SummandImage));
    }
    AddClear(image, gbk);
  }





  void ApplyGeneral(RingElem& image, ConstRefRingElem arg, const RingHom CoeffHom, const vector<RingElem>& IndetImages)
  {
    ring S = owner(image);
    const long NumInd = len(IndetImages);
    for (SparsePolyIter i=BeginIter(arg); !IsEnded(i); ++i)
    {
      RingElem SummandImage = CoeffHom(coeff(i));
      CoCoA_ASSERT(owner(SummandImage) == S);
      if (IsZero(SummandImage)) continue; // efficiency hack????
      ConstRefPPMonoidElem t(PP(i));
      for (long ind=0; ind < NumInd; ++ind)
      {
        const long d = exponent(t, ind); // ??? should we compute exponents?
        if (d == 0) continue;
        SummandImage *= power(IndetImages[ind], d);
      }
      S->myAdd(raw(image), raw(image), raw(SummandImage));
    }
  }
}  // end of anonymous namespace


  void SparsePolyRingBase::HomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    RingElem ans(myCodomain);  // Putting result into ans is exception safe and avoids aliasing problems.
    if ( IsSparsePolyRing(myCodomain) )
      ApplySPRCodomain(ans, RingElemAlias(myDomain, rawarg), myCoeffHom, myIndetImages);
    else
      ApplyGeneral(ans, RingElemAlias(myDomain, rawarg), myCoeffHom, myIndetImages);
    myCodomain->mySwap(rawimage, raw(ans));
  }


  void SparsePolyRingBase::HomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    const SparsePolyRing P = myDomain;
    if (NumIndets(P) == 0) return;
    out << " sending "
        << indets(myDomain) << " to " << myIndetImages;    
  }


  RingHom SparsePolyRingBase::myCoeffEmbeddingHomCtor() const
  {
    return RingHom(new CoeffEmbeddingHomImpl(SparsePolyRing(this)));
  }


  RingHom SparsePolyRingBase::myHomCtor(const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages) const
  {
    // Args already sanity checked by PolyRingHom (see PolyRing.C)
// DON'T KNOW IF I REALLY WANT TO MAKE THIS CHECK...
//       // Check to see if we're building an identity homomorphism
//       if (ring(this) == codomain && IsIdentity(CoeffHom))
//       {
//         bool IndetsFixed = true;
//         for (long i=0; i < myNumIndetsValue; ++i)
//           IndetsFixed &= (myIndetVector[i] == IndetImages[i]);
//         if (IndetsFixed) return IdentityHom(ring(this));
//       }
      // General case
    return RingHom(new HomImpl(SparsePolyRing(this), codomain, CoeffHom, IndetImages));
  }


  RingHom SparsePolyRingBase::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    vector<RingElem> IndetImages;
    for (long var=0; var < myNumIndets(); ++var)
      IndetImages.push_back(phi(theta(myIndets()[var])));

    return myHomCtor(codomain(phi), phi(theta(myCoeffEmbeddingHomCtor())), IndetImages);
  }


  bool SparsePolyRingBase::myImageLiesInSubfield(const RingHom& phi) const
  {
    CoCoA_THROW_ERROR1(ERR::NYI);
    return false;
  }


  //-- CoeffEmbeddingHomImpl ----------------------------------------

  //---------------------------------------------------------------------------
  // Functions for the class SparsePolyRingBase::CoeffEmbeddingHomImpl


  SparsePolyRingBase::CoeffEmbeddingHomImpl::CoeffEmbeddingHomImpl(const SparsePolyRing& P):
    RingHomEmbeddingBase(CoeffRing(P), P)
  {}


  void SparsePolyRingBase::CoeffEmbeddingHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    const SparsePolyRing P = myCodomain;
    RingElem ans(P);  // don't use image here for aliasing
    // ??? ANNA profile this:  (probably better to have myMonomial)
    if (!myDomain->myIsZero(rawarg))
      ans = monomial(P, RingElemAlias(myDomain, rawarg), one(PPM(P)));
    P->mySwap(rawimage, raw(ans));
  }


} // end of namespace CoCoA
