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
  // ??? appropriate use of inheritance here?  this is getting pretty hugly

  // assume image==0
  void ApplySPRCodomain(RingElem& image, ConstRefRingElem arg, const RingHom CoeffHom, const vector<RingElem>& IndetImages)
  {
    const SparsePolyRing S = owner(image);
    geobucket gbk(S);

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
        if (/*IsPPMonoidOv(myCodomain) && */IsMatrixOrdering(ordering(PPM(S))))
          if (d >= 32749)  // BUG BUG: nasty hack to avoid exp overflow!!!
            CoCoA_THROW_ERROR(ERR::ExpTooBig, "ApplySPRCodomain");
        SummandImage *= power(IndetImages[ind], d);
      }
      //        SparsePolyRingPtr(S)->myAddClear(raw(ans), raw(SummandImage));
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
    CoCoA_THROW_ERROR(ERR::NYI, "SparsePolyRingBase::myImageLiesInSubfield");
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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-RingHom.C,v 1.11 2024/01/19 15:45:21 bigatti Exp $
// $Log: SparsePolyOps-RingHom.C,v $
// Revision 1.11  2024/01/19 15:45:21  bigatti
// Summary: Redmine #1623: modified printing for RingHom
//
// Revision 1.10  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.9  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.8  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.7  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.6  2020/01/26 14:37:24  abbott
// Summary: Removed useless include (of NumTheory)
//
// Revision 1.5  2018/07/26 15:34:05  abbott
// Summary: Added comment
//
// Revision 1.4  2018/07/26 15:21:05  bigatti
// -- checking exponent overflow in ApplySPRCodomain
//
// Revision 1.3  2018/05/18 12:23:50  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.2  2018/05/17 15:47:54  bigatti
// -- added include SparsePolyIter
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.1  2018/04/16 21:47:36  bigatti
// - added SparsePolyOps-RingHom
//
