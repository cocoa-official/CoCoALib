//   Copyright (c)  2018  John Abbott and Anna M. Bigatti
//   Authors:  2018  John Abbott and Anna M. Bigatti

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


// Source code for abstract class SparsePolyRing and friends

#include "CoCoA/SparsePolyOps-involutive.H"

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
#include "CoCoA/RingZZ.H" // for IsZZ
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SparsePolyOps-MinPoly.H" // for MinPolyDef
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpUniversalInvolutiveBasisContainer.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"     // for myGcd
#include "CoCoA/matrix.H" // for OrdMat, myDivMod
#include "CoCoA/random.H" // for RandomLongStream
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

  //-- namespace Involutive functions --------------------------------

  const std::vector<RingElem>& Involutive::JanetBasis(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myJanetBasis();
  }


  bool Involutive::IsDeltaRegular(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->InvIamDeltaRegular();
  }


  bool Involutive::IsMonomial(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->InvIamMonomial();
  }


  bool Involutive::IsHomogeneous(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->InvIamHomogeneous();
  }


  bool Involutive::IsCohenMacaulay(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->InvIamCohenMacaulay();
  }


  std::map<PPMonoidElem, std::vector<bool> > Involutive::MultVars(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    std::map<PPMonoidElem, std::vector<bool> > res;
    ptrI->InvMultVars(res);
    return res;
  }


  std::map<PPMonoidElem, std::vector<bool> > Involutive::NonMultVars(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    std::map<PPMonoidElem, std::vector<bool> > res;
    ptrI->InvNonMultVars(res);
    return res;
  }


  RingElem Involutive::HilbertPol(const ideal& I, ConstRefRingElem var)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    RingElem res;
    ptrI->InvHilbertPol(res, var);
    return res;
  }


  RingElem Involutive::HilbertSeries(const ideal& I, ConstRefRingElem var)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    RingElem res;
    ptrI->InvHilbertSeries(res, var);
    return res;
  }


  FGModule Involutive::FirstSyzygy(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    FGModule syz = NewFreeModule(RingOf(I), 1); // fake initialization
    ptrI->InvFirstSyzygy(syz);
    return syz;
  }


  long Involutive::Dimension(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->InvDimension();
  }


  std::vector<std::pair<PPMonoidElem, std::vector<bool> > > Involutive::ComplementaryDecomposition(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    std::vector<std::pair<PPMonoidElem, std::vector<bool> > > res;
    ptrI->InvComplementaryDecomposition(res);
    return res;
  }


  long Involutive::Depth(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->InvDepth();
  }


  long Involutive::ProjDim(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->InvProjDim();
  }


  std::vector<RingElem> Involutive::Socle(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    std::vector<RingElem> res;
    ptrI->InvSocle(res);
    return res;
  }


  std::map<std::pair<long, long>, long> Involutive::ExtremalBettiNumbers(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    std::map<std::pair<long, long>, long> res;
    ptrI->InvExtremalBettiNumbers(res);
    return res;
  }


  std::vector<RingElem> Involutive::RegularSequence(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    std::vector<RingElem> res;
    ptrI->InvRegularSequence(res);
    return res;
  }


  std::vector<RingElem> Involutive::MaximalStronglyIndependentSet(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    std::vector<RingElem> res;
    ptrI->InvMaximalStronglyIndependentSet(res);
    return res;
  }


  long Involutive::Regularity(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->InvRegularity();
  }


  long Involutive::Satiety(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->InvSatiety();
  }


  std::vector<RingElem> Involutive::Saturation(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRingBase::IdealImpl* const ptrI =
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    std::vector<RingElem> res;
    ptrI->InvSaturation(res);
    return res;
  }

  ////-- member functions ---------------------------------------------

  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myJanetBasis() const
  { return myInvBasisContainerPtr->myJanetBasis(); }


  bool SparsePolyRingBase::IdealImpl::InvIamDeltaRegular() const
  { return myInvBasisContainerPtr->IamDeltaRegular(); }


  bool SparsePolyRingBase::IdealImpl::InvIamMonomial() const
  {return myInvBasisContainerPtr->IamMonomial(); }


  bool SparsePolyRingBase::IdealImpl::InvIamHomogeneous() const
  { return myInvBasisContainerPtr->IamHomogeneous(); }


  bool SparsePolyRingBase::IdealImpl::InvIamCohenMacaulay() const
  { return myInvBasisContainerPtr->IamCohenMacaulay(); }


  void SparsePolyRingBase::IdealImpl::InvMultVars(std::map<PPMonoidElem, std::vector<bool> >& MultVars) const
  { MultVars = myInvBasisContainerPtr->myMultVars(); }


  void SparsePolyRingBase::IdealImpl::InvNonMultVars(std::map<PPMonoidElem, std::vector<bool> >& NonMultVars) const
  { NonMultVars = myInvBasisContainerPtr->myNonMultVars(); }


  void SparsePolyRingBase::IdealImpl::InvHilbertPol(RingElem& HilbertPol, ConstRefRingElem var) const
  { HilbertPol = myInvBasisContainerPtr->myHilbertPol(var); }


  void SparsePolyRingBase::IdealImpl::InvHilbertSeries(RingElem& HilbertSeries, ConstRefRingElem var) const
  { HilbertSeries = myInvBasisContainerPtr->myHilbertSeries(var); }


  void SparsePolyRingBase::IdealImpl::InvFirstSyzygy(FGModule& syz) const
  { syz = myInvBasisContainerPtr->myFirstSyzygy(); }


  long SparsePolyRingBase::IdealImpl::InvDimension() const
  { return myInvBasisContainerPtr->myDimension(); }


  void SparsePolyRingBase::IdealImpl::InvComplementaryDecomposition(std::vector<std::pair<PPMonoidElem, std::vector<bool> > >& CompDecomp) const
  { CompDecomp = myInvBasisContainerPtr->myComplementaryDecomposition(); }


  long SparsePolyRingBase::IdealImpl::InvDepth() const
  { return myInvBasisContainerPtr->myDepth(); }


  long SparsePolyRingBase::IdealImpl::InvProjDim() const
  { return myInvBasisContainerPtr->myProjDim(); }


  void SparsePolyRingBase::IdealImpl::InvSocle(std::vector<RingElem>& socle) const
  { socle = myInvBasisContainerPtr->mySocle(); }


  void SparsePolyRingBase::IdealImpl::InvExtremalBettiNumbers(std::map<std::pair<long, long>, long>& ExtremalBettiNumbers) const
  { ExtremalBettiNumbers = myInvBasisContainerPtr->myExtremalBettiNumbers(); }


  void SparsePolyRingBase::IdealImpl::InvRegularSequence(std::vector<RingElem>& RegSeq) const
  { RegSeq = myInvBasisContainerPtr->myRegularSequence(); }


  void SparsePolyRingBase::IdealImpl::InvMaximalStronglyIndependentSet(std::vector<RingElem>& MaxSet) const
  { MaxSet = myInvBasisContainerPtr->myMaximalStronglyIndependentSet(); }


  long SparsePolyRingBase::IdealImpl::InvRegularity() const
  { return myInvBasisContainerPtr->myRegularity(); }


  long SparsePolyRingBase::IdealImpl::InvSatiety() const
  { return myInvBasisContainerPtr->mySatiety(); }


  void SparsePolyRingBase::IdealImpl::InvSaturation(std::vector<RingElem>& saturation) const
  { saturation = myInvBasisContainerPtr->mySaturation(); }



} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-involutive.C,v 1.10 2024/08/01 16:49:18 bigatti Exp $
// $Log: SparsePolyOps-involutive.C,v $
// Revision 1.10  2024/08/01 16:49:18  bigatti
// Summary: error codes: Req(Elem)SparsePolyRing;
// now using CoCoA_THROW_ERROR1/2
//
// Revision 1.9  2024/03/18 20:57:28  abbott
// Summary: Removed const from return value type (redmine 1606)
//
// Revision 1.8  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.7  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.6  2020/01/26 14:37:24  abbott
// Summary: Removed useless include (of NumTheory)
//
// Revision 1.5  2018/05/18 12:23:50  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.4  2018/05/17 15:50:34  bigatti
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.3  2018/04/10 14:20:47  bigatti
// -- fixed includes
//
// Revision 1.2  2018/03/15 14:18:06  bigatti
// -- added files SparsePolyOps-ideal.H and SparsePolyOps-involutive.H
//
// Revision 1.1  2018/03/12 14:39:07  bigatti
// -- Renamed SparsePoly-ideal/involutive into SparsePolyOps-ideal/involutive
//
// Revision 1.1  2018/02/27 17:25:39  bigatti
// -- split from SparsePolyRing
//
// Revision 1.1  2018/02/27 17:12:37  bigatti
// -- Renamed SparsePolyRing_ideal.C --> SparsePoly-ideal.C
//
