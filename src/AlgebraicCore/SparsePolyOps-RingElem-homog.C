//   Copyright (c)  2005-2018,2020  John Abbott and Anna M. Bigatti
//   Authors:  2005-2018,2020  John Abbott and Anna M. Bigatti

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


// Source code for fns related to homogeneity in SparsePolyRing

#include "CoCoA/SparsePolyOps-RingElem.H"

#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/degree.H"
#include "CoCoA/geobucket.H"


namespace CoCoA
{

  namespace // anonymous
  {

    // These 2 are duplicates of fns in SparsePolyOps-RingElem.C
    
    inline void CheckCompatible(ConstRefRingElem x, ConstRefRingElem y, const char* const FnName)
    {
      if (owner(x) != owner(y))  CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    }


    inline void CheckElemSparsePolyRing(ConstRefRingElem f, const char* const FnName)
    {
      if (!IsSparsePolyRing(owner(f))) CoCoA_THROW_ERROR(ERR::ReqElemSparsePolyRing, FnName);
    }

  } // end of namespace anonymous


  degree wdeg(ConstRefRingElem f)
  {
    CheckElemSparsePolyRing(f, "wdeg(f)");
    if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    return wdeg(LPP(f)); // yes! see [KR] introduction Sec.4.3
  }


  int CmpWDeg(ConstRefRingElem f1, ConstRefRingElem f2)
  {
    CheckCompatible(f1, f2, "CmpWDeg(f1,f2)");
    CheckElemSparsePolyRing(f1, "CmpWDeg(f1,f2)");
    if (IsZero(f1) || IsZero(f2))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    // Now we know that both f1 and f2 are non-zero
    //    return SparsePolyRingPtr(owner(f1))->myCmpWDeg(raw(f1), raw(f2));
    return CmpWDeg(LPP(f1), LPP(f2));
  }


  int CmpWDegPartial(ConstRefRingElem f1, ConstRefRingElem f2, long i)  // assumes 0 <= i
  {
    CheckCompatible(f1, f2, "CmpWDegPartial(f1,f2,i)");
    CheckElemSparsePolyRing(f1, "CmpWDegPartial(f1,f2,i)");
    if (IsZero(f1) || IsZero(f2))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    // Now we know that both f1 and f2 are non-zero
    //    return SparsePolyRingPtr(owner(f1))->myCmpWDeg(raw(f1), raw(f2));
    return CmpWDegPartial(LPP(f1), LPP(f2), i);
  }


  bool IsHomog(ConstRefRingElem f)
  {
    CheckElemSparsePolyRing(f, "IsHomog(f)");
    if (GradingDim(owner(f))==0)  CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    return SparsePolyRingPtr(owner(f))->myIsHomog(raw(f));
  }


  bool IsHomogPartial(ConstRefRingElem f, long n)  // assumes n >= 0
  {
    CheckElemSparsePolyRing(f, "IsHomogPartial(f,n)");
    return SparsePolyRingPtr(owner(f))->myIsHomogPartial(raw(f), n);
  }


  RingElem homog(ConstRefRingElem f, ConstRefRingElem h)
  {
    const char* const FnName = "homog(RingElem, RingElem)";
    CheckCompatible(f, h, FnName);
    const SparsePolyRing P = owner(f);
    if ( GradingDim(P)!=1 )  CoCoA_THROW_ERROR2(ERR::BadArg, "GrDim must be 1");
    if ( !IsIndet(h) )  CoCoA_THROW_ERROR2(ERR::ReqIndet, "arg 2");
    if ( wdeg(h)[0]!=1 )  CoCoA_THROW_ERROR2(ERR::BadArg, "degree of hom.indet must be 1");
    RingElem fHom(P);
    P->myHomog(raw(fHom), raw(f), raw(h));
    return fHom;
  }


  RingElem LF(ConstRefRingElem f)
  {
    if (!IsSparsePolyRing(owner(f)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRing& P = owner(f);
    if (GradingDim(P) == 0)  CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);

    RingElem LeadingForm(P);
    ConstRefPPMonoidElem LPPf(LPP(f));
    for (SparsePolyIter it=BeginIter(f) ; !IsEnded(it) ; ++it )
    {
      if (CmpWDeg(PP(it), LPPf) != 0) break;
      PushBack(LeadingForm, coeff(it), PP(it));
    }
    return LeadingForm;
  }


  RingElem CutLF(RingElem& f) // MODIFIES f
  {
    if (!IsSparsePolyRing(owner(f)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const SparsePolyRing& P = owner(f);
    if (GradingDim(P) == 0)  CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);

    RingElem ans(P);
    do
    {
      P->myMoveLMToBack(raw(ans), raw(f));
    }
    while (!IsZero(f) && (CmpWDeg(LPP(f), LPP(ans)) == 0));
    return ans;
  }


  bool SparsePolyRingBase::myIsHomogPartial(ConstRawPtr rawf, long n) const  // assumes 0 <= n <= GrDim
  {
    CoCoA_ASSERT(0 <= n && n <= myGradingDim());
    if (myIsZero(rawf)) { return true; }
    SparsePolyIter itf=myBeginIter(rawf);
    const PPMonoidElem FirstPP=PP(itf);
    for (++itf; !IsEnded(itf); ++itf)
    {
      CoCoA_ASSERT( cmp(FirstPP, PP(itf)) > 0 ); // assert f is correctly sorted
      if ( CmpWDegPartial(FirstPP, PP(itf), n) != 0 )  return false;
    }
    return true;
  }

  bool SparsePolyRingBase::myIsHomog(ConstRawPtr rawf) const
  {
    if (myGradingDim()==0)  CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    if (myIsZero(rawf))  return true;
    SparsePolyIter itf=myBeginIter(rawf);
    const PPMonoidElem FirstPP=PP(itf);
    for (++itf; !IsEnded(itf); ++itf)
    {
      CoCoA_ASSERT( cmp(FirstPP, PP(itf))>0 ); // assert f is correctly sorted
      if ( CmpWDeg(FirstPP, PP(itf))!=0 )  return false;
    }
    return true;
  }


  void SparsePolyRingBase::myHomog(RawPtr rawfHom, ConstRawPtr rawf, ConstRawPtr rawh) const
  {
    const SparsePolyRing P(this);
    RingElemAlias h(P,rawh);
    CoCoA_ASSERT( myGradingDim()==1 );
    CoCoA_ASSERT( IsIndet(h) );
    CoCoA_ASSERT( IsOne(wdeg(h)[0]) );
    geobucket gbk(P);  // accumulate result into a geobucket for speed
    RingElemAlias f(P,rawf);
    ConstRefPPMonoidElem PPh = LPP(h);
    const BigInt d = wdeg(f)[0];
    for (SparsePolyIter it=BeginIter(f) ; !IsEnded(it) ; ++it )
    {
      RingElem term = monomial(P, coeff(it),PP(it) * power(PPh,d-wdeg(PP(it))[0]));
      gbk.myAddClear(term,1);
    }
    RingElem tmp(P);  // for exception safety
    AddClear(tmp, gbk);
    mySwap(rawfHom, raw(tmp));
  }


  RingElem HomogCompt(ConstRefRingElem f, long d)
  {
    const ring& P = owner(f);
    if (!IsSparsePolyRing(P))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (GradingDim(P) != 1)  CoCoA_THROW_ERROR2(ERR::BadArg, "GradingDim must be 1");
    if (d < 0)  CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "arg 2");
    RingElem form(P);
    for (auto it=BeginIter(f); !IsEnded(it); ++it)
    {
      const long DegPP = ConvertTo<long>(wdeg(PP(it))[0]);
//      if (DegPP != d) continue; // ASSUME term-ord is wdeg compatible!!
      if (DegPP > d) continue;
      if (DegPP < d) break;
      PushBack(form, coeff(it), PP(it));
    }
    return form;
  }


} // end of namespace CoCoA
