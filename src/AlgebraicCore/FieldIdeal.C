//   Copyright (c)  2005-2007,2009,2021  John Abbott and Anna M. Bigatti

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


#include "CoCoA/FieldIdeal.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"


//#include <vector>
using std::vector;

namespace CoCoA
{

  class FieldIdealImpl: public IdealBase
  {
  private:
    friend ideal NewFieldIdeal(const ring& k, const std::vector<RingElem>& gens);
    FieldIdealImpl(const ring& k, const vector<RingElem>& gens);
    // Default copy ctor works fine -- used by myClone(), needed by IdealBase::MakeUnique
    // Assignment disabled
    virtual ~FieldIdealImpl() {}
  public: // disable assignment
    FieldIdealImpl& operator=(const FieldIdealImpl&) = delete;

  public:
    virtual FieldIdealImpl* myClone() const override;
//???    virtual void swap(ideal& other);

    // functions every ideal must implement
    const ring& myRing() const override  { return myR; }
    const std::vector<RingElem>& myGens() const override  { return myGensValue;}
    const std::vector<RingElem>& myTidyGens(const CpuTimeLimit&) const override  { return myTidyGensValue; }

    bool IamZero() const override;
    bool IamOne() const override;
    bool IhaveElem(RingElemConstRawPtr rawx) const override;
    void myReduceMod(RingElemRawPtr rawx) const override;
    void myAdd(const ideal&) override;
    void myMul(const ideal&) override;
    void myIntersect(const ideal&) override;
    void myColon(const ideal&) override;
    bool myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const override; // lhs = num/den modulo the ideal  (lhs = 0 if quotient does not exist or ideal = ideal(1))

//???    void myOutputSelf(OpenMath::OutputChannel&) const override;
    static const FieldIdealImpl* ourGetPtr(const ideal&);
  protected: // more functions every ideal must implement
    void myTestIsMaximal() const override;
    void myTestIsPrimary() const override;
    void myTestIsPrime() const override;
    void myTestIsRadical() const override;

  private: // data members
    const ring myR;
    vector<RingElem> myGensValue;
    vector<RingElem> myTidyGensValue;
  };

  //---------------------------------------------------------------------------
  // Functions to do with FieldIdealImpl

  FieldIdealImpl::FieldIdealImpl(const ring& k, const vector<RingElem>& gens):
      myR(k),
      myGensValue(gens)
  {
    // NewFieldIdeal has already checked that the args are good.
    bool AllZero = true;
    const long n = len(gens);
    for (long i=0; i < n; ++i)
      if (!IsZero(gens[i])) { AllZero = false; break; }
    if (!AllZero) myTidyGensValue.push_back(one(myR));
    IamPrime3Flag = IamMaximal3Flag = AllZero;
  }


  FieldIdealImpl* FieldIdealImpl::myClone() const
  { return new FieldIdealImpl(*this); }


  inline const FieldIdealImpl* FieldIdealImpl::ourGetPtr(const ideal& J)
  { return dynamic_cast<const FieldIdealImpl*>(J.myIdealPtr()); }


  bool FieldIdealImpl::IamZero() const
  { return myTidyGensValue.empty(); }


  bool FieldIdealImpl::IamOne() const
  { return !myTidyGensValue.empty(); }


  void FieldIdealImpl::myTestIsMaximal() const
  {
    IamMaximal3Flag = IamZero();
    IamPrime3Flag = IamZero();
    IamPrimary3Flag = IamZero();
  }


  void FieldIdealImpl::myTestIsPrimary() const
  { myTestIsMaximal(); }


  void FieldIdealImpl::myTestIsPrime() const
  { myTestIsMaximal(); }


  void FieldIdealImpl::myTestIsRadical() const
  { myTestIsMaximal(); }


  bool FieldIdealImpl::IhaveElem(const RingElemConstRawPtr rawx) const
  {
    if (myR->myIsZero(rawx)) return true;
    return IamOne();
  }


  void FieldIdealImpl::myReduceMod(RingElemRawPtr rawx) const
  {
    if (IamZero()) return;
    myR->myAssignZero(rawx);
  }


  void FieldIdealImpl::myAdd(const ideal& J)
  {
    if (IamZero() && !ourGetPtr(J)->IamZero())
      myTidyGensValue.push_back(one(myR));
    // clever insert (skip 0s, deal with 1s) & check what this really does (Anna)
    myGensValue.insert(myGensValue.end(), gens(J).begin(), gens(J).end());
    IamPrime3Flag = IamMaximal3Flag = IamZero();
  }


  void FieldIdealImpl::myMul(const ideal& /*J*/)
  {
    if (IamZero()) return;
    CoCoA_THROW_ERROR(ERR::NYI, "FieldIdealImpl::myMul");
  }


  void FieldIdealImpl::myIntersect(const ideal& J)
  {
    if (IamOne() && ourGetPtr(J)->IamZero())
      myTidyGensValue.clear();
    myGensValue = myTidyGensValue;
    IamPrime3Flag = IamMaximal3Flag = IamZero();
  }


  void FieldIdealImpl::myColon(const ideal& J)
  {
    if (IamZero() && ourGetPtr(J)->IamZero())
      myTidyGensValue.push_back(one(myR));
    myGensValue = myTidyGensValue;
    IamPrime3Flag = IamMaximal3Flag = IamZero();
  }


  bool FieldIdealImpl::myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const
  {
    if (!IamZero()) myR->myAssignZero(rawlhs);
    else myR->myDiv(rawlhs, rawnum, rawden);
    return true;
  }


//   void FieldIdealImpl::myOutputSelf_OM(OpenMath::OutputChannel& OMOut) const
//   {
//     OMOut->OutputApplyStart();
//     OMOut->OutputSymbol(OpenMath::symbol("cocoa", "ideal"));
//     const long NumGens = len(myGensValue);
//     OMOut->OutputInteger(NumGens);   // Number of gens, should be an attribute???
//     for (long i=0; i < NumGens; ++i) // To be reconsidered ???
//       OMOut << myGensValue[i];       // ???
//     OMOut->OutputApplyEnd();
//   }


  ideal NewFieldIdeal(const ring& k, const std::vector<RingElem>& gens)
  {
    // Check that k is indeed a field, and that all gens belong to k
    if (!IsField(k))
      CoCoA_THROW_ERROR(ERR::ReqField, "NewFieldIdeal");
    const long n = len(gens);
    for (long i=0; i < n; ++i)
      if (owner(gens[i]) != k)
        CoCoA_THROW_ERROR(ERR::MixedRings, "NewFieldIdeal");
    return ideal(new FieldIdealImpl(k, gens));
  }

} // end of namespace CoCoA


// RCS header/log below
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/FieldIdeal.C,v 1.24 2024/07/11 15:42:42 bigatti Exp $
// $Log: FieldIdeal.C,v $
// Revision 1.24  2024/07/11 15:42:42  bigatti
// *** empty log message ***
//
// Revision 1.23  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.22  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.21  2021/10/30 19:36:35  abbott
// Summary: Added in some missing "override" (according to clang)
//
// Revision 1.20  2021/10/30 11:53:48  abbott
// Summary: Used keyword override (redmine 1625); a little tidying too
//
// Revision 1.19  2021/01/07 15:07:02  abbott
// Summary: Corrected copyright
//
// Revision 1.18  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.17  2020/02/12 09:01:47  bigatti
// -- changed myTestIsMaximal etc to return void (and consequences)
//
// Revision 1.16  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.15  2018/03/29 09:36:40  bigatti
// -- added member functions myTestIsRadical, myTestIsPrimary and flags
//
// Revision 1.14  2018/03/20 15:19:25  bigatti
// -- minor fix
//
// Revision 1.13  2018/03/20 11:38:08  bigatti
// -- changed iAm***Test --> myTestIs***;  and it returns bool
//
// Revision 1.12  2014/07/30 14:04:36  abbott
// Summary: Changed myAmbientRing into myRing
// Author: JAA
//
// Revision 1.11  2014/06/17 10:07:39  abbott
// Summary: Commented out unusued param; corrected error mesg
// Author: JAA
//
// Revision 1.10  2014/04/30 16:07:00  abbott
// Summary: Replaced size_t by long; Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.9  2013/06/28 17:03:51  abbott
// Modified semantics of IdealBase::myDivMod;
// it now returns a boolean.
// Several consequential changes.
//
// Revision 1.8  2012/05/30 16:04:55  bigatti
// -- applied "3" convention on bool3 functions and member fields
//
// Revision 1.7  2012/01/26 16:47:00  bigatti
// -- changed back_inserter into insert
//
// Revision 1.6  2010/06/10 08:00:02  bigatti
// -- fixed naming conventions
//
// Revision 1.5  2010/02/04 09:57:11  bigatti
// -- added "mul" for ideals.  Implemented only for SparsePolyRing
//
// Revision 1.4  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.3  2009/09/24 14:26:45  abbott
// Removed pointless include directive.
// Removed some unnecessary "std::" prefixes.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.3  2007/01/15 13:35:14  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.2  2006/11/27 16:18:33  cocoa
// -- moved classes declarations from .H to .C (DenseMatrix, DiagMatrix,
//    FieldIdeal, SpecialMatrix)
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.3  2006/03/21 09:43:14  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
//
// Revision 1.2  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.4  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.3  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.2  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.1  2004/11/05 15:30:57  cocoa
// Separated FieldIdealImpl from the "abstract" classes ideal and
// IdealBase.
//
