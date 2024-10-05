//   Copyright (c)  2010,2021  John Abbott and Anna M. Bigatti

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


#include "CoCoA/PPMonoidHom.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
#include<vector>
using std::vector; // used only in IsInKer

namespace CoCoA
{

  PPMonoidElem PPMonoidHom::operator()(ConstRefPPMonoidElem x) const
  {
    if (owner(x) != domain(*this))
      CoCoA_THROW_ERROR(ERR::BadPPMonoidHomArg, "Applying PPMonoidHom to PPMonoidElem");
    PPMonoidElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(x));
    return ans;
  }


  PPMonoidHom PPMonoidHom::operator()(const PPMonoidHom& theta) const
  {
    if (codomain(theta) != domain(*this))
      CoCoA_THROW_ERROR(ERR::BadCompose, "PPMonoidHom(PPMonoidHom)  i.e. PPMonoidHom composition");
    CoCoA_THROW_ERROR(ERR::NYI,"COMPOSE for PPMonoidHom");return theta;
//    return domain(theta)->myCompose(*this, theta);
  }


  void PPMonoidHomBase::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "PPMonoidHom(" << myDomain << " --> " << myCodomain;
    myOutputSelfDetails(out);
    out << ")";
  }


  void PPMonoidHomBase::myOutputSelfDetails(std::ostream& /*out*/) const
  {
    // Default definition does nothing (as there are no extra details to print).
    // SHOULD THIS BE PURE VIRTUAL (with no default defn)???
  }


  std::ostream& operator<<(std::ostream& out, const PPMonoidHom& phi)
  {
    if (!out) return out;  // short-cut for bad ostreams
    phi->myOutputSelf(out);
    return out;
  }



  //---------------------------------------------------------------------------

  class IdentityPPMonoidHomImpl: public PPMonoidHomBase
  {
  private:
    explicit IdentityPPMonoidHomImpl(const PPMonoid& R);
    friend PPMonoidHom IdentityHom(const PPMonoid& R); // The only function that calls the ctor.
  public: // fns every PPMonoidHom must implement
    void myApply(PPMonoidElemRawPtr image, PPMonoidElemConstRawPtr arg) const override;
    void myOutputSelfDetails(std::ostream& out) const override;
  };


  IdentityPPMonoidHomImpl::IdentityPPMonoidHomImpl(const PPMonoid& R):
      PPMonoidHomBase(R, R)
  {}


  void IdentityPPMonoidHomImpl::myApply(PPMonoidElemRawPtr image, PPMonoidElemConstRawPtr arg) const
  {
    myDomain->myAssign(image, arg);
  }


  void IdentityPPMonoidHomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << ": the identity";
  }


  PPMonoidHom IdentityHom(const PPMonoid& R)
  {
    return PPMonoidHom(new IdentityPPMonoidHomImpl(R));
  }


  //-----------------------------------------------------------------------------

  class GeneralPPMonoidHomImpl:  public PPMonoidHomBase
  {
  private:
    explicit GeneralPPMonoidHomImpl(const PPMonoid& domain, const std::vector<PPMonoidElem>& images);
    friend PPMonoidHom GeneralHom(const PPMonoid& PPM, const std::vector<PPMonoidElem>& images); // The only function that calls the ctor.
  public: // fns every PPMonoidHom must implement
    void myApply(PPMonoidElemRawPtr image, PPMonoidElemConstRawPtr arg) const override;
    void myOutputSelfDetails(std::ostream& out) const override;
  private: // data members
    const vector<PPMonoidElem> myImages;
  };

  GeneralPPMonoidHomImpl::GeneralPPMonoidHomImpl(const PPMonoid& PPM, const vector<PPMonoidElem>& images):
      PPMonoidHomBase(PPM, owner(images[0])),
      myImages(images)
  {
    CoCoA_ASSERT(len(images) == NumIndets(PPM));
  }


  void GeneralPPMonoidHomImpl::myApply(PPMonoidElemRawPtr image, PPMonoidElemConstRawPtr arg) const
  {
    PPMonoidElem ans(myCodomain);
    const long N = NumIndets(myDomain);
    vector<long> exp(N);
    myDomain->myExponents(exp, arg);
    for (long i=0; i < N; ++i)
    {
      if (/*IsPPMonoidOv(myCodomain)&&*/
          IsMatrixOrdering(ordering(myCodomain)))
        if (exp[i] >= 32749)  // BUG BUG BUG  now matrix ordering allows higher exponents!!!
          CoCoA_THROW_ERROR(ERR::ExpTooBig, "GeneralPPMonoidHomImpl::myApply");
      ans *= power(myImages[i], exp[i]);
    }
    myCodomain->myAssign(image, raw(ans));
  }


  void GeneralPPMonoidHomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << ": ";
    const long N = NumIndets(myDomain);
    for  (long i=0; i < N; ++i)
    {
      if (i > 0) out << ", ";
      out << indet(myDomain, i) << " |--> " << myImages[i];
    }
  }


  PPMonoidHom GeneralHom(const PPMonoid& PPM, const std::vector<PPMonoidElem>& images)
  {
    if (len(images) != NumIndets(PPM))
      CoCoA_THROW_ERROR(ERR::BadArg, "PPMonoid GeneralHom -- wrong number of images");
    return PPMonoidHom(new GeneralPPMonoidHomImpl(PPM, images));
  }



  // ----------------------------------------------------------------------
  class RestrictionPPMonoidHomImpl:  public PPMonoidHomBase
  {
  private:
    explicit RestrictionPPMonoidHomImpl(const PPMonoid& domain, const std::vector<bool>& mask);
    friend PPMonoidHom RestrictionHom(const PPMonoid& PPM, const std::vector<long>& IndetIndices); // The only function that calls the ctor.
  public: // fns every PPMonoidHom must implement
    void myApply(PPMonoidElemRawPtr image, PPMonoidElemConstRawPtr arg) const override;
    void myOutputSelfDetails(std::ostream& out) const override;
  private: // data members
    const vector<bool> myMask;
  };


  RestrictionPPMonoidHomImpl::RestrictionPPMonoidHomImpl(const PPMonoid& PPM, const vector<bool>& mask):
      PPMonoidHomBase(PPM, PPM),
      myMask(mask)
  {
    CoCoA_ASSERT(len(mask) == NumIndets(PPM));
  }


  void RestrictionPPMonoidHomImpl::myApply(PPMonoidElemRawPtr image, PPMonoidElemConstRawPtr arg) const
  {
    PPMonoidElem ans(myCodomain);
    const long N = NumIndets(myDomain);
    vector<long> exp(N);
    myDomain->myExponents(exp, arg);  // BUG BUG BUG   should use BigExponents   BUG BUG BUG!!!!
    for (long i=0; i < N; ++i)
      if (myMask[i] == false)
        exp[i] = 0;
    myCodomain->myAssign(image, exp);
  }


  void RestrictionPPMonoidHomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << ": restriction hom retaining";
    const long N = NumIndets(myDomain);
    for  (long i=0; i < N; ++i)
    {
      if (myMask[i])
        out << " " << indet(myDomain, i);
    }
  }



  PPMonoidHom RestrictionHom(const PPMonoid& PPM, const std::vector<long>& IndetIndices)
  {
    // Check that entries in IndetIndices are in range -- currently we do not accept duplicates.
    const long nvars = NumIndets(PPM);
    vector<bool> mask(nvars);
    const long n = len(IndetIndices);
    for (long j=0; j < n; ++j)
    {
      const long i = IndetIndices[j];
      if (i < 0 || i >= nvars || mask[i] == true)
        CoCoA_THROW_ERROR(ERR::BadIndex, "RestrictionHom pseudo ctor");
      mask[i] = true;
    }
    return PPMonoidHom(new RestrictionPPMonoidHomImpl(PPM, mask));
  }


} // end of namespace CoCoA


// RCS header/log
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PPMonoidHom.C,v 1.15 2024/03/08 19:59:59 abbott Exp $
// $Log: PPMonoidHom.C,v $
// Revision 1.15  2024/03/08 19:59:59  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.14  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.13  2021/10/30 11:53:48  abbott
// Summary: Used keyword override (redmine 1625); a little tidying too
//
// Revision 1.12  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.11  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.10  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.9  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.8  2019/10/14 16:00:08  bigatti
// -- added check if exponent too big (IsMatrixOrdering)
//
// Revision 1.7  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.6  2012/02/24 13:09:56  abbott
// Added (missing) const.
//
// Revision 1.5  2012/02/14 15:15:43  bigatti
// -- fixed index bug in RestrictionHom
//
// Revision 1.4  2012/02/10 17:08:06  abbott
// Added new pseudo-ctor for RestrictionHom (& related impl class).
//
// Revision 1.3  2011/11/09 14:09:53  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.2  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.1  2010/07/09 17:03:42  abbott
// First simple implementation of PPMonoid homs.
//
//
