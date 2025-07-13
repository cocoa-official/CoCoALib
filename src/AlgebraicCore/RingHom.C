//   Copyright (c)  2005,2021  John Abbott and Anna M. Bigatti

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


#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"
#include "CoCoA/matrix.H"

#include <iostream>
using std::ostream;
#include<vector>
using std::vector; // used only in IsInKer

namespace CoCoA
{

  //-------------------------------------------------------
  // CLASS PartialRingHom

  RingElem PartialRingHom::operator()(ConstRefRingElem x) const
  {
    if (owner(x) != domain(*this))
      CoCoA_THROW_ERROR1(ERR::BadRingHomArg);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(x));
    return ans;
  }


  RingElem PartialRingHom::operator()(const MachineInt& n) const
  {
    RingElem arg(domain(*this), n);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }

  RingElem PartialRingHom::operator()(const BigInt& N) const
  {
    RingElem arg(domain(*this), N);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }

  RingElem PartialRingHom::operator()(const BigRat& q) const
  {
    RingElem arg(domain(*this), q);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }


  //PROTOTYPE IMPL!!!!!!!!
  //Result ought to be dense/sparse/diag according as M is?
  matrix PartialRingHom::operator()(ConstMatrixView M) const
  {
    CoCoA_ASSERT(domain(*this) == RingOf(M));
    RingElem ans(codomain(*this));
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    matrix NewM(NewDenseMat(codomain(*this), nrows, ncols));
    for (long i=0; i < nrows; ++i)
      for (long j=0; j < ncols; ++j)
      {
        mySmartPtr->myApply(raw(ans), raw(M(i,j)));
        SetEntry(NewM, i, j, ans);
      }
    return NewM;
  }


  PartialRingHom PartialRingHom::operator()(const RingHom& theta) const
  {
    if (codomain(theta) != domain(*this))
      CoCoA_THROW_ERROR1(ERR::BadCompose);
// BUG BUG BUG can be clever here!!!
    return sequential(*this, theta);
  }


  PartialRingHom PartialRingHom::operator()(const PartialRingHom& theta) const
  {
    if (codomain(theta) != domain(*this))
      CoCoA_THROW_ERROR1(ERR::BadCompose);
    return sequential(*this, theta);
  }


  void RingHomBase::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "PartialRingHom(" << myDomain << " -> " << myCodomain;
    myOutputSelfDetails(out);
    out << ")";
  }


  void RingHomBase::myOutputSelfDetails(std::ostream& /*out*/) const
  {
    // Default definition does nothing (as there are no extra details to print).
    // SHOULD THIS BE PURE VIRTUAL (with no default defn)???
  }


  std::ostream& operator<<(std::ostream& out, const PartialRingHom& phi)
  {
    if (!out) return out;  // short-cut for bad ostreams
/// anna: should it be PartialRingHom?
    out << "RingHom(" << domain(phi) << " -> " << codomain(phi);
    phi->myOutputSelfDetails(out);
    out << ")";
    return out;
  }


  // better to use indices or iterators in this fn???
// USE AN ALGORITHM!!!
  bool IsInKer(const ideal& I, const PartialRingHom& phi)
  {
//    const auto& G = gens(I);
//    return std::all_of(G.begin(),G.end(), [](const RingElem& g){return (IsZero(phi(g)));}); // needs <algorithm>
    for (const RingElem& g: gens(I))
    {
      if (!IsZero(phi(g))) return false;
    }
    return true;
  }


  //-------------------------------------------------------
  // CLASS RingHom

  RingElem RingHom::operator()(ConstRefRingElem x) const
  {
    if (owner(x) != domain(*this))
      CoCoA_THROW_ERROR1(ERR::BadRingHomArg);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(x));
    return ans;
  }


  RingElem RingHom::operator()(const MachineInt& n) const
  {
    RingElem arg(domain(*this), n);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }

  RingElem RingHom::operator()(const BigInt& N) const
  {
    RingElem arg(domain(*this), N);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }

  RingElem RingHom::operator()(const BigRat& q) const
  {
    RingElem arg(domain(*this), q);
    RingElem ans(codomain(*this));
    mySmartPtr->myApply(raw(ans), raw(arg));
    return ans;
  }

  // simple but ugly impl
  std::vector<RingElem> RingHom::operator()(const std::vector<RingElem>& v) const
  {
    const long n = len(v);
    const ring& R = domain(*this);
    std::vector<RingElem> ans(n,zero(codomain(*this)));
    for (long i=0; i < n; ++i)
    {
      if (owner(v[i]) != R)  CoCoA_THROW_ERROR1(ERR::BadRingHomArg);
      mySmartPtr->myApply(raw(ans[i]),raw(v[i]));
    }
    return ans; // std::move?
  }


  //PROTOTYPE IMPL!!!!!!!!  (see above PartialRingHom)
  //Result ought to be dense/sparse/diag according as M is?
  matrix RingHom::operator()(ConstMatrixView M) const
  {
    CoCoA_ASSERT(domain(*this) == RingOf(M));
    RingElem ans(codomain(*this));
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    matrix NewM(NewDenseMat(codomain(*this), nrows, ncols));
    for (long i=0; i < nrows; ++i)
      for (long j=0; j < ncols; ++j)
      {
        mySmartPtr->myApply(raw(ans), raw(M(i,j)));
        SetEntry(NewM, i, j, ans);
      }
    return NewM;
  }


  RingHom RingHom::operator()(const RingHom& theta) const
  {
    if (codomain(theta) != domain(*this))
      CoCoA_THROW_ERROR1(ERR::BadCompose);
    return domain(theta)->myCompose(*this, theta);
  }


  PartialRingHom RingHom::operator()(const PartialRingHom& theta) const
  {
    if (codomain(theta) != domain(*this))
      CoCoA_THROW_ERROR1(ERR::BadCompose);
    return sequential(*this, theta);
  }


//   void RingHomBase::myOutputSelf(std::ostream& out) const
//   {
//     out << "RingHom(" << myDomain << " --> " << myCodomain;
//     myOutputSelfDetails(out);
//     out << ")";
//   }


//   void RingHomBase::myOutputSelfDetails(std::ostream& /*out*/) const
//   {
//     // Default definition does nothing (as there are no extra details to print).
//     // SHOULD THIS BE PURE VIRTUAL (with no default defn)???
//   }

////// why there are 3 printing functions?
  std::ostream& operator<<(std::ostream& out, const RingHom& phi)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "RingHom(" << domain(phi) << " -> " << codomain(phi);
    phi->myOutputSelfDetails(out);
    out << ")";
    return out;
  }


  // better to use indices or iterators in this fn???
// USE AN ALGORITHM!!!
  bool IsInKer(const ideal& I, const RingHom& phi)
  {
//    const auto& G = gens(I);
//    return std::all_of(G.begin(),G.end(), [](const RingElem& g){return (IsZero(phi(g)));}); // needs <algorithm>
    for (const RingElem& g: gens(I))
    {
      if (!IsZero(phi(g))) return false;
    }
    return true;
  }



  //---------------------------------------------------------------------------

  class IdentityRingHomImpl: public RingHomBase
  {
  private:
    explicit IdentityRingHomImpl(const ring& R);
    friend RingHom IdentityHom(const ring& R); // The only function that calls the ctor.
  public: // fns every RingHom must implement
    void myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const override;
    bool IamPartial() const override  { return false; }
    void myOutputSelfDetails(std::ostream& out) const override;
  };


  IdentityRingHomImpl::IdentityRingHomImpl(const ring& R):
      RingHomBase(R, R)
  {}


  void IdentityRingHomImpl::myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const
  {
    myDomain->myAssign(image, arg);
  }


  void IdentityRingHomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << " the identity";
  }


  RingHom IdentityHom(const ring& R)
  {
    return RingHom(new IdentityRingHomImpl(R));
  }


  //---------------------------------------------------------------------------

  class SequentialRingHomImpl: public RingHomBase
  {
  private:
    SequentialRingHomImpl(const PartialRingHom& phi, const PartialRingHom& theta);
    friend PartialRingHom sequential(const PartialRingHom& phi, const PartialRingHom& theta); // The only function that calls the ctor.
  public: // fns every RingHom must implement
    void myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const override;
    bool IamPartial() const override  { return IsPartial(myPhi) || IsPartial(myTheta); }
    void myOutputSelfDetails(std::ostream& out) const override;
  private: // data members
    PartialRingHom myPhi;
    PartialRingHom myTheta;
  };


  SequentialRingHomImpl::SequentialRingHomImpl(const PartialRingHom& phi, const PartialRingHom& theta):
      RingHomBase(domain(theta), domain(phi)),
      myPhi(phi),
      myTheta(theta)
  {
    CoCoA_ASSERT(domain(myPhi) == codomain(myPhi));
  }


  void SequentialRingHomImpl::myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const
  {
//    RingElem tmp = myPhi(myTheta(arg));
    RingElem tmp = myPhi(myTheta(RingElemAlias(domain(myTheta),arg)));
    myDomain->myAssign(image, raw(tmp));
  }


  void SequentialRingHomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << " composite of " << myPhi << " applied to " << myTheta;
  }


  PartialRingHom sequential(const PartialRingHom& phi, const PartialRingHom& theta)
  {
    if (domain(phi) != codomain(theta))
      CoCoA_THROW_ERROR1(ERR::BadCompose);
    return PartialRingHom(new SequentialRingHomImpl(phi, theta));
  }


  bool ImageLiesInSubfield(const RingHom& phi)
  {
    if (IsField(domain(phi)) || IsField(codomain(phi))) return true;
    return codomain(phi)->myImageLiesInSubfield(phi);
  }


  //---------------------------------------------------------------------------

  void RingHomEmbeddingBase::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << " canonical embedding";
  }

  //---------------------------------------------------------------------------

  RingHomInducedBase::RingHomInducedBase(const ring& NewDomain, const RingHom& InducingHom):
      RingHomBase(NewDomain, codomain(InducingHom)),
      myInducingHom(InducingHom)
  {}


  void RingHomInducedBase::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << " induced by " << myInducingHom;
  }


} // end of namespace CoCoA
