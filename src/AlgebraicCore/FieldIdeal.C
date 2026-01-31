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
    for (const auto& g: gens)
      if (!IsZero(g)) { AllZero = false; break; }
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
    CoCoA_THROW_ERROR1(ERR::NYI);
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
      CoCoA_THROW_ERROR1(ERR::ReqField);
    for (const auto& g: gens)
      if (owner(g) != k)
        CoCoA_THROW_ERROR1(ERR::MixedRings);
    return ideal(new FieldIdealImpl(k, gens));
  }

} // end of namespace CoCoA
