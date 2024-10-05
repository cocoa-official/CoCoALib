//   Copyright (c)  2009,2021  John Abbott &  Anna M. Bigatti

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

#include "CoCoA/SugarDegree.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/ReductionCog.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/TmpGPoly.H"
#include "CoCoA/assert.H"

#include <algorithm>
using std::max; // for std sugar
#include <iostream>
//using std::ostream;
//#include <memory> // already included by SugarDegree.H
using std::unique_ptr;


namespace CoCoA
{
  std::ostream& operator<<(std::ostream& out, const SugarDegree& s)
  {
    if (!out) return out;  // short-cut for bad ostreams
    if (IsInitialized(s)) s->myOutput(out);
    else out << "uninitialized";
    return out;
  }


  SugarDegree::SugarDegree(const SugarDegree& rhs)
  {
    if (IsInitialized(rhs)) myPtr.reset((rhs.myPtr)->myClone());
  }


  SugarDegree& SugarDegree::operator=(const SugarDegree& rhs)
  {
    if (this == &rhs) return *this;
    if (IsInitialized(rhs)) myPtr.reset((rhs.myPtr)->myClone());
    else                    myPtr.reset(nullptr);
    return *this;
  }
  

  void SugarDegreeBase::myUpdate(ReductionCog F, const GPoly& g)
  { myUpdate(ActiveLPP(F)/LPPForOrd(g), g); }
  

  namespace SugarTypes
  {
    //---  StdDegBase and WDegBase abstract classes -------------

    class StdDegBase: public SugarDegreeBase
    {
    public:
      StdDegBase(long s): SugarDegreeBase() { myValue = s; }

      // functions every SugarDegree must implement
      const degree& myWSugar() const override;  ///< returns error
      long myStdSugar() const override{ return myValue; }
      int myCmp(const SugarDegreeBase& s) const override;  // this <=> s ? <0,=0,>0

      void myOutput(std::ostream& out) const override;
      void myMul(ConstRefPPMonoidElem pp) override;
      void myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g) override;
      virtual long myDeg(ConstRefPPMonoidElem pp) const = 0;
      long myDeg0(ConstRefRingElem f) const; // based on myDeg; also defined on 0
    protected:
      long myValue;
    };


    class WDegBase: public SugarDegreeBase
    {
    public:
      WDegBase(const degree& s): SugarDegreeBase(), myValue(s) {}

      // functions every SugarDegree must implement
      const degree& myWSugar() const override { return myValue; }
      long myStdSugar() const override;   ///< returns error
      int myCmp(const SugarDegreeBase& s) const override;  // this <=> s ? <0,=0,>0

      void myOutput(std::ostream& out) const override;
      void myMul(ConstRefPPMonoidElem pp) override;
      void myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g) override;
      virtual degree myDeg(ConstRefPPMonoidElem pp) const = 0;
      degree myDeg0(ConstRefRingElem f) const; // based on myDeg; also defined on 0
    protected:
      degree myValue;
    };


    //---  StdDegBase concrete classes -------------

    class StdDegPlain: public StdDegBase
    {
    public:
      StdDegPlain(long s): StdDegBase(s) {}
      StdDegPlain(ConstRefRingElem f): StdDegBase(0) { myValue=myDeg0(f); }
      StdDegPlain* myClone() const override  { return new StdDegPlain(myValue); }
      long myDeg(ConstRefPPMonoidElem pp) const override;
    };


    class StdDegNoIdx: public StdDegBase
    {
    public:
      StdDegNoIdx(long s, long NoIdx): StdDegBase(s) { myNoIdxValue=NoIdx;}
      StdDegNoIdx(ConstRefRingElem f, long NoIdx): StdDegBase(0) { myNoIdxValue=NoIdx; myValue=myDeg0(f); }
      StdDegNoIdx* myClone() const override  { return new StdDegNoIdx(myValue, myNoIdxValue); }
      long myDeg(ConstRefPPMonoidElem pp) const override;
    private:
      long myNoIdxValue;
    };


    class StdDegNoIdxSat: public StdDegBase
    {
    public:
      StdDegNoIdxSat(long s, long NoIdx): StdDegBase(s) { myNoIdxValue=NoIdx;}
      StdDegNoIdxSat(ConstRefRingElem f, long NoIdx): StdDegBase(0) { myNoIdxValue=NoIdx; myValue=myDeg0(f); }
      StdDegNoIdxSat* myClone() const override  { return new StdDegNoIdxSat(myValue, myNoIdxValue); }
      long myDeg(ConstRefPPMonoidElem pp) const override;
    private:
      long myNoIdxValue;
    };


    class StdDegSat: public StdDegBase
    {
    public:
      StdDegSat(long s): StdDegBase(s) {}
      StdDegSat(ConstRefRingElem f): StdDegBase(0) { myValue=myDeg0(f); }
      StdDegSat* myClone() const override  { return new StdDegSat(myValue); }
      long myDeg(ConstRefPPMonoidElem pp) const override;
    };


    //---  WDegBase concrete classes -------------

    class WDegPlain: public WDegBase
    {
    public:
      WDegPlain(const degree& s): WDegBase(s) {}
      WDegPlain(ConstRefRingElem f);
      WDegPlain* myClone() const override  { return new WDegPlain(myValue); }
      degree myDeg(ConstRefPPMonoidElem pp) const override;
    };


    class WDeg1CompTmp: public WDegBase
    {
    public:
      WDeg1CompTmp(const degree& s): WDegBase(s) {}
      WDeg1CompTmp(ConstRefRingElem f);
      WDeg1CompTmp* myClone() const override  { return new WDeg1CompTmp(myValue); }
      degree myDeg(ConstRefPPMonoidElem pp) const override;
      void myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g) override; ///< no default
    };


    class WDegConst: public WDegBase
    {
    public:
      WDegConst(const degree& s): WDegBase(s) {}
      WDegConst(ConstRefRingElem f);
      WDegConst* myClone() const override  { return new WDegConst(myValue); }
      degree myDeg(ConstRefPPMonoidElem pp) const override;
      void myUpdate(ConstRefPPMonoidElem /*CofactorPP*/, const GPoly& /*g*/) override {}; ///< no default
    };


    class WDegSat: public WDegBase ///< wip
    {
    public:
      WDegSat(const degree& s): WDegBase(s) {}
      WDegSat(ConstRefRingElem f);
      WDegSat* myClone() const override  { return new WDegSat(myValue); }
      degree myDeg(ConstRefPPMonoidElem pp) const override;
      void myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g) override;
    };



    //--  StdDegBase  ----------------------------------

    const degree& StdDegBase::myWSugar() const
    {
      CoCoA_THROW_ERROR("Wrong type for this sugar", "SugarTypes::StdDegBase");
      static degree NeverUsed(0); // just to keep the compiler quiet
      return NeverUsed;
    }

    void StdDegBase::myOutput(std::ostream& out) const
    { out << myValue; }

    void StdDegBase::myMul(ConstRefPPMonoidElem pp)
    { myValue += myDeg(pp); }

    void StdDegBase::myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g)
    { myValue = max(myValue, myDeg(CofactorPP) + sugar(g)->myStdSugar()); }


    int StdDegBase::myCmp(const SugarDegreeBase& s) const
    {
      CoCoA_ASSERT((dynamic_cast<const StdDegBase*>(&s)) != nullptr);
      const StdDegBase& ss = static_cast<const StdDegBase&>(s);
      if (myValue != ss.myValue)
      {
        if (myValue > ss.myValue) return 1; else return -1;
      }
      return 0;
    }

    long StdDegBase::myDeg0(ConstRefRingElem f) const
    {
      if (IsZero(f)) return 0;
      return myDeg(LPP(f));  // ANNA: wrong! the right one follows
      // after finishing the testing phase use the correct definition
      long degf = 0;
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
        degf = max(degf, myDeg(PP(it)));
      return degf;
    }
    
    //--  WDegBase  ----------------------------------

    long WDegBase::myStdSugar() const
    {
      CoCoA_THROW_ERROR("Wrong type for this sugar", "SugarTypes::WDegBase");
      return 0;  // just to keep the compiler quiet
    }

    void WDegBase::myOutput(std::ostream& out) const
    { out << myValue; }

    void WDegBase::myMul(ConstRefPPMonoidElem pp)
    { myValue += myDeg(pp); }

    void WDegBase::myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g)
    { myValue = top(myValue, myDeg(CofactorPP) + sugar(g)->myWSugar()); }

    int WDegBase::myCmp(const SugarDegreeBase& s) const
    {
      CoCoA_ASSERT(dynamic_cast<const WDegBase*>(&s) != nullptr);
      const WDegBase& ws = static_cast<const WDegBase&>(s);
      return FastCmp(myValue, ws.myValue); // cmp = compatibility test + FastCmp
    }

    degree WDegBase::myDeg0(ConstRefRingElem f) const
    {
      if (IsZero(f)) return myDeg(one(PPM(owner(f))));
      return myDeg(LPP(f));
      // ANNA: BUG BUG BUG return above is wrong!!!! the right one follows
      SparsePolyIter it = BeginIter(f);
      degree degf = myDeg(PP(it));
      for (++it; !IsEnded(it); ++it)
        degf = top(degf, myDeg(PP(it)));
      return degf;
    }

    //--  end abstract classes  ----------------------------------

    //--  StdDegPlain  ----------------------------------

    long StdDegPlain::myDeg(ConstRefPPMonoidElem pp) const
    { return StdDeg(pp); }

    //--  StdDegSat  ----------------------------------

    long StdDegSat::myDeg(ConstRefPPMonoidElem pp) const
    { return StdDeg(pp) - exponent(pp, NumIndets(owner(pp))-1); }

    //--  StdDegNoIdx  ----------------------------------

    long StdDegNoIdx::myDeg(ConstRefPPMonoidElem pp) const
    { return StdDeg(pp) - exponent(pp, myNoIdxValue); }

    //--  StdDegNoIdxSat  ----------------------------------

    long StdDegNoIdxSat::myDeg(ConstRefPPMonoidElem pp) const
    {
      return StdDeg(pp)
        - exponent(pp, myNoIdxValue)
        - exponent(pp, NumIndets(owner(pp))-1);
    }

    //--  WDegPlain  ----------------------------------

    WDegPlain::WDegPlain(ConstRefRingElem f): WDegBase(degree(GradingDim(owner(f)))) // initialize
    { myValue = myDeg0(f); }

    degree WDegPlain::myDeg(ConstRefPPMonoidElem pp) const
    { return wdeg(pp); }

    //--  WDegConst  ----------------------------------

    WDegConst::WDegConst(ConstRefRingElem f): WDegBase(degree(GradingDim(owner(f))))
    { if (!IsZero(f)) myValue=wdeg(LPP(f)); }

    degree WDegConst::myDeg(ConstRefPPMonoidElem pp) const
    { return wdeg(pp); }

    // inline myUpdate(pp, g) does *nothing* (unlike default implementation)

    //--  WDeg1CompTmp  ----------------------------------

    WDeg1CompTmp::WDeg1CompTmp(ConstRefRingElem f): WDegBase(degree(GradingDim(owner(f))))
    { myValue = myDeg0(f); }

    degree WDeg1CompTmp::myDeg(ConstRefPPMonoidElem pp) const
    { return wdeg(pp); }

    void WDeg1CompTmp::myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g)
    { // very temporary!!  only for checking backward compatibility
      // moreover it uses max instead of top (a>b ? a : b)
      // could it be useful with only one component?
      if (g.myGRingInfo().IsMyGradingPosPlus())
        CoCoA_THROW_ERROR("TMP: use StdDeg for PosTO", "WDeg1CompTmp::myUpdate");
      myValue = max(myValue, wdeg(CofactorPP) + sugar(g)->myWSugar());
    }

    //--  WDegSat  ----------------------------------

    WDegSat::WDegSat(ConstRefRingElem f): WDegBase(degree(GradingDim(owner(f))))
    { myValue = myDeg0(f); }


    degree WDegSat::myDeg(ConstRefPPMonoidElem pp) const
    {
      degree tmp = wdeg(pp);
      // bi-grading trick
      // tmp.mySetComponent(0, tmp[1]);//WARNING: this works only for GrDim=1
      // TEMPORARY
      const BigInt d = tmp[0] - exponent(pp, NumIndets(owner(pp))-1); // only for GrDim=1
      SetComponent(tmp, 0, d);
      return tmp;
    }


    void WDegSat::myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g)
    {
      if (g.myGRingInfo().IsMyGradingPosPlus())
        CoCoA_THROW_ERROR("TMP: use StdDeg for PosTO", "WDeg1CompTmp::myUpdate");
      myValue = max(myValue, myDeg(CofactorPP) + sugar(g)->myWSugar());
    }

    //-------

  }  // namespace SugarTypes

  //----------------------------------------------------------------------
  // Here are the pseudo-constructors:

  SugarDegree NewStdSugar(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::StdDegPlain(f)); }

  SugarDegree NewStdSugarNoIdx(ConstRefRingElem f, long NoIdx)
  { return SugarDegree(new SugarTypes::StdDegNoIdx(f, NoIdx)); }

  SugarDegree NewStdSugarSat(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::StdDegSat(f)); }

  SugarDegree NewStdSugarNoIdxSat(ConstRefRingElem f, long NoIdx)
  { return SugarDegree(new SugarTypes::StdDegNoIdxSat(f, NoIdx)); }


  SugarDegree NewWSugar(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::WDegPlain(f)); }

  SugarDegree NewWDeg1CompTmp(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::WDeg1CompTmp(f)); }

  SugarDegree NewWSugarConst(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::WDegConst(f)); }

  SugarDegree NewWSugarSat(ConstRefRingElem f)
  { return SugarDegree(new SugarTypes::WDegSat(f)); }



} // namespace CoCoA



//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SugarDegree.C,v 1.29 2022/02/18 14:11:59 abbott Exp $
// $Log: SugarDegree.C,v $
// Revision 1.29  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.28  2021/10/29 19:47:59  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.27  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.26  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.25  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.24  2019/10/15 11:54:09  abbott
// Summary: Changed 0 into nullptr (where appropriate)
//
// Revision 1.23  2019/03/04 10:28:00  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.22  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.21  2018/05/18 12:54:02  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.20  2018/05/17 15:51:48  bigatti
// -- added include SparsePolyIter
//
// Revision 1.19  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.18  2014/07/07 12:45:30  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.17  2014/05/09 09:29:03  bigatti
// added  using SugarDegreeBase::myUpdate; // disables warnings of overloading
//
// Revision 1.16  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.15  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.14  2010/05/21 16:11:37  bigatti
// -- removed CheckSugar
//
// Revision 1.13  2010/05/14 09:53:09  bigatti
// -- removed empty ctor for SugarDegree
// -- added marker for SugarDegree(uninitialized)
// -- SugarDegree for GBasis input is initialized by myPrepareGBasis
//
// Revision 1.12  2009/10/27 17:15:55  bigatti
// -- fixed: using sugar(g)->myWSugar() insted of wsugar(g)
//
// Revision 1.11  2009/09/25 12:36:47  bigatti
// -- added CheckSugar (temporary)
//
// Revision 1.10  2009/07/24 15:13:50  bigatti
// -- some fixes
//
// Revision 1.9  2009/07/20 14:28:44  bigatti
// -- modified interface for pseudo constructors (now with RingElem argument)
//
// Revision 1.8  2009/04/24 16:37:41  bigatti
// -- added myDeg and pseudo-ctor
// -- improved inheritance
//
// Revision 1.7  2009/03/20 14:01:34  bigatti
// -- minor changes (names)
//
// Revision 1.6  2009/03/18 16:37:13  bigatti
// -- "almost-final" cleanup, tested against some current sugar definitions
//
// Revision 1.5  2009/03/18 15:12:02  abbott
// Minor cosmetic changes.
//
// Revision 1.4  2009/03/16 16:39:48  bigatti
// -- now use auto_ptr instead of SmartPtrIRC
// -- added copy ctor, assignment, cmp
//
// Revision 1.3  2009/03/16 07:27:20  bigatti
// -- added necessary "const"
// -- added WDeg1CompTmp (temporary, for testing)
//
// Revision 1.2  2009/02/20 15:45:25  bigatti
// -- sugar.H --> SugarDegree.H
//
// Revision 1.1  2009/02/20 13:27:15  bigatti
// -- renamed from "sugar.[CH]"
//
// Revision 1.3  2009/02/20 11:01:15  bigatti
// -- added NewHomogWSugar for graded=homoeneous case (just constant)
//
// Revision 1.2  2009/02/20 09:53:27  bigatti
// -- changed name: sweetener --> SugarDegree
// -- introduced abstract classes SugarTypes::StdDegBase, SugarTypes::WDegBase
//
// Revision 1.1  2009/02/09 13:57:05  bigatti
// -- first import
//
