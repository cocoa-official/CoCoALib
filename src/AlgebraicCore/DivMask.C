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


#include "CoCoA/DivMask.H"
#include "CoCoA/assert.H"

#include <algorithm>
using std::min;
#include <iostream>
using std::ostream;

namespace CoCoA
{

  // 2021-02-03  JAA no measurable speed difference between std::min and this impl
  // namespace
  // {
  //   // our own defn of min -- does not need addrs of its args
  //   inline long min(long a, long b) { if (a < b) return a; else return b; }
  // }

  // Need this line so that the call to std::min works (sigh)
  constexpr const long DivMask::ourMaskWidth; // value specified in DivMask.H

  std::ostream& operator<<(std::ostream& out, const DivMask& dm)
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << "DivMask(" << bits(dm) << ")";
  }


  namespace DvMskRule
  {

    // Next come the definitions of the various concrete DivMask rules: these are
    // in the .C file rather than the .H because they do not need to be in the .H
    // file, and it seems pointless cluttering the .H with unnecessary
    // implementation details.


    //-- class DivMaskNullImpl --------------------------------------------

    class NullImpl: public DivMaskRuleBase
    {
    public:
      // Default ctor and dtor are fine; nobody uses copy ctor or assignment.
      void myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const override;
      void myOutputSelf(std::ostream& out) const override;
    };


    void NullImpl::myAssignFromExpv(DivMask& dm, const SmallExponent_t* /*expv*/, long /*NumIndets*/) const
    {
      myBits(dm).reset();
    }


    void NullImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "DivMaskNull";
    }


    //-- class SingleBitImpl -----------------------------------------

    class SingleBitImpl: public DivMaskRuleBase
    {
    public:
      // Default ctor and dtor are fine; nobody uses copy ctor or assignment.
      void myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const override;
      void myOutputSelf(std::ostream& out) const override;
    };


    void SingleBitImpl::myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const
    {
      CoCoA_ASSERT(NumIndets > 0);
      myBits(dm).reset();
      const long imax = min(NumIndets, DivMask::ourMaskWidth);
      for (long i=0; i < imax; ++i)
        if (expv[i] > 0)
          myBits(dm).set(i);
    }


    void SingleBitImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "DivMaskSingleBit";
    }


    //-- class SingleBitWrapImpl -----------------------------------------

    class SingleBitWrapImpl: public DivMaskRuleBase
    {
    public:
      // Default ctor and dtor are fine; nobody uses copy ctor or assignment.
      void myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const override;
      void myOutputSelf(std::ostream& out) const override;
    };


    void SingleBitWrapImpl::myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const
    {
      CoCoA_ASSERT(NumIndets > 0);
      myBits(dm).reset();
      for (long indet=0; indet < NumIndets; ++indet)
      {
        if (expv[indet] > 0)
          // NB: i&(ourMaskWidth-1) = i%ourMaskWidth (which is a power of 2)
          myBits(dm).set(indet&(DivMask::ourMaskWidth-1));
      }
    }


    void SingleBitWrapImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "DivMaskSingleBitWrap";
    }


    //-- class EvenPowersImpl -----------------------------------------

    class EvenPowersImpl: public DivMaskRuleBase
    {
    public:
      // Default ctor and dtor are fine; nobody uses copy ctor or assignment.
      void myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const override;
      void myOutputSelf(std::ostream& out) const override;
    };


    void EvenPowersImpl::myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const
    {
      CoCoA_ASSERT(NumIndets > 0);
      myBits(dm).reset();
      const long NumBitsPerIndet = DivMask::ourMaskWidth/NumIndets;
      if (NumBitsPerIndet<2)
      {
        for (long i = 0; i < NumIndets; ++i)
          // NB: i&(ourMaskWidth-1) = i%ourMaskWidth (which is a power of 2)
          if (expv[i]) myBits(dm).set(i&(DivMask::ourMaskWidth-1));
        for (long i = NumIndets, indt = 0; i < DivMask::ourMaskWidth; ++i, ++indt )
          if (expv[indt]>1)  myBits(dm).set(i);
      }
      else
      {
        long i = 0;
        unsigned long b = 0;
        for (long indt = 0; i < NumBitsPerIndet*NumIndets; ++i, ++indt)
        {
          if (indt==NumIndets) { indt=0; ++b;}
          if (expv[indt]>2*b)  myBits(dm).set(i);
        }
        for (long indt = 0; i < DivMask::ourMaskWidth; ++i, ++indt)
          if (expv[indt]>1)  myBits(dm).set(i);
      }
    }


    void EvenPowersImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "DivMaskEvenPowers";
    }



    //-- class HashingImpl -----------------------------------------

    class HashingImpl: public DivMaskRuleBase
    {
    public:
      // Default ctor and dtor are fine; nobody uses copy ctor or assignment.
      void myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const override;
      void myOutputSelf(std::ostream& out) const override;
    private: // impl detail
      void myAdjoin(DivMask& dm, long var, SmallExponent_t exp) const;
    };


    // This impl prefers simplicity over speed.  Might get fixed later; might not.
    void HashingImpl::myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const
    {
      CoCoA_ASSERT(NumIndets > 0);
      myBits(dm).reset();
      for (long indet=0; indet < NumIndets; ++indet)
      {
        if (expv[indet] > 0)
          myAdjoin(dm, indet, expv[indet]);
      }
    }


    void HashingImpl::myAdjoin(DivMask& dm, long indet, SmallExponent_t exp) const
    {
      CoCoA_ASSERT(indet >= 0);
      CoCoA_ASSERT(exp > 0);
      const unsigned long W = DivMask::ourMaskWidth; // W is just shorthand
      unsigned long index = indet%W;
      const long step = (24*(indet/W)+13)%W; // just a heuristic
      for (unsigned long k=0; k < W && k*k < exp; ++k) // limit max value of k??
      {
        myBits(dm).set(index);
        index += step;
        if (index >= W) index -= W;
      }
    }


    void HashingImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "DivMaskHashing";
    }

  }  // end of namespace DvMskRule


  //----------------------------------------------------------------------

  // This simply prints out the name of the divmask rule.
  std::ostream& operator<<(std::ostream& out, const DivMaskRule& DMR)
  {
    if (!out) return out;  // short-cut for bad ostreams
    DMR->myOutputSelf(out);
    return out;
  }


  //----------------------------------------------------------------------
  // Here are the pseudo-constructors:

  DivMaskRule NewDivMaskNull()
  {
    return DivMaskRule(new DvMskRule::NullImpl());
  }

  DivMaskRule NewDivMaskSingleBit()
  {
    return DivMaskRule(new DvMskRule::SingleBitImpl());
  }

  DivMaskRule NewDivMaskSingleBitWrap()
  {
    return DivMaskRule(new DvMskRule::SingleBitWrapImpl());
  }

  DivMaskRule NewDivMaskEvenPowers()
  {
    return DivMaskRule(new DvMskRule::EvenPowersImpl());
  }

  DivMaskRule NewDivMaskHashing()
  {
    return DivMaskRule(new DvMskRule::HashingImpl());
  }


} // end of namespace CoCoA
