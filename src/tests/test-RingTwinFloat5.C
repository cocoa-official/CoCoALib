//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/BuildInfo.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"


#include <iostream>
using std::cerr;
using std::endl;



namespace CoCoA
{

// This routine deliberately amplifies any deviation from 1 in r
  void IsTrulyOne(RingElem r)
  {
    RingElem one(owner(r),1); one /= one; /*DELIBERATE SELF ASSIGN*/
    const RingElem nine = one+one+one+one+one+one+one+one+one;
    const RingElem ten = nine+one;

    for(int i=0; i <= 99; ++i)
      r = ten*r-nine;
    CoCoA_ASSERT_ALWAYS(r == 1);
  }


  void trial(const ring& R)
  {
    // We cannot conduct this test if 3 is not invertible
    if (!IsInvertible(RingElem(R,3))) return;

    RingElem one(R,1); one /= one; /*DELIBERATE SELF ASSIGN*/
    const RingElem three = one+one+one;
    try
    {
      IsTrulyOne(RingElem(R,1));
    }
    catch (const RingTwinFloat::InsufficientPrecision&) {}

    try
    {
      IsTrulyOne(one);
    }
    catch (const RingTwinFloat::InsufficientPrecision&) {}

    try
    {
      IsTrulyOne(3*(one/3));
    }
    catch (const RingTwinFloat::InsufficientPrecision&) {}

    try
    {
      IsTrulyOne(three*(one/three));
    }
    catch (const RingTwinFloat::InsufficientPrecision&) {}
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    for (int prec=10; prec < 333; ++prec)
      trial(NewRingTwinFloat(prec));

    trial(RingQQ());
    trial(NewZZmod(13));
  }

} // end of namespace CoCoA


// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
