//   Copyright (c)  2016  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingFqLog.H"
#include "CoCoA/RingFqVec.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::endl;



namespace CoCoA
{

  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    GlobalManager CoCoAFoundations;

    int p = 101;
    const int DegMinPoly = 3;
    const long CARD = SmallPower(p,DegMinPoly);

    ring FqVec = NewRingFqVec(p,DegMinPoly);

    CoCoA_ASSERT_ALWAYS(characteristic(FqVec) == p);
    CoCoA_ASSERT_ALWAYS(LogCardinality(FqVec) == DegMinPoly);
    CoCoA_ASSERT_ALWAYS(IsQuotientRing(FqVec));
    CoCoA_ASSERT_ALWAYS(IsPolyRing(BaseRing(FqVec)));
    CoCoA_ASSERT_ALWAYS(IsRingFp(CoeffRing(BaseRing(FqVec))));

    ring FqLog = NewRingFqLog(p,DegMinPoly);

    CoCoA_ASSERT_ALWAYS(characteristic(FqLog) == p);
    CoCoA_ASSERT_ALWAYS(LogCardinality(FqLog) == DegMinPoly);
    CoCoA_ASSERT_ALWAYS(IsQuotientRing(FqLog));
    CoCoA_ASSERT_ALWAYS(IsPolyRing(BaseRing(FqVec)));
    CoCoA_ASSERT_ALWAYS(IsRingFp(CoeffRing(BaseRing(FqVec))));

    {
      ring Fpx = BaseRing(FqVec);
      RingElem x = indet(Fpx,0);
      RingElem gen = QuotientingHom(FqVec)(x);
      gen = power(gen,29641);
      RingElem pwr = one(FqVec);
      for (int i=0; i < CARD-1; ++i)
      {
        CoCoA_ASSERT_ALWAYS(i == 0 || !IsOne(pwr));
        pwr *= gen;
      }
      CoCoA_ASSERT_ALWAYS(IsOne(pwr));
    }


    {
      ring FpxLog = BaseRing(FqLog);
      RingElem x = indet(FpxLog,0);
      RingElem gen = QuotientingHom(FqLog)(x);
      gen = power(gen,29641);
      RingElem pwr = one(FqLog);
      CoCoA_ASSERT_ALWAYS(IsOne(pwr));
      for (int i=0; i < CARD-1; ++i)
      {
        CoCoA_ASSERT_ALWAYS(i == 0 || !IsOne(pwr));
        pwr *= gen;
      }
      CoCoA_ASSERT_ALWAYS(IsOne(pwr));
    }

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
