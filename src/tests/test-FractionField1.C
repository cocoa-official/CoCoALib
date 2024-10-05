//   Copyright (c)  2013  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/BigRat.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"


#include <iostream>
using std::cerr;
using std::endl;



namespace CoCoA
{

  // Test for FractionField: ctor, basic arith, InducedHom
  void program()
  {
    GlobalManager CoCoAFoundations;

    ring ZZ = RingZZ();
    FractionField QQ = RingQQ();
    ring Q = NewFractionField(ZZ);
    CoCoA_ASSERT_ALWAYS(Q == QQ);
    RingElem q(QQ,2);
    q /= 3;
    CoCoA_ASSERT_ALWAYS(num(q) == 2 && den(q) == 3);

    ring QQx = NewPolyRing(QQ, symbols("x"));
    FractionField FrFQQx = NewFractionField(QQx);
    RingElem x(FrFQQx, symbol("x"));
    RingElem f = (x-1)/(x+1);
    RingHom embed = EmbeddingHom(FrFQQx);
    CoCoA_ASSERT_ALWAYS(embed(num(f)) == x-1);
    CoCoA_ASSERT_ALWAYS(embed(den(f)) == x+1);

    RingElem f1 = deriv(f, x);
    CoCoA_ASSERT_ALWAYS(f1 == 2/power(x+1,2));

    RingHom ZZtoQQx = CanonicalHom(ZZ, QQx);
    RingHom QQtoQQx = InducedHom(QQ, ZZtoQQx);
    RingElem image = QQtoQQx(q);
    CoCoA_ASSERT_ALWAYS(image == BigRat(2,3));
    RingHom QQtoFrF = InducedHom(QQ, CanonicalHom(ZZ, FrFQQx));
    CoCoA_ASSERT_ALWAYS(QQtoFrF(q) == BigRat(2,3));
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
