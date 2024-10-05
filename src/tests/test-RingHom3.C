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
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ring.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for RingHom
// ZZxy->FpxII, FpxII->Fpy, FpxII->FpxII
//----------------------------------------------------------------------
namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    const QuotientRing Fp = NewZZmod(101);
    const PolyRing ZZxy = NewPolyRing(RingZZ(), symbols("x,y"));
    const PolyRing Fpy = NewPolyRing(Fp, SymbolRange("y",0,1));
    const PolyRing FpxI = NewPolyRing_DMPI(Fp, SymbolRange("z",0,1));
    const PolyRing FpxII = NewPolyRing_DMPII(Fp, SymbolRange("x",0,1));

//   RingHom ZZToFpxI = ZZEmbeddingHom(FpxI);
//   RingHom ZZToFpxII = ZZEmbeddingHom(FpxII);
    RingHom ZZToFpxI = CanonicalHom(RingZZ(), FpxI);
    RingHom ZZToFpxII = CanonicalHom(RingZZ(), FpxII);

    const vector<RingElem>& x = indets(FpxII);
    vector<RingElem> imx;
    imx.push_back(x[0]*x[0]);
    imx.push_back(x[1]*x[1]);

    const vector<RingElem>& z = indets(FpxI);
    vector<RingElem> imx1;
    imx1.push_back(z[0]*z[0]);
    imx1.push_back(z[1]*z[1]);

    const vector<RingElem>& y = indets(Fpy);
    vector<RingElem> imx2;
    imx2.push_back(y[0]*y[0]);
    imx2.push_back(y[1]*y[1]);

    RingHom ZZxyToFpxII = PolyRingHom(ZZxy, FpxII, ZZToFpxII, imx);
    RingHom FpxIIToFpxII = PolyRingHom(FpxII, FpxII, CoeffEmbeddingHom(FpxII), imx);
    RingHom FpxIIToFpy = PolyRingHom(FpxII, Fpy, CoeffEmbeddingHom(Fpy), imx2);

    RingHom ZZxyToFpxI = PolyRingHom(ZZxy, FpxI, ZZToFpxI, imx1);
//   RingHom FpxIToFpxI = PolyRingHom(FpxI, FpxI, CoeffEmbeddingHom(FpxI), imx1);
//   RingHom FpxIToFpy = PolyRingHom(FpxI, Fpy, CoeffEmbeddingHom(Fpy), imx2);
    RingHom FpxIToFpxI = PolyAlgebraHom(FpxI, FpxI, imx1);
    RingHom FpxIToFpy = PolyAlgebraHom(FpxI, Fpy, imx2);

    cout << "Simple test on polynomial ring hom" << endl
         << ZZxy << "  -->  " << FpxII << endl
         << "sending " << indet(ZZxy,0) << "  to  " << imx[0] << "  and" << endl
         << "sending " << indet(ZZxy,1) << "  to  " << imx[1] << endl;

    RingElem f = indet(ZZxy,0) + 2*indet(ZZxy,1);
    RingElem g = ZZxyToFpxII(f);
    RingElem h = ZZxyToFpxI(f);

    cout << f << " in " << owner(f) << "\n  maps to " << g << " in " << owner(g) << endl << endl;
    cout << g << " in " << owner(g) << "\n  maps to " << FpxIIToFpxII(g) << " in " << owner(FpxIIToFpxII(g)) << endl << endl;
    cout << g << " in " << owner(g) << "\n  maps to " << FpxIIToFpy(g) << " in " << owner(FpxIIToFpy(g)) << endl << endl;

    cout << f << " in " << owner(f) << "\n  maps to " << h << " in " << owner(h) << endl << endl;
    cout << h << " in " << owner(h) << "\n  maps to " << FpxIToFpxI(h) << " in " << owner(FpxIToFpxI(h)) << endl << endl;
    cout << h << " in " << owner(h) << "\n  maps to " << FpxIToFpy(h) << " in " << owner(FpxIToFpy(h)) << endl << endl;

    RingHom EvalFpxI = PolyAlgebraHom(FpxI, Fp, "13,31");
    RingHom EvalFpxII = PolyAlgebraHom(FpxII, Fp, "13,31");
    CoCoA_ASSERT_ALWAYS(EvalFpxII(g) == EvalFpxI(h));
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
