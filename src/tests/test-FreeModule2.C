//   Copyright (c)  2007-2013  John Abbott, Anna M. Bigatti

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
#include "CoCoA/FractionField.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/ModuleOrdering.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/degree.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

namespace CoCoA
{

  //----------------------------------------------------------------------
  // Test for FreeModules over SparsePolyRing with many different orderings
  // functions:  LPosn(v),  LPP(v),  wdeg(v),  IsHomog(v)
  // environments: Q[x,y,z], WDegPosnOrd/OrdPosn, shifts
  //----------------------------------------------------------------------

  void trial(const FreeModule& F)
  {
    const vector<ModuleElem>& e = gens(F);
    const SparsePolyRing P = RingOf(F);
    const vector<RingElem>& x = indets(P);

    const ModuleElem u = power(x[1],5)*e[0] + x[0]*e[1];
    const ModuleElem v = x[2]*x[2]*e[1] + x[1]*x[2]*e[2];
    const RingElem a = 2 * one(RingOf(F));  // RingOf(F) is R

    cout << "---- F = " << F << " ----" << endl;
    cout << "-- " << ordering(F) << endl;
    //  cout << "-- " << ordering(PPM(P)) << endl;
    cout << "u[0] = " << u[0] << "  and  u[1] = " << u[1] << endl;
    cout << "u = " << u << endl;
    cout << "LPosn(u) = " << LPosn(u) << std::endl;
    cout << "LPP(u) = " << LPP(u) << endl;
    cout << "wdeg(u) = " << wdeg(u) << endl;
    if (GradingDim(P)==0)
      cout << "IsHomog(u) undefined because GradingDim is 0" << endl;
    else
      cout << "IsHomog(u) = " << IsHomog(u) << std::endl;
    cout << "v = " << v << endl;
    cout << "LPosn(v) = " << LPosn(v) << std::endl;

    cout << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;
    std::cout << std::boolalpha;

    ring QQ = RingQQ();
    SparsePolyRing PLex = NewPolyRing(QQ, symbols("x,y,z"), lex);
    SparsePolyRing PDegLex = NewPolyRing(QQ, symbols("x,y,z"), StdDegLex);
    long n = 4;
    std::vector<degree> sh(n, wdeg(one(PDegLex)));
    sh[1] = wdeg(power(indet(PDegLex,0),4));
    trial(NewFreeModule(PLex, n));
    trial(NewFreeModule(PDegLex, n));
    //  trial(NewFreeModule(PDegLex, n, PosnOrd));
    trial(NewFreeModule(PDegLex, n, OrdPosn));
    trial(NewFreeModule(PLex, n, WDegPosnOrd));
    trial(NewFreeModule(PDegLex, n, WDegPosnOrd));
    trial(NewFreeModule(PDegLex, sh, WDegPosnOrd));
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
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
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
