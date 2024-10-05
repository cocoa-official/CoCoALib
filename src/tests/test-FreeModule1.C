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
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/degree.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;
namespace CoCoA
{
  //----------------------------------------------------------------------
  // Test for GradedFreeModules in the basic case
  // functions:  LPosn(v),  LPP(v),  wdeg(v),  IsHomog(v)
  // environments: Q[x,y,z]
  //----------------------------------------------------------------------


  void trial(const ring& R)
  {
    const long dim = 4;
    const FreeModule F = NewFreeModule(R, dim);
    const vector<ModuleElem>& e = gens(F);
    const vector<RingElem>& x = indets(R);
    CoCoA_ASSERT_ALWAYS(len(e) == dim);
// This does not yet compile
//   for (long i=0; i < dim; ++i)
//     for (long j=0; j < dim; ++j)
//       CoCoA_ASSERT_ALWAYS(Fgens[i][j] == (i==j));
    CoCoA_ASSERT_ALWAYS(IsFGModule(F));


    ModuleElem v1 = (x[0]-x[1])*e[0] + (x[0]*x[0]-x[1]*x[2])*e[1]
      + x[0]*x[0]*e[2] + (x[1]*x[2])*e[3];
    ModuleElem v2 = 4*e[0] + 3*e[1] + 2*e[2] + 1*e[3];
    //  cout << "v1 = " << v1 << std::endl;
    CoCoA_ASSERT_ALWAYS( LPosn(v1) == 1 );
    CoCoA_ASSERT_ALWAYS( LPP(v1) == LPP(x[0]*x[0]) );
    CoCoA_ASSERT_ALWAYS( wdeg(v1)[0] == 2 );
    CoCoA_ASSERT_ALWAYS( !IsHomog(v1) );
    CoCoA_ASSERT_ALWAYS( IsHomog(v2) );
    
    //  CmpWDeg(v1, v2);
    //  CmpLT(v1, v2);
    //  CmpPos(v1, v2);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    trial(NewPolyRing(RingQQ(), symbols("x,y,z")));
    trial(NewPolyRing(NewZZmod(101), SymbolRange("x",1,9)));
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
