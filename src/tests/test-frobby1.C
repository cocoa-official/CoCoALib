//   Copyright (c)  2009 Anna Bigatti

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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"

#ifdef CoCoA_WITH_FROBBY
#include "CoCoA/ExternalLibs-Frobby.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ideal.H"
#endif


#include <iostream>
using std::cerr;
using std::endl;
#ifdef CoCoA_WITH_FROBBY
#include <vector>
using std::vector;
#endif
//----------------------------------------------------------------------
// First test for Frobby library.
// Trivial computation to test proper integration.
//----------------------------------------------------------------------
namespace CoCoA
{

void program()
{
  GlobalManager CoCoAFoundations;

#ifdef CoCoA_WITH_FROBBY
  PolyRing P = NewPolyRing(RingQQ(), symbols("x,y,z"));
  RingElem x(indet(P,0));
  RingElem y(indet(P,1));
  RingElem z(indet(P,2));
  
  //I := Ideal(x^2, x*y, y^2, z^2);
  ideal I = ideal(x*x, x*y, y*y, z*z);
  
  ideal AD = FrbAlexanderDual(I, LPP(x*x*y*y*z*z));// = Ideal(x^2yz, xy^2z);
  CoCoA_ASSERT_ALWAYS(AD == ideal(x*x*y*z, x*y*y*z));

  ideal MSM = FrbMaximalStandardMonomials(I);// = Ideal(yz, xz);
  CoCoA_ASSERT_ALWAYS(MSM == ideal(y*z, x*z));

  vector<ideal> ID;
  FrbIrreducibleDecomposition(ID, I);// = [Ideal(x, y^2, z^2), Ideal(x^2, y, z^2)];
  CoCoA_ASSERT_ALWAYS(ID.size() == 2);
  CoCoA_ASSERT_ALWAYS(ID[0] == ideal(x, y*y, z*z));
  CoCoA_ASSERT_ALWAYS(ID[1] == ideal(x*x, y, z*z));

  vector<ideal> PD;
  FrbPrimaryDecomposition(PD, I);
  CoCoA_ASSERT_ALWAYS(PD.size() == 1);
  CoCoA_ASSERT_ALWAYS(PD[0] == I);

  vector<ideal> AP;
  FrbAssociatedPrimes(AP, I);
  CoCoA_ASSERT_ALWAYS(AP.size() == 1);
  CoCoA_ASSERT_ALWAYS(AP[0] == ideal(x, y, z));

  CoCoA_ASSERT_ALWAYS(FrbDimension(ideal(x, y)) == 1);
  CoCoA_ASSERT_ALWAYS(FrbDimension(ideal(zero(P))) == 3);
  CoCoA_ASSERT_ALWAYS(FrbDimension(ideal(one(P))) == -1);

  CoCoA_ASSERT_ALWAYS(FrbMultigradedHilbertPoincareNumerator(ideal(x,y)) ==
			  1 - x - y + x * y);
  CoCoA_ASSERT_ALWAYS(FrbMultigradedHilbertPoincareNumerator(ideal(zero(P))) == 1);
  CoCoA_ASSERT_ALWAYS(FrbMultigradedHilbertPoincareNumerator(ideal(one(P))) == 0);

  CoCoA_ASSERT_ALWAYS(FrbTotalDegreeHilbertPoincareNumerator(ideal(x,y), x + 2 * y) ==
			  1 - 2 * (x + 2 * y) + (x + 2 * y) * (x + 2 * y));
  CoCoA_ASSERT_ALWAYS(FrbTotalDegreeHilbertPoincareNumerator(ideal(zero(P)), x) == 1);
  CoCoA_ASSERT_ALWAYS(FrbTotalDegreeHilbertPoincareNumerator(ideal(one(P)), y) == 0);

  const RingElem p = FrbTotalDegreeHilbertPoincareNumerator(ideal(x,y));
  const RingElem pIndet = indets(owner(p))[0];
  CoCoA_ASSERT_ALWAYS(p == 1 - 2 * pIndet + pIndet * pIndet);
  CoCoA_ASSERT_ALWAYS(FrbTotalDegreeHilbertPoincareNumerator(ideal(zero(P))) == 1);
  CoCoA_ASSERT_ALWAYS(FrbTotalDegreeHilbertPoincareNumerator(ideal(one(P))) == 0);

#endif
}

} // end of namespace CoCoA


//----------------------------------------------------------------------
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
