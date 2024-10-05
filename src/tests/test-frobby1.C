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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-frobby1.C,v 1.17 2018/09/28 15:54:06 abbott Exp $
// $Log: test-frobby1.C,v $
// Revision 1.17  2018/09/28 15:54:06  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.16  2018/05/18 16:46:17  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.15  2016/09/16 16:36:41  abbott
// Summary: Changed TEST_ASSERT into CoCoA_ASSERT_ALWAYS; removed all include assert.H lines
//
// Revision 1.14  2015/05/07 15:00:40  abbott
// Summary: Moved test code into namespace CoCoA
// Author: JAA
//
// Revision 1.13  2014/07/07 13:39:55  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.12  2013/06/28 12:07:48  abbott
// Updated Frobby fn names.
//
// Revision 1.11  2012/02/10 11:57:12  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.10  2012/02/08 17:36:48  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.9  2011/07/29 14:58:58  bigatti
// -- added (temporarily?) "Frobby" suffix to Frobby functions
//
// Revision 1.8  2011/07/05 15:15:33  bigatti
// -- changed AlexanderDual --> AlexanderDualFrobby
// -- changed PrimaryDecomposition --> PrimaryDecompositionFrobby
//
// Revision 1.7  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.6  2010/11/26 15:32:16  bigatti
// -- renamed Dimension into dimension (coding conventions)
//
// Revision 1.5  2010/11/22 17:39:15  bigatti
// -- updated "TmpFrobby" --> "ExternalLibs-Frobby"
//
// Revision 1.4  2010/02/04 09:38:28  bigatti
// -- new syntax for frobby (more CoCoALib-like)
//
// Revision 1.3  2010/02/03 18:00:22  bigatti
// -- more functions from frobby
//
// Revision 1.2  2009/07/30 15:41:25  bigatti
// -- now using the new nice constructors for ideals
//
// Revision 1.1  2009/02/11 15:08:26  bigatti
// -- first import
//

