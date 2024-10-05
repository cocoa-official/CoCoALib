//   Copyright (c)  2014  John Abbott,  Anna Bigatti
//   Orig author  2012-2014  Christof Soeger

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
#include "CoCoA/BigRatOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/ExternalLibs-Normaliz.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuasiPoly.H"
#include "CoCoA/SparsePolyOps-hilbert.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


// First test for Normaliz library.
// Simple computation to test proper integration.
// This tests is for the functions working directly on cone.
//----------------------------------------------------------------------

namespace CoCoA
{

  vector<BigInt> BigIntVec(const int* CVector, int len)
  {
    vector<BigInt> v(len);
    for (int i=0; i<len; ++i)  v[i] = (BigInt(CVector[i]));
    return v;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

#ifdef CoCoA_WITH_NORMALIZ
    using namespace CoCoA::Normaliz;

// this is the polytope example from the Normaliz examples
    const int M[4][4] = {{0, 0, 0, 1},
                         {2, 0, 0, 1},
                         {0, 3, 0, 1},
                         {0, 0, 5, 1}};
    const int g[4] = {0, 0, 0, 1};
    vector<vector<BigInt> > l;
    for (int i=0; i<4; ++i)
      l.push_back(BigIntVec(M[i], 4));

//  cout << "l -> " << len(l) << endl;

    std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > m;
    m[libnormaliz::Type::integral_closure] = l;
    m[libnormaliz::Type::grading] = vector<vector<BigInt> >(1,BigIntVec(g, 4));

    vector<vector<BigInt> > res_m;
    vector<BigInt>  res_v;

    cone C(m);
    vector<vector<BigInt> > hb = HilbertBasis(C);

    CoCoA_ASSERT_ALWAYS(len(hb) == 19);
    for (long i=0; i<len(hb); ++i)
    {
      CoCoA_ASSERT_ALWAYS(len(hb[i]) == 4);
      //check if the points are inside the cone by evaluating the support hyperplanes
      CoCoA_ASSERT_ALWAYS(hb[i][0] >= 0);
      CoCoA_ASSERT_ALWAYS(hb[i][1] >= 0);
      CoCoA_ASSERT_ALWAYS(hb[i][2] >= 0);
      CoCoA_ASSERT_ALWAYS(30*hb[i][3]-15*hb[i][0]-10*hb[i][1]-6*hb[i][2] >= 0);
    }

    const vector< vector<BigInt> > deg1 = Deg1Elements(C);
    CoCoA_ASSERT_ALWAYS(len(deg1) == 18);

    const vector< vector<BigInt> > sh = SupportHyperplanes(C);
    CoCoA_ASSERT_ALWAYS(len(sh) == 4);

    const vector< vector<BigInt> > rays = ExtremeRays(C);
    CoCoA_ASSERT_ALWAYS(len(rays) == 4);

    const vector< vector<BigInt> > equ = Equations(C);
    CoCoA_ASSERT_ALWAYS(len(equ) == 0);

    const vector< vector<BigInt> > cong = Congruences(C);
    CoCoA_ASSERT_ALWAYS(len(cong) == 0);

    CoCoA_ASSERT_ALWAYS(IsPointed(C));
    CoCoA_ASSERT_ALWAYS(!IsInhomogeneous(C));
    CoCoA_ASSERT_ALWAYS(!IsIntegrallyClosed(C));
    CoCoA_ASSERT_ALWAYS(!IsDeg1HilbertBasis(C));

    CoCoA_ASSERT_ALWAYS(multiplicity(C) == 30);

    HPSeries hs = HilbertSeries(C);
    // compare it with:  (1 + 14*t + 15*t^2) / (1-t)^4
    PolyRing QQt = RingQQt(1);
    const RingElem t = indet(QQt,0);
    RingElem ref_num(RingElem(QQt, "1 + 14*t + 15*t^2"));
    CoCoA_ASSERT_ALWAYS(num(hs) == ref_num);
    factorization<RingElem> den = DenFactors(hs);
    CoCoA_ASSERT_ALWAYS(den.myFactors() == vector<RingElem>(1,RingElem(QQt, "1-t")));
    CoCoA_ASSERT_ALWAYS(den.myMultiplicities() == std::vector<long>(1,4));
    CoCoA_ASSERT_ALWAYS(den.myRemainingFactor() == one(QQt));

    RingElem hp = HilbertPoly(C);
    // compare with Hilbert polynomial:   1 +4*t +8*t^2 +5*t^3
    RingElem ref_hp = RingElem(QQt, "1 + 4*t + 8*t^2 +5*t^3");
    CoCoA_ASSERT_ALWAYS(hp == ref_hp);

    // test QuasiPoly of period 1
    const QuasiPoly qp = HilbertQuasiPoly(C);
    vector<RingElem> qpv = constituents(qp);
    CoCoA_ASSERT_ALWAYS(len(qpv) == 1);
    CoCoA_ASSERT_ALWAYS(qpv[0] == ref_hp);
    CoCoA_ASSERT_ALWAYS(qp(BigInt(0)) == 1);
    CoCoA_ASSERT_ALWAYS(qp(BigInt(1)) == 18); // deg1
    CoCoA_ASSERT_ALWAYS(qp(BigInt(2)) == 81);

//----------------------------------------------------------------------
    // new very simple example not generated in deg1

    std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > m2;
    vector<vector<BigInt> > l2 = vector<vector<BigInt> >(2,vector<BigInt>(2));
    l2[0][0] = 1; l2[0][1] = 2;
    l2[1][0] = 2; l2[1][1] = 1;
    m2[libnormaliz::Type::integral_closure] = l2;
    m2[libnormaliz::Type::grading] = vector< vector<BigInt> >(1,vector<BigInt>(2,BigInt(1)));


    cone C2(m2);

    const QuasiPoly qp2 = HilbertQuasiPoly(C2);
    //cout << qp2;
    const vector<RingElem>& qpv2 = constituents(qp2);
    CoCoA_ASSERT_ALWAYS(len(qpv2) == 3);
    CoCoA_ASSERT_ALWAYS(qp2(BigInt(0)) == 1);
    CoCoA_ASSERT_ALWAYS(qp2(BigInt(1)) == 0); // deg1
    CoCoA_ASSERT_ALWAYS(qp2(BigInt(2)) == 1);
    CoCoA_ASSERT_ALWAYS(qp2(BigInt(3)) == 2);
    CoCoA_ASSERT_ALWAYS(qp2(BigInt(4)) == 1);
    CoCoA_ASSERT_ALWAYS(qp2(BigInt(5)) == 2);
    CoCoA_ASSERT_ALWAYS(qp2(BigInt(6)) == 3);

//----------------------------------------------------------------------
    // example with positive shift

    std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > m3;
    vector<vector<BigInt> > l3 = vector<vector<BigInt> >(2,vector<BigInt>(2));
    l3[0][0] = 1; l3[0][1] = 2;
    l3[1][0] = 2; l3[1][1] = 1;
    m3[libnormaliz::Type::integral_closure] = l3;
    m3[libnormaliz::Type::grading] = vector< vector<BigInt> >(1,vector<BigInt>(2,BigInt(1)));
//    m3[libnormaliz::Type::strict_inequalities] = vector< vector<BigInt> >(1,vector<BigInt>(2,BigInt(1)));
    m3[libnormaliz::Type::excluded_faces] = vector< vector<BigInt> >(1,vector<BigInt>(2,BigInt(1)));

    cone C3(m3);

    HPSeries hs3 = HilbertSeries(C3);
    RingElem ref_num3(RingElem(QQt, "t^2 + t^3 - t^4"));
    CoCoA_ASSERT_ALWAYS(num(hs3) == ref_num3);
    factorization<RingElem> den3 = DenFactors(hs3);
    vector<RingElem> ref_den_factors3(2);
    ref_den_factors3[0] = RingElem(QQt, "1-t");
    ref_den_factors3[1] = RingElem(QQt, "1-t^3");
    CoCoA_ASSERT_ALWAYS(den3.myFactors() == ref_den_factors3);
    CoCoA_ASSERT_ALWAYS(den3.myMultiplicities() == std::vector<long>(2,1));
    CoCoA_ASSERT_ALWAYS(den3.myRemainingFactor() == one(QQt));

    const QuasiPoly qp3 = HilbertQuasiPoly(C3);
    const vector<RingElem>& qpv3 = constituents(qp3);
    CoCoA_ASSERT_ALWAYS(len(qpv3) == 3);
//    CoCoA_ASSERT_ALWAYS(qp3(BigInt(0)) == 0);
    CoCoA_ASSERT_ALWAYS(qp3(BigInt(1)) == 0); // deg1
    CoCoA_ASSERT_ALWAYS(qp3(BigInt(2)) == 1);
    CoCoA_ASSERT_ALWAYS(qp3(BigInt(3)) == 2);
    CoCoA_ASSERT_ALWAYS(qp3(BigInt(4)) == 1);
    CoCoA_ASSERT_ALWAYS(qp3(BigInt(5)) == 2);
    CoCoA_ASSERT_ALWAYS(qp3(BigInt(6)) == 3);

//----------------------------------------------------------------------
    // example with negative shift, normaliz example "InhomIneq.in"

    std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > m4;
    vector<vector<BigInt> > l4 = vector<vector<BigInt> >(3,vector<BigInt>(3));
    l4[0][0] = 0; l4[0][1] = 2; l4[0][2] = 1;
    l4[1][0] = 0; l4[1][1] =-2; l4[1][2] = 3;
    l4[2][0] = 2; l4[2][1] =-2; l4[2][2] = 3;
    vector< vector<BigInt> > grading4(1,vector<BigInt>(2,BigInt(0)));
    grading4[0][0] = 1;
    m4[libnormaliz::Type::inhom_inequalities] = l4;
    m4[libnormaliz::Type::grading] = grading4;

    cone C4(m4);

// here we need to represent a negative shift!!
    HPSeries hs4 = HilbertSeries(C4);
    factorization<RingElem> den4 = DenFactors(hs4);
    RingElem ref_num4(RingElem(QQt, "1 + t"));
    CoCoA_ASSERT_ALWAYS(num(hs4) == ref_num4);
    vector<RingElem> ref_den_factors4(2);
    ref_den_factors4[0] = RingElem(QQt, "1-t");
    ref_den_factors4[1] = RingElem(QQt, "t");
    CoCoA_ASSERT_ALWAYS(den4.myFactors() == ref_den_factors4);
    CoCoA_ASSERT_ALWAYS(den4.myMultiplicities() == std::vector<long>(2,1));
// this could be another representation of the series
//    RingElem ref_num4(RingElem(QQt, "t^(-1) + 1"));
//    CoCoA_ASSERT_ALWAYS(num(hs4) == ref_num4);
//    CoCoA_ASSERT_ALWAYS(den4.myFactors() == vector<RingElem>(1,RingElem(QQt, "1-t")));
//    CoCoA_ASSERT_ALWAYS(den4.myMultiplicities() == std::vector<long>(1,1));
    CoCoA_ASSERT_ALWAYS(den4.myRemainingFactor() == one(QQt));

    const QuasiPoly qp4 = HilbertQuasiPoly(C4);
    const vector<RingElem>& qpv4 = constituents(qp4);
    CoCoA_ASSERT_ALWAYS(len(qpv4) == 1);
    CoCoA_ASSERT_ALWAYS(qp4(BigInt(0)) == 2);
    CoCoA_ASSERT_ALWAYS(qp4(BigInt(1)) == 2); // deg1
    CoCoA_ASSERT_ALWAYS(qp4(BigInt(2)) == 2);
    CoCoA_ASSERT_ALWAYS(qp4(BigInt(3)) == 2);
    CoCoA_ASSERT_ALWAYS(qp4(BigInt(4)) == 2);
    CoCoA_ASSERT_ALWAYS(qp4(BigInt(5)) == 2);
    CoCoA_ASSERT_ALWAYS(qp4(BigInt(6)) == 2);

#endif // #ifdef CoCoA_WITH_NORMALIZ
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-normaliz1.C,v 1.30 2022/02/11 09:54:00 abbott Exp $
// $Log: test-normaliz1.C,v $
// Revision 1.30  2022/02/11 09:54:00  abbott
// Summary: Updated copyright (redmine 855)
//
// Revision 1.29  2018/05/22 14:16:41  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.28  2018/05/18 12:35:19  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.27  2018/04/10 15:19:41  bigatti
// -- fixed inlcudes
//
// Revision 1.26  2017/11/08 13:44:28  abbott
// Summary: Replaced ReadExpr by RingElem(...); removed all includes of RingElemInput.H
//
// Revision 1.25  2017/04/27 14:06:35  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.24  2016/09/16 16:36:42  abbott
// Summary: Changed TEST_ASSERT into CoCoA_ASSERT_ALWAYS; removed all include assert.H lines
//
// Revision 1.23  2015/09/03 10:19:55  bigatti
// -- changes by Christof Soeger (in Aarhus)
//
// Revision 1.22  2015/05/11 13:00:53  abbott
// Summary: Put all code into namespace CoCoA (redmine #642)
// Author: JAA
//
// Revision 1.21  2014/07/14 13:17:13  abbott
// Summary: Added two consts
// Author: JAA
//
// Revision 1.20  2014/07/14 12:24:04  abbott
// Summary: Added a const
// Author: JAA
//
// Revision 1.19  2014/07/14 10:08:23  abbott
// Summary: Christof added test for quasipolys
// Author: JAA
//
// Revision 1.18  2014/06/17 10:22:07  abbott
// Summary: Removed pointless call to AsPolyRing (on RingQQt(...))
// Author: JAA
//
// Revision 1.17  2014/05/12 10:20:42  abbott
// Summary: Added some consts
// Author: JAA
//
// Revision 1.16  2014/05/09 14:56:58  bigatti
// -- new fn by Christof Soeger
//
// Revision 1.15  2012/10/08 13:53:44  bigatti
// -- more cleaning and updates by Christof Soeger
//
// Revision 1.14  2012/08/03 16:31:22  bigatti
// -- changed: procedural --> functional (by C.Soeger)
//
// Revision 1.13  2012/07/19 17:15:28  abbott
// NEEDS TO BE REWRIITEN.
// Replaced calls to NewConeLong and NewConeBigInt by calls to unified pseudo-ctor NewCone.
//
// Revision 1.12  2012/06/19 16:11:00  bigatti
// -- changed Ht1 --> Deg1 (by C.Soeger)
//
// Revision 1.11  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.10  2012/03/12 11:36:42  abbott
// Added two "const"s and improved the indentation.
//
// Revision 1.9  2011/11/07 11:29:59  bigatti
// -- changed syntax for HilbertBasis
//
// Revision 1.8  2011/10/12 15:49:21  abbott
// Simplified use of CPP macro CoCoA_WITH_NORMALIZ, but compilation may be a little slower.
//
// Revision 1.7  2011/09/30 15:57:54  bigatti
// -- moved namespace CoCoA::Normaliz inside #ifdef
//
// Revision 1.6  2011/09/30 12:55:44  bigatti
// -- introduced namespace "Normaliz" and removed Normaliz from function names
// -- input of Normaliz functions in CoCoA-5 is now a matrix instead of
//    vector<vector<BigInt>>
//
// Revision 1.5  2011/08/23 06:40:30  bigatti
// -- fixed with new name "BigInt" for old "ZZ"
//
// Revision 1.4  2011/07/20 13:49:37  bigatti
// -- added "Normaliz" postfix to Normaliz function calls
//
// Revision 1.3  2011/07/20 12:45:12  bigatti
// -- new normaliz interface (not yet public)
//
// Revision 1.2  2011/07/20 10:12:17  bigatti
// -- fixed compilation without normaliz
//
// Revision 1.1  2011/07/19 16:24:16  bigatti
// -- first import
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

