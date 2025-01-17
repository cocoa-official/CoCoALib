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
