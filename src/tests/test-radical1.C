// Copyright (c) 2006,2023  John Abbott,  Anna M. Bigatti
// original authors: Nicolas Jagersma (testData6Nov.C), Alice Moallemy

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
#include "CoCoA/QuotientRing.H"  // for NewZZmod
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-ideal.H"  // for radical
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/time.H"
#include "CoCoA/verbose.H"  // for SetVerbosityLevel

//#include <getopt.h>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

//----------------------------------------------------------------------
// This test the function radical 0-dim or not
//----------------------------------------------------------------------

namespace CoCoA
{
  
  void TestRadical (ideal& RadI, ideal in)
  {
    VerboseLog VERBOSE("TestRadical");
    VERBOSE(20) << "IsZeroDim?  " << IsZeroDim(in) << endl;
    double t0;
    if (IsZeroDim(in))
    {
      t0 = CpuTime();
      /*[[maybe_unused]]*/ const bool IsRad = IsRadical_tmp(in); // attribute needs C++17
      VERBOSE(20) << "Time for IsRadical: " << CpuTime()-t0 << "s " << endl;
    }
    t0 = CpuTime();
    RadI = radical(in);
    VERBOSE(20) << "Time for radical:  " << CpuTime()-t0 << "s " << endl;
    t0 = CpuTime();
    ideal radOfRad = radical(ideal(gens(RadI)));  // force the computation
    VERBOSE(20) << "Time for radOfRad: " << CpuTime()-t0 << "s " << endl;
    VERBOSE(20) << "Checking rad = rad(rad)? ......... ";
    if (RadI != radOfRad)  CoCoA_THROW_ERROR2(ERR::ShouldNeverGetHere, "radI != radOfRad");
    VERBOSE(20) << " done!" << endl;
  }

  
  void FillExampleList(std::vector<ideal>& examp)
  {
    ring QQ = RingQQ();

     SparsePolyRing P1 = NewPolyRing(QQ, symbols("x"));
     SparsePolyRing P2 = NewPolyRing(QQ, symbols("x, y"));
     SparsePolyRing P3 = NewPolyRing(QQ, symbols("x, z, y"));
     SparsePolyRing P4 = NewPolyRing(QQ, symbols("x, y, z, t"));


     ring ZZMod101 = NewZZmod(101);
     SparsePolyRing ZZ101_4 = NewPolyRing(ZZMod101, symbols("x, y, z, t"));

     //Example 6.1 in paper
     RingElem F1 = ReadExpr(ZZ101_4, "x*y*z*t +83*y^3 +73*x^2 -85*z^2 -437*t");
     RingElem F2 = ReadExpr(ZZ101_4, "y^3*z*t +z -t");
     RingElem F3 = ReadExpr(ZZ101_4, "t^4 +z*t^2 -324*z^3 +94*x^2 +76*y");
     RingElem F4 = ReadExpr(ZZ101_4, "x^11*z +26*t^3 +625*y");
      
     ideal I = ideal(F1, F2, F3, F4);
     examp.push_back(I); // case 0

     
     SparsePolyRing ZZ101_6 = NewPolyRing(ZZMod101, symbols("a_1, a_2, a_3, a_4, a_5, a_6"));
     //Example 6.2 in paper
     std::vector<RingElem> SymPoly_mod;
     SymPoly_mod.push_back(ReadExpr(ZZ101_6, "a_1 +a_2 +a_3 +a_4 +a_5 +a_6"));
     SymPoly_mod.push_back(ReadExpr(ZZ101_6, "a_1*a_2 +a_1*a_3 +a_2*a_3 +a_1*a_4 +a_2*a_4 +a_3*a_4 +a_1*a_5 +a_2*a_5 +a_3*a_5 +a_4*a_5 +a_1*a_6 +a_2*a_6 +a_3*a_6 +a_4*a_6 +a_5*a_6"));
     SymPoly_mod.push_back(ReadExpr(ZZ101_6, "a_1*a_2*a_3 +a_1*a_2*a_4 +a_1*a_3*a_4 +a_2*a_3*a_4 +a_1*a_2*a_5 +a_1*a_3*a_5 +a_2*a_3*a_5 +a_1*a_4*a_5 +a_2*a_4*a_5 +a_3*a_4*a_5 +a_1*a_2*a_6 +a_1*a_3*a_6 +a_2*a_3*a_6 +a_1*a_4*a_6 +a_2*a_4*a_6 +a_3*a_4*a_6 +a_1*a_5*a_6 +a_2*a_5*a_6 +a_3*a_5*a_6 +a_4*a_5*a_6"));
     SymPoly_mod.push_back(ReadExpr(ZZ101_6, "a_1*a_2*a_3*a_4 +a_1*a_2*a_3*a_5 +a_1*a_2*a_4*a_5 +a_1*a_3*a_4*a_5 +a_2*a_3*a_4*a_5 +a_1*a_2*a_3*a_6 +a_1*a_2*a_4*a_6 +a_1*a_3*a_4*a_6 +a_2*a_3*a_4*a_6 +a_1*a_2*a_5*a_6 +a_1*a_3*a_5*a_6 +a_2*a_3*a_5*a_6 +a_1*a_4*a_5*a_6 +a_2*a_4*a_5*a_6 +a_3*a_4*a_5*a_6"));
     SymPoly_mod.push_back(ReadExpr(ZZ101_6, "a_1*a_2*a_3*a_4*a_5 +a_1*a_2*a_3*a_4*a_6 +a_1*a_2*a_3*a_5*a_6 +a_1*a_2*a_4*a_5*a_6 +a_1*a_3*a_4*a_5*a_6 +a_2*a_3*a_4*a_5*a_6-7"));
     SymPoly_mod.push_back(ReadExpr(ZZ101_6, "a_1*a_2*a_3*a_4*a_5*a_6-1"));

     ideal I_6 = ideal(SymPoly_mod);
     //ideal r_6 = radical(I_6);
     //cout << "Radical = " << r_6 << endl;
     examp.push_back(I_6); // case 1

     
     ring ZZMod32003 = NewZZmod(32003);
     SparsePolyRing ZZ32003 = NewPolyRing(ZZMod32003, symbols("x, y, z, t"));

     // Example 6.3 paper
     ideal I_32003 = ideal(ReadExpr(ZZ32003, "x*y*z*t^5 +x^3 +73*y^2 -z^2 -2*t"), ReadExpr(ZZ32003, "x^6*z*t+z-t"),
			    ReadExpr(ZZ32003, "z*t^2 +7*x*y^2 -34*z^3 -2*t^4"), ReadExpr(ZZ32003, "y^4*z +x +26*t^3"));
     //ideal rad_I_32003 = radical(I_32003);
     examp.push_back(I_32003); // case 2
     //cout << "Radical = " << rad_I_32003 << endl;

     
     SparsePolyRing ZZ101_3 = NewPolyRing(ZZMod101, symbols("x, y, z"));

     // Example 6.4 paper
     ideal J1 = ideal(ReadExpr(ZZ101_3, "(x^7-y-3*z)^2"), ReadExpr(ZZ101_3, "x*y^5 -7*z^2 -2"), ReadExpr(ZZ101_3, "y*z^6-x-z+14"));
     ideal J2_t = ideal(ReadExpr(ZZ101_3, "x"), ReadExpr(ZZ101_3, "y"), ReadExpr(ZZ101_3, "z"));
     ideal J2 = J2_t * J2_t;

     ideal J = intersect(J1, J2);

     //examp.push_back(ideal(one(RingOf(J2))));
     examp.push_back(J); // case 3
     //ideal radJ = radical(J);
     //cout << "Radical = " << radJ << endl;


     ring ZZ23 = NewZZmod(23);
     SparsePolyRing ZZ23P = NewPolyRing(ZZ23, symbols("x, y, z"));
     // Example 6.5 paper
     ideal I_ZZ23 = ideal(ReadExpr(ZZ23P, "x^16 +8*x^15 -6*x^14 -8*x^13 +4*x^12 -4*x^11 +5*x^10 +8*x^9 +5*x^8 -4*x^7 +5*x^6 +2*x^5 -7*x^4 +4*x^3 +10*x^2 +3*x +8"), ReadExpr(ZZ23P,"y^5 -7*y^4 +2*y^3 +11*y^2 -y +5"), ReadExpr(ZZ23P, "z^11 +9*z^10 -9*z^9 +7*z^8 -8*z^7 -4*z^6 +9*z^5 +z^4 -5*z^3 +7*z^2 +z +10"));

     //ideal radI_ZZ23 = radical(I_ZZ23);
     examp.push_back(I_ZZ23); // case 4
     //cout << "rad = " << radI_ZZ23 << endl;


     // Example in Q
     // Example 6.7 paper doesn't work?
     ideal IQ = ideal(ReadExpr(P4, "x*y*z*t+83*x^3+73*y^2-85*z^2-437*t"), ReadExpr(P4, "x^3*z*t+z-t"), ReadExpr(P4, "z*t^2+76*x+94*y^2-324*z^3-255*t^4"), ReadExpr(P4, "y^2*z+625*x+26*t^3"));
     //ideal rad_IQ = radical(IQ);
     examp.push_back(IQ); // case 5
     //cout << "rad = " << rad_IQ << endl;

     
     // Example 6.8 paper doesn't work?
     ideal IQQ = ideal(ReadExpr(P4, "x^4 +83*x^3 +73*y^2 -85*z^2 -437*t"), ReadExpr(P4, "y^3-x"), ReadExpr(P4, "z^3 +z -t"), ReadExpr(P4, "t^3-324*z^2 +94*y^2 +76*x"));
     examp.push_back(IQQ); // case 6
     //ideal radQQ = radical(IQQ);
     //cout << "rad = " << radQQ << endl;

     
     //Bsp 6.9 paper doesn't work??
     ideal IQQQ = ideal(ReadExpr(P4, "x^4 +83*z^3 +73*y^2 -t^2 -437*t"), ReadExpr(P4, "y^3-z-t"), ReadExpr(P4, "z^3 +x -t"), ReadExpr(P4, "t^4-12*z^2 +77*y^2 +15*x"));
     examp.push_back(IQQQ); // case 7
     //ideal radQQQ = Radical0dim(IQQQ);
     //cout << "rad = " << radQQQ << endl;

     
     //Example 6.10 paper
     SparsePolyRing PP11 = NewPolyRing(QQ, symbols("a_1, a_2, a_3, a_4, a_5"));

     std::vector<RingElem> SymPolyMod2;
     SymPolyMod2.push_back(ReadExpr(PP11, "a_1 +a_2 +a_3 +a_4 +a_5"));
     SymPolyMod2.push_back(ReadExpr(PP11, "a_1*a_2 +a_1*a_3 +a_2*a_3 +a_1*a_4 +a_2*a_4 +a_3*a_4 +a_1*a_5 +a_2*a_5 +a_3*a_5 +a_4*a_5"));
     SymPolyMod2.push_back(ReadExpr(PP11, "a_1*a_2*a_3 +a_1*a_2*a_4 +a_1*a_3*a_4 +a_2*a_3*a_4 +a_1*a_2*a_5 +a_1*a_3*a_5 +a_2*a_3*a_5 +a_1*a_4*a_5 +a_2*a_4*a_5 +a_3*a_4*a_5"));
     SymPolyMod2.push_back(ReadExpr(PP11, "a_1*a_2*a_3*a_4 +a_1*a_2*a_3*a_5 +a_1*a_2*a_4*a_5 +a_1*a_3*a_4*a_5 +a_2*a_3*a_4*a_5 +1"));
     SymPolyMod2.push_back(ReadExpr(PP11, "a_1*a_2*a_3*a_4*a_5 -2"));

     ideal symI2 = ideal(SymPolyMod2);
     examp.push_back(symI2); // case 8
     //     ideal radSymI2 = radical(symI2);
     //cout << "rad = " << radSymI2 << endl;

     
     // Example 6.11 paper
     SparsePolyRing symRing = NewPolyRing(QQ, symbols("a_1, a_2, a_3, a_4, a_5, a_6"));
     std::vector<RingElem> SymPoly_mod3;
     SymPoly_mod3.push_back(ReadExpr(symRing, "a_1 +a_2 +a_3 +a_4 +a_5 +a_6"));
     SymPoly_mod3.push_back(ReadExpr(symRing, "a_1*a_2 +a_1*a_3 +a_2*a_3 +a_1*a_4 +a_2*a_4 +a_3*a_4 +a_1*a_5 +a_2*a_5 +a_3*a_5 +a_4*a_5 +a_1*a_6 +a_2*a_6 +a_3*a_6 +a_4*a_6 +a_5*a_6"));
     SymPoly_mod3.push_back(ReadExpr(symRing, "a_1*a_2*a_3 +a_1*a_2*a_4 +a_1*a_3*a_4 +a_2*a_3*a_4 +a_1*a_2*a_5 +a_1*a_3*a_5 +a_2*a_3*a_5 +a_1*a_4*a_5 +a_2*a_4*a_5 +a_3*a_4*a_5 +a_1*a_2*a_6 +a_1*a_3*a_6 +a_2*a_3*a_6 +a_1*a_4*a_6 +a_2*a_4*a_6 +a_3*a_4*a_6 +a_1*a_5*a_6 +a_2*a_5*a_6 +a_3*a_5*a_6 +a_4*a_5*a_6"));
     SymPoly_mod3.push_back(ReadExpr(symRing, "a_1*a_2*a_3*a_4 +a_1*a_2*a_3*a_5 +a_1*a_2*a_4*a_5 +a_1*a_3*a_4*a_5 +a_2*a_3*a_4*a_5 +a_1*a_2*a_3*a_6 +a_1*a_2*a_4*a_6 +a_1*a_3*a_4*a_6 +a_2*a_3*a_4*a_6 +a_1*a_2*a_5*a_6 +a_1*a_3*a_5*a_6 +a_2*a_3*a_5*a_6 +a_1*a_4*a_5*a_6 +a_2*a_4*a_5*a_6 +a_3*a_4*a_5*a_6"));
     SymPoly_mod3.push_back(ReadExpr(symRing, "a_1*a_2*a_3*a_4*a_5 +a_1*a_2*a_3*a_4*a_6 +a_1*a_2*a_3*a_5*a_6 +a_1*a_2*a_4*a_5*a_6 +a_1*a_3*a_4*a_5*a_6 +a_2*a_3*a_4*a_5*a_6-7"));
     SymPoly_mod3.push_back(ReadExpr(symRing, "a_1*a_2*a_3*a_4*a_5*a_6-1"));

     ideal I_62 = ideal(SymPoly_mod3);
     //ideal r_62 = radical(I_62);
     //cout << "rad = " << r_62 << endl;
     examp.push_back(I_62); // case 9

     
     // Example 6.13 paper  
     ideal II13 = ideal(ReadExpr(P3, "x^7-y-3*z"), ReadExpr(P3, "x*y^5-5057*z^2-2"), ReadExpr(P3, "y*z^6-x-z+14"));
     examp.push_back(II13); // case 10
     //cout << "test" << endl;   
     //ideal radII13 = radical(II13);
     //cout << "rad = " << radII13 << endl;

     
     /*
     // radical of polynomial
     RingElem f11 = ReadExpr(P4, "x^2*t +z^2*t +1");
     RingElem f12 = ReadExpr(P4, "y^5-y");
     RingElem f21 = ReadExpr(P3,"y^2 +x*z +y*z +x^2 -y^2 +x*z +z^2 +x +y");
     RingElem f31 = ReadExpr(P4,"x^2*t +z^2*t +1 +x^3 +t^3 +y");
     RingElem f41 = ReadExpr(P2, "x^4");
     RingElem f51 = ReadExpr(P2, "x^2 + 2*x*y + y^2");
     
     CoCoA_ASSERT_ALWAYS(radical(f11) == ReadExpr(P4, "x^2*t +z^2*t +1"));
     CoCoA_ASSERT_ALWAYS(radical(f12) == ReadExpr(P4, "y^5 -y"));
     CoCoA_ASSERT_ALWAYS(radical(f21) == ReadExpr(P3, "x^2 +2*x*z +y*z +z^2 +x +y"));
     CoCoA_ASSERT_ALWAYS(radical(f31) == ReadExpr(P4, "x^3 +x^2*t +z^2*t +t^3 +y +1"));
     CoCoA_ASSERT_ALWAYS(radical(f41) == ReadExpr(P2, "x")); 
     CoCoA_ASSERT_ALWAYS(radical(f51) == ReadExpr(P2, "x +y"));

     //ideals in QQ[x,y]:
     ideal J11 = ideal(ReadExpr(P2, "y^2 + 2"), ReadExpr(P2, "x^3 + x*y^2"));
     ideal J22 = ideal(ReadExpr(P2, "x^3 + x^2 + 1"), ReadExpr(P2, "y^3 + y"));
     ideal J3 = ideal(ReadExpr(P2, "x^3 + y^3 + 1"), ReadExpr(P2, "y^3 + y^2"));
     ideal J4 = ideal(ReadExpr(P2, "x^3 + x + 1"), ReadExpr(P2, "y^3 + y^2"));
     ideal J5 = ideal(ReadExpr(P2, "y^3 + x^2 + 1"), ReadExpr(P2, "x*y^2 + y^2"));
     ideal J6 = ideal(ReadExpr(P2, "x^2 + y^2 + 1"), ReadExpr(P2, "x*y^2 + y^2"));
     */
     
     ideal N1 = ideal(ReadExpr(P2, "x^4"));
     ideal N2 = ideal(ReadExpr(P2, "x^2 -x"));

     examp.push_back(N1); // case 11
     examp.push_back(N2); // case 12
     
     /*
     CoCoA_ASSERT_ALWAYS(radical(J11) == ideal(ReadExpr(P2, "y^2 +2"), ReadExpr(P2, "x^3 +x*y^2")));
     CoCoA_ASSERT_ALWAYS(radical(J22) == ideal(ReadExpr(P2, "x^3 +x^2 +1"), ReadExpr(P2, "y^3 +y")));
     CoCoA_ASSERT_ALWAYS(radical(J3) == ideal(ReadExpr(P2, "y^3 +y^2 +x*y^2"), ReadExpr(P2, "x^3 -y^2 +1"), ReadExpr(P2, "y^2 +y")));
     CoCoA_ASSERT_ALWAYS(radical(J4) == ideal(ReadExpr(P2, "y^3 +y^2"), ReadExpr(P2, "x^3 +x +1"), ReadExpr(P2, "y^2 +y")));
     CoCoA_ASSERT_ALWAYS(radical(J5) == ideal(ReadExpr(P2, "y^3 +x^2 +1"), ReadExpr(P2, "x*y^2 +y^2"), ReadExpr(P2, "x^3 +x^2 +x +1"), ReadExpr(P2, "-x^2*y +y")));
     CoCoA_ASSERT_ALWAYS(radical(J6) == ideal(ReadExpr(P2, "x^2 +y^2 +1"), ReadExpr(P2, "x*y^2 +y^2"), ReadExpr(P2, "y^4 +2*y^2"), ReadExpr(P2, "y^3 +2*y")));
     
     CoCoA_ASSERT_ALWAYS(radical(N1) == ideal(ReadExpr(P2, "x")));
     CoCoA_ASSERT_ALWAYS(radical(N2) == N2);
         

     //ideals in QQ[x, y, z]:
     ideal A = ideal(ReadExpr(P3,"x"), ReadExpr(P3, "z^2"), ReadExpr(P3, "y"));
     ideal C = ideal(ReadExpr(P3,"y-1"), ReadExpr(P3, "z^2-z"), ReadExpr(P3, "x*z-x"), ReadExpr(P3, "x^2-x"));
     ideal CC = ideal(ReadExpr(P3, "x^3 +y^4 +x*y*z"), ReadExpr(P3, "y^5 +y^2*z^2"), ReadExpr(P3, "z^5"));
     */

     
     ideal B = ideal(ReadExpr(P3,"3*y^4 + 2*z*x"), ReadExpr(P3, "x*y^3*z^2"), ReadExpr(P3, "z+y"));
     examp.push_back(B); // case 13
     //cout << " Radical of B = " << radical(B) << endl;

     
     /*
     CoCoA_ASSERT_ALWAYS(radical(A) == ideal(ReadExpr(P3, "y"), ReadExpr(P3, "x"), ReadExpr(P3, "z^2"), ReadExpr(P3, "z")));
     CoCoA_ASSERT_ALWAYS(radical(C) == ideal(ReadExpr(P3, "y -1"), ReadExpr(P3, "z^2 -z"), ReadExpr(P3, "x*z -x"), ReadExpr(P3, "x^2 -x")));
     CoCoA_ASSERT_ALWAYS(radical(CC) == ideal(ReadExpr(P3, "y"), ReadExpr(P3, "x"), ReadExpr(P3, "z^5"), ReadExpr(P3, "z")));
     CoCoA_ASSERT_ALWAYS(radical(B) == ideal(ReadExpr(P3, "z"), ReadExpr(P3, "y")));
     */

     
     //ideals in QQ[x, y, z, t]:
     ideal K1 = ideal(ReadExpr(P4, "t^3 + y^2 +1"), ReadExpr(P4, "y*z^2 + z*y"), ReadExpr(P4, "y*t^2 + y*t"), ReadExpr(P4, "x^3 +1"));
     ideal K2 = ideal(ReadExpr(P4,"x^2*t +z^2*t +1"), ReadExpr(P4,"x^3 +t^3"), ReadExpr(P4,"x*y^2+x"), ReadExpr(P4,"y*z^2 +y*z"));
     ideal K3 = ideal(ReadExpr(P4, "y*z^2 + z*y"), ReadExpr(P4, "y*t^2 + y*t"));

     examp.push_back(K1); // case 14
     examp.push_back(K2); // case 15
     examp.push_back(K3); // case 16

     
     /*
     std::vector<RingElem> resK1;
     resK1.push_back(ReadExpr(P4, "y*t"));
     resK1.push_back(ReadExpr(P4, "y^3 +y"));
     resK1.push_back(ReadExpr(P4, "t^3 +y^2 +1"));
     resK1.push_back(ReadExpr(P4, "y*z^2 +y*z"));
     resK1.push_back(ReadExpr(P4, "x^3 +1"));
     
     CoCoA_ASSERT_ALWAYS(radical(K1) == ideal(resK1));
     CoCoA_ASSERT_ALWAYS(radical(K2) == ideal(ReadExpr(P4, "x^2*t +z^2*t +1"), ReadExpr(P4, "x^3 +t^3"), ReadExpr(P4, "x*y^2 +x"), ReadExpr(P4, "y*z^2 +y*z")));
     CoCoA_ASSERT_ALWAYS(radical(K3) == ideal(ReadExpr(P4, "y*t^2 +y*t"), ReadExpr(P4, "y*z^2 +y*z")));  
     */

     
     // Butcher
     SparsePolyRing P4b = NewPolyRing(QQ, symbols("x,y,z,t,u,v,w"));         

     RingElem f_1 = ReadExpr(P4b, "2*t*w + 2*u*z + 2*v*y - 2*w^2 -w - 1");
     RingElem f_2 = ReadExpr(P4b, "- 3*t*w^2 - t + 3*u^2*z + 3*v^2*y + 3*w^3 + 3*w^2 + 4*w");
     RingElem f_3 = ReadExpr(P4b, "- 6*t*w^2 - 3*t*w -t + 6*v*x*z + 6*w^3 + 6*w^2 + 4*w");
     RingElem f_4 = ReadExpr(P4b, "4*t*w^3 + 4*t*w + 4*u^3*z + 4*v^3*y - 4*w^4 - 6*w^3 - 10*w^2 - w - 1");
     RingElem f_5 = ReadExpr(P4b, "8*t*w^3 + 4*t*w^2 + 4*t*w + 8*u*v*x*z - 8*w^4 - 12*w^3 - 14*w^2 - 3*w - 1");
     RingElem f_6 = ReadExpr(P4b, "12*t*w^3 + 12*t*w^2 + 8*t*w + 12*v^2*x*z - 12*w^4 - 18*w^3 - 14*w^2 - w - 1");
     RingElem f_7 = ReadExpr(P4b, "- 24*t*w^3 - 24*t*w^2 - 8*t*w + 24*w^4 + 36*w^3 + 26*w^2 + 7*w + 1"); 
     std::vector<RingElem> fis;
     fis.push_back(f_1);
     fis.push_back(f_2);
     fis.push_back(f_3);       
     fis.push_back(f_4); 
     fis.push_back(f_5);  
     fis.push_back(f_6); 
     fis.push_back(f_7); 
     ideal I_1 = ideal(fis);
     examp.push_back(I_1); // case 17
     /*
     //timing cmd
       ideal rI_1 = radical(I_1);
     cout << "here the radical of butcher = " << rI_1 << endl;
     */

     
     // hairer1
     SparsePolyRing P5h = NewPolyRing(QQ, symbols("a,b,c,d,e,f,g,h"));

     RingElem g_1 = ReadExpr(P5h, "a-f");
     RingElem g_2 = ReadExpr(P5h, "b-g-h");
     RingElem g_3 = ReadExpr(P5h, "c + d + e-1");
     RingElem g_4 = ReadExpr(P5h, "2*a*d + 2*b*c-1");
     RingElem g_5 = ReadExpr(P5h, "3*a^2*d + 3*b^2*c-1");
     RingElem g_6 = ReadExpr(P5h, "6*a*c*g-1");

     std::vector<RingElem> gis;
     gis.push_back(g_1);
     gis.push_back(g_2);
     gis.push_back(g_3);
     gis.push_back(g_4);
     gis.push_back(g_5);
     gis.push_back(g_6);
     ideal Ih_2 = ideal(gis);
     examp.push_back(Ih_2); // case 18

     
     /*
     RingElem gg1 = ReadExpr(P5h, "f^2*g -f*g +(1/3)*g^2 +(-1/3)*f*h +(2/3)*g*h +(1/3)*h^2");
     RingElem gg2 = ReadExpr(P5h, "e*g^2 +2*e*g*h +e*h^2 +(4/3)*d*g +2*e*g +(-3/2)*f*g -g^2 +(2/3)*e*h -2*g*h -h^2 +(1/3)*h +1/3");
     RingElem gg3 = ReadExpr(P5h, "e*f*h +(-2/3)*d*g -e*g +(-1/3)*e*h -f*h +(1/2)*f +g +(5/6)*h -2/3");
     RingElem gg4 = ReadExpr(P5h, "e*f*g +(2/3)*d*g +e*g -f*g +(1/3)*e*h +(-1/2)*g +(-1/3)*h +1/3");
     RingElem gg5 = ReadExpr(P5h, "d*g^2 +d*g*h -e*g*h -e*h^2 -2*d*g -3*e*g +(3/2)*f*g -e*h +g*h +h^2 +g -1/2");
     RingElem gg6 = ReadExpr(P5h, "d*f -d*g -e*g -d*h -e*h +g +h -1/2");
     RingElem gg7 = ReadExpr(P5h, "d*e*g*h +e^2*g*h +d*e*h^2 +e^2*h^2 +(2/3)*d*e*g +e^2*g +(1/3)*e^2*h -d*g*h -2*e*g*h -d*h^2 -2*e*h^2 +(1/3)*d*g -e*g +d*h +(1/3)*e*h +g*h +h^2 +(-1/3)*d +(1/6)*e +(-2/3)*h +1/12");
     RingElem gg8 = ReadExpr(P5h, "d^2*g +(5/2)*d*e*g +(3/2)*e^2*g +(1/2)*d*e*h +(1/2)*e^2*h +(-7/4)*d*g +(-9/4)*e*g +(-1/2)*d*h -e*h +(1/2)*d +(1/4)*e +(3/4)*g +(1/2)*h -1/4");
     RingElem gg9 = ReadExpr(P5h, "d^2*e*h^2 +2*d*e^2*h^2 +e^3*h^2 +(1/3)*d*e^2*g +(1/2)*e^3*g +(1/6)*e^3*h -d^2*h^2 +(-13/4)*d*e*h^2 +(-9/4)*e^2*h^2 +(-2/3)*d*e*g +(-5/4)*e^2*g +d^2*h +(3/2)*d*e*h +(1/2)*e^2*h +(5/4)*d*h^2 +(3/2)*e*h^2 +(-1/3)*d^2 +(-2/3)*d*e +(1/12)*e^2 +(1/3)*d*g +e*g +(-3/4)*d*h +(-3/4)*e*h +(-1/4)*h^2 +(1/6)*d +(1/12)*e +(-1/4)*g +(1/12)*h +1/48");
     RingElem gg10 = ReadExpr(P5h, "c +d +e -1");
     RingElem gg11 = ReadExpr(P5h, " b -g -h");
     RingElem gg12 = ReadExpr(P5h, "a -f");

     std::vector<RingElem> resOfIdrI_2;
     resOfIdrI_2.push_back(gg1);
     resOfIdrI_2.push_back(gg2);
     resOfIdrI_2.push_back(gg3);
     resOfIdrI_2.push_back(gg4);
     resOfIdrI_2.push_back(gg5);
     resOfIdrI_2.push_back(gg6);
     resOfIdrI_2.push_back(gg7);
     resOfIdrI_2.push_back(gg8);
     resOfIdrI_2.push_back(gg9);
     resOfIdrI_2.push_back(gg10);
     resOfIdrI_2.push_back(gg11);
     resOfIdrI_2.push_back(gg12);

     ideal resOfIDRI_2 = ideal(resOfIdrI_2);

     //timing cmd
     ideal rI_2 = radical(Ih_2);
     CoCoA_ASSERT_ALWAYS(radical(Ih_2) == resOfIDRI_2);
     //cout << "Ih_2 = resOfIDRI_2? " << (radical(Ih_2) == resOfIDRI_2) << endl;
     */
     

     // hcyclic5
     // sorting incorrect :(

     SparsePolyRing P6h = NewPolyRing(QQ, symbols("x_1,x_2,x_3,x_4,x_5,w"));
     RingElem h_1 = ReadExpr(P6h, "x_1 + x_2 + x_3 + x_4 + x_5");
     RingElem h_2 = ReadExpr(P6h, "x_1*x_2 + x_1*x_5 + x_2*x_3 + x_3*x_4 + x_4*x_5");
     RingElem h_3 = ReadExpr(P6h, "x_1*x_2*x_3 + x_1*x_2*x_5 + x_1*x_4*x_5 + x_2*x_3*x_4 + x_3*x_4*x_5");
     RingElem h_4 = ReadExpr(P6h, "x_1*x_2*x_3*x_4 + x_1*x_2*x_3*x_5 + x_1*x_2*x_4*x_5 + x_1*x_3*x_4*x_5 + x_2*x_3*x_4*x_5");
     RingElem h_5 = ReadExpr(P6h, "-w^5 + x_1*x_2*x_3*x_4*x_5");

     std::vector<RingElem> his;
     his.push_back(h_1);
     his.push_back(h_2);
     his.push_back(h_3);
     his.push_back(h_4);
     his.push_back(h_5);
     
     ideal Ih_3 = ideal(his);  
     examp.push_back(Ih_3); // case 19
     /*
     //timing cmd
     ideal rI_3 = radical(Ih_3);
     cout << "hcyclics " << rI_3 << endl;
     */

     
     // hemmecke
     SparsePolyRing Ph7 = NewPolyRing(QQ, symbols("x,y,z,w"));
     RingElem c_1 = ReadExpr(Ph7, "x^2-z^20-z^10");
     RingElem c_2 = ReadExpr(Ph7, "x*y^3-z^30-z^10");
     RingElem c_3 = ReadExpr(Ph7, "-w^40*x^4 + y^6");
     ideal Ih_4 = ideal(c_1, c_2, c_3);

     examp.push_back(Ih_4); // case 20
     
     /*
     //timing cmd
     ideal rI_4 = radical(Ih_4);
     cout << " hier hemmecke = " << rI_4 << endl;
     */

     
     // this thing has an entirely different problem...
     // quadfor2
     SparsePolyRing P8q = NewPolyRing(QQ, symbols("w_1,w_2,w_3,x_1,x_2"));
     RingElem d_1 = ReadExpr(P8q, "w_1 + w_2-1");
     RingElem d_2 = ReadExpr(P8q, "w_1*x_1 + w_2*x_2");
     RingElem d_3 = ReadExpr(P8q, "3*w_1*x_1^2 + 3*w_2*x_2^2-2");
     RingElem d_4 = ReadExpr(P8q, "w_1*x_1^3 + w_2*x_2^3");
     ideal I_5q = ideal(d_1, d_2, d_3, d_4);
     examp.push_back(I_5q); // case 21
     
     /*
     //timing cmd
     ideal rI_5 = radical(I_5q);
     cout << "quadfor2 = " << rI_5 << endl;
     */

     
     // hietarinta1
     SparsePolyRing P9h = NewPolyRing(QQ, symbols("a,c,j,l,m,n,p,v,g,h"));

     RingElem q_0 = ReadExpr(P9h, "-a*c + a + c^2*p + c*l-h-p");
     RingElem q_1 = ReadExpr(P9h, "a*c*h + c^2 + c*n + h*m");
     RingElem q_2 = ReadExpr(P9h, "-a^2*c + a^2 + a*c*g + a*c*l-a*h-a*m + 2*c^2-c*h*l + h*l-2");
     RingElem q_3 = ReadExpr(P9h, "-a*c^2 + a*c*j + a*c*n-c^2*m + c^2*v-c*h*n-c*m + h*n");
     RingElem q_4 = ReadExpr(P9h, "-(a*l + c*g*l + j + 1)");
     RingElem q_5 = ReadExpr(P9h, "-c*g*n + c*g + c*h-c*l-c*m + j*m");
     RingElem q_6 = ReadExpr(P9h, "-a*j-a*n + a + c*g-c*j*l + g + j*l-v");
     RingElem q_7 = ReadExpr(P9h, "-c*j*n + c*j-c*n + j*n");
     RingElem q_8 = ReadExpr(P9h, "a*l + a*p + c*m*p + c + l*m-l*n*p-l*p-n-2");
     RingElem q_9 = ReadExpr(P9h, "-c*h + c*l + c*m + c*p + 2*m*n + m-n^2*p");
     RingElem q_10 = ReadExpr(P9h, "a*c-a-c*g + 2*c*m-h*l*m + h-l*n + m*n + p");
     RingElem q_11 = ReadExpr(P9h, "a*m*n + c^2-c*j-c*l*m-c*m^2 + c*m*v-h*m*n-m^2*n-n^2 + n");
     RingElem q_12 = ReadExpr(P9h, "a^2*l-a*l^2-a*n + 2*a-c*l + c*n*p-g*l*m + g-l + m-v");
     RingElem q_13 = ReadExpr(P9h, "a*c*l + a*m*n + c*h*l-c*l^2 + c*n + 2*c-g*m*n-l*m*n-m^2 + m*v-n^2");
     RingElem q_14 = ReadExpr(P9h, "a*l + c*g*l + c*n-j*l*m + j-l*n*v + n^2 + 1");
     RingElem q_15 = ReadExpr(P9h, "a*n^2 + c*j*l-c*l*n + c*n*v-j*m*n-m*n^2-n^2*v + n*v");
     RingElem q_16 = ReadExpr(P9h, "a*g + a*h-a*l-a*m + c*m*p + c-g*p + h*n*p-j*l*p-j + n-1");
     RingElem q_17 = ReadExpr(P9h, "c*g-c*l + h*l*m + h*n + j*m-j*n*p");
     RingElem q_18 = ReadExpr(P9h, "a*c-a + c*h-c*l + c*m-g*n-g-h^2*l + h*l^2 + h*n + h-j*l + j*m-l");
     RingElem q_19 = ReadExpr(P9h, "a*j*m + c^2-c*g*m-c*h*m + c*m*v-c*n-h^2*n + h*l*n + h*n*v-j*m^2-2*j*n");
     RingElem q_20 = ReadExpr(P9h, "-g*h*l-g*n + h-j*l + j*n*p-m");
     RingElem q_21 = ReadExpr(P9h, "-a*g*n + g^2*n + g*l + g*m-g*v-h*j*l + j*l^2-j*l*v + 2*j*n + j-1");
     RingElem q_22 = ReadExpr(P9h, "-c*g*n + g*j*n-h*j*n + j*l*n-j*m*n + j*m");
     RingElem q_23 = ReadExpr(P9h, "-a^2*p + a*l*p + 3*a-l*p*v + m^2*p + 2*m-v");
     RingElem q_24 = ReadExpr(P9h, "-a*c*p-a*m + c*l*p + 2*c-h*m + h*n*p + l*m + m*v-n*p*v + 2*n");
     RingElem q_25 = ReadExpr(P9h, "-a*c*p-a*m + c-g*m + g*n*p-h*l-h*p + l^2-l*m + l*p-l*v + m^2 + m*v + n");
     RingElem q_26 = ReadExpr(P9h, "-c^2*p-2*c*m-h*n-j*m + j*n*p + l*n-m*n");
     RingElem q_27 = ReadExpr(P9h, "-a*l-a*p-g*l-g*p-l*m-l*v + m*n*p + n*p*v");
     RingElem q_28 = ReadExpr(P9h, "-c*l-c*p-g*n + h*l*m + h*n-m");
     RingElem q_29 = ReadExpr(P9h, "-a*n + g*l*m + g*n-j*l-j*p + l^2*v-l*n-l*v^2 + m*n + 2*n*v-v");
     RingElem q_30 = ReadExpr(P9h, "-c*n + j*l*m + l*n*v-n^2");
     RingElem q_31 = ReadExpr(P9h, "-a*c*p-a*h + c*h*p + 2*c-g*h-g*p + h*j*p + j-1");
     RingElem q_32 = ReadExpr(P9h, "-a^2*h-2*a*c + a*h^2 + a*j-c*g + c*h-c*m + c*v + g*h*m + h-j*p");
     RingElem q_33 = ReadExpr(P9h, "a*c-a-c^2*p-c*l + g*h*l-g*j-2*g + h*j + m");
     RingElem q_34 = ReadExpr(P9h, "-a*c*h-c^2-c*n + g*h*n + h*j*v-h*m-j^2-j");
     RingElem q_35 = ReadExpr(P9h, "-c*g-c*p-2*g*j-g-h + j^2*p + l");
     RingElem q_36 = ReadExpr(P9h, "-a*g*j-a*h + c*g^2-c*g*v-2*c + g*h*j + g*j*m + h^2-h*l + j^2-j");
     RingElem q_37 = ReadExpr(P9h, "-a*g*j-c*j + g^2*j + g^2 + g*h + g*j*l-g*v + j^2 + n-1");
     RingElem q_38 = ReadExpr(P9h, "-a*j^2-c*j*v + g*j^2 + g*j*n + h*j-h*n + j^2*v-j*v");
     RingElem q_39 = ReadExpr(P9h, "a*g*p-a*m*p-g*l*p-2*g + h*m*p + h-l + 2*m");
     RingElem q_40 = ReadExpr(P9h, "a*h-a*m + c*g*p + c + g*m-h^2 + h*j*p + h*l + h*m-j*l*p-m^2 + m*v-n*p*v + n");
     RingElem q_41 = ReadExpr(P9h, "-c*m*p-c + g*j*p + g*p + j-m*n*p-n + 1");
     RingElem q_42 = ReadExpr(P9h, "c*h-c*m + g*n-h*l*m + j^2*p-j*l + j*v-2*m*n-n*v");
     RingElem q_43 = ReadExpr(P9h, "a*g-a*l + g^2-g*l-g*m-g*v-h*l + h*n*p + j*p*v-j + l^2-l*n*p-m*p-1");
     RingElem q_44 = ReadExpr(P9h, "c*g-c*l + g*j + h^2*l + h*j-h*l^2 + h*n + h-j*l + j*v-l*n-m*n-m-n*v");
     RingElem q_45 = ReadExpr(P9h, "g*h*l + 2*g*j + g + h*n-j*m + j*v-l-n^2*p-n*v");
     RingElem q_46 = ReadExpr(P9h, "g*n*v + h*j*l-h*l*n + h*n*v + 2*j^2-j*l*v-j*m*v + j*v^2-2*n^2-n*v^2");
     RingElem q_47 = ReadExpr(P9h, "a^2*p-a*h*p-3*a-g^2*p-2*g + h*p*v + v");
     RingElem q_48 = ReadExpr(P9h, "a*c*p + a*h + c*m*p + g*h-g*j*p + h*m + h*v-j*p*v");
     RingElem q_49 = ReadExpr(P9h, "a*g + a*p-c*h*p + c*l*p-g^2 + g*h + g*m-g*v-h^2 + h*l + h*v-j*m*p-j-1");
     RingElem q_50 = ReadExpr(P9h, "a*j + c*n*p + c*v-g*h*m-g*j-h^2*v + h*j + h*n + h*v^2-j*m-2*j*v");
     RingElem q_51 = ReadExpr(P9h, "a*g + a*p-g*h + g*l-g*v-h*p-j*l*p + j*p*v-2*j-2");
     RingElem q_52 = ReadExpr(P9h, "c*g + c*p-g*h*l + h-j*l + j*m");
     RingElem q_53 = ReadExpr(P9h, "g*j + g*n + 2*g-h*j + j*l-j*n*p + p");
     RingElem q_54 = ReadExpr(P9h, "-g*h*n-h*j*v + j^2 + j");
     RingElem q_55 = ReadExpr(P9h, "-a*h*p-a + h*p*v + m^2*p + 2*m-p*v^2 + 3*v");
     RingElem q_56 = ReadExpr(P9h, "-g*h*p + g*p*v-2*g-h + l*m*p + l-m*p*v + 2*m");
     RingElem q_57 = ReadExpr(P9h, "-c*h*p-c-h*j*p-h*m-h*v + j*p*v-2*j + m*n*p + n");
     RingElem q_58 = ReadExpr(P9h, "a*l*p + a-g^2*p-2*g-l*p*v + p*v^2-3*v");
     RingElem q_59 = ReadExpr(P9h, "c*l*p + c-g*j*p-g*v-h*p + h*v-j-l*v + m*n*p + m*v + n-1");
     RingElem q_60 = ReadExpr(P9h, "-g*j*p + g*l-j + l*n*p + l*p + l*v-n*p*v + 2*n + 1");
     RingElem q_61 = ReadExpr(P9h, "-h*n-j^2*p + j*l-j*v + n^2*p + n*v");

     std::vector<RingElem> hia;
     hia.push_back(q_0);
     hia.push_back(q_1);
     hia.push_back(q_2);
     hia.push_back(q_3);
     hia.push_back(q_4);
     hia.push_back(q_5);
     hia.push_back(q_6);
     hia.push_back(q_7);
     hia.push_back(q_8);
     hia.push_back(q_9);
     hia.push_back(q_10);
     hia.push_back(q_11);
     hia.push_back(q_12);
     hia.push_back(q_13);
     hia.push_back(q_14);
     hia.push_back(q_15);
     hia.push_back(q_16);
     hia.push_back(q_17);
     hia.push_back(q_18);
     hia.push_back(q_19);
     hia.push_back(q_20);
     hia.push_back(q_21);
     hia.push_back(q_22);
     hia.push_back(q_23);
     hia.push_back(q_24);
     hia.push_back(q_25);
     hia.push_back(q_26);
     hia.push_back(q_27);
     hia.push_back(q_28);
     hia.push_back(q_29);
     hia.push_back(q_30);
     hia.push_back(q_31);
     hia.push_back(q_32);
     hia.push_back(q_33);
     hia.push_back(q_34);
     hia.push_back(q_35);
     hia.push_back(q_36);
     hia.push_back(q_37);
     hia.push_back(q_38);
     hia.push_back(q_39);
     hia.push_back(q_40);
     hia.push_back(q_41);
     hia.push_back(q_42);
     hia.push_back(q_43);
     hia.push_back(q_44);
     hia.push_back(q_45);
     hia.push_back(q_46);
     hia.push_back(q_47);
     hia.push_back(q_48);
     hia.push_back(q_49);
     hia.push_back(q_50);
     hia.push_back(q_51);
     hia.push_back(q_52);
     hia.push_back(q_53);
     hia.push_back(q_54);
     hia.push_back(q_55);
     hia.push_back(q_56);
     hia.push_back(q_57);
     hia.push_back(q_58);
     hia.push_back(q_59);
     hia.push_back(q_60);
     hia.push_back(q_61);
     
     ideal Ih_6 = ideal(hia);
     examp.push_back(Ih_6); // case 22
     
     /*
     //timing cmd
     ideal rI_6 = radical(Ih_6);
     cout << "hiatarinta = " << rI_6 << endl;
     */

     
     // hf855
     // works maybe
     SparsePolyRing P10 = NewPolyRing(QQ, symbols("X8, X7, X6, X5, X4, X3, X2, x8, x7, x6, x5, x4, x3, x2, h"));

     RingElem y_1 = ReadExpr(P10, "h*x2*x4*X3 + h*x2*x5*X3 + h*x2*x5*X4 + h*x2*x6*X3 + h*x2*x6*X4 + h*x2*x6*X5 + h*x2*x7*X3 + h*x2*x7*X4 + h*x2*x7*X5 + h*x2*x7*X6 + h*x2*x8*X3 + h*x2*x8*X4 + h*x2*x8*X5 + h*x2*x8*X6 + h*x2*x8*X7 + h*x3*x5*X4 + h*x3*x6*X4 + h*x3*x6*X5 + h*x3*x7*X4 + h*x3*x7*X5 + h*x3*x7*X6 + h*x3*x8*X4 + h*x3*x8*X5 + h*x3*x8*X6 + h*x3*x8*X7-h*x3*X2*X4-h*x3*X2*X5-h*x3*X2*X6-h*x3*X2*X7-h*x3*X2*X8 + h*x4*x6*X5 + h*x4*x7*X5 + h*x4*x7*X6 + h*x4*x8*X5 + h*x4*x8*X6 + h*x4*x8*X7-h*x4*X2*X5-h*x4*X2*X6-h*x4*X2*X7-h*x4*X2*X8-h*x4*X3*X5-h*x4*X3*X6-h*x4*X3*X7-h*x4*X3*X8 + h*x5*x7*X6 + h*x5*x8*X6 + h*x5*x8*X7-h*x5*X2*X6-h*x5*X2*X7-h*x5*X2*X8-h*x5*X3*X6-h*x5*X3*X7-h*x5*X3*X8-h*x5*X4*X6-h*x5*X4*X7-h*x5*X4*X8 + h*x6*x8*X7-h*x6*X2*X7-h*x6*X2*X8-h*x6*X3*X7-h*x6*X3*X8-h*x6*X4*X7-h*x6*X4*X8-h*x6*X5*X7-h*x6*X5*X8-h*x7*X2*X8-h*x7*X3*X8-h*x7*X4*X8-h*x7*X5*X8-h*x7*X6*X8-2*x2*x4*X3*X5-2*x2*x4*X3*X6-2*x2*x4*X3*X7-2*x2*x4*X3*X8-2*x2*x5*X3*X6-2*x2*x5*X3*X7-2*x2*x5*X3*X8-2*x2*x5*X4*X6-2*x2*x5*X4*X7-2*x2*x5*X4*X8-2*x2*x6*X3*X7-2*x2*x6*X3*X8-2*x2*x6*X4*X7-2*x2*x6*X4*X8-2*x2*x6*X5*X7-2*x2*x6*X5*X8-2*x2*x7*X3*X8-2*x2*x7*X4*X8-2*x2*x7*X5*X8-2*x2*x7*X6*X8 + 2*x3*x5*X2*X4-2*x3*x5*X4*X6-2*x3*x5*X4*X7-2*x3*x5*X4*X8 + 2*x3*x6*X2*X4 + 2*x3*x6*X2*X5-2*x3*x6*X4*X7-2*x3*x6*X4*X8-2*x3*x6*X5*X7-2*x3*x6*X5*X8 + 2*x3*x7*X2*X4 + 2*x3*x7*X2*X5 + 2*x3*x7*X2*X6-2*x3*x7*X4*X8-2*x3*x7*X5*X8-2*x3*x7*X6*X8 + 2*x3*x8*X2*X4 + 2*x3*x8*X2*X5 + 2*x3*x8*X2*X6 + 2*x3*x8*X2*X7 + 2*x4*x6*X2*X5 + 2*x4*x6*X3*X5-2*x4*x6*X5*X7-2*x4*x6*X5*X8 + 2*x4*x7*X2*X5 + 2*x4*x7*X2*X6 + 2*x4*x7*X3*X5 + 2*x4*x7*X3*X6-2*x4*x7*X5*X8-2*x4*x7*X6*X8 + 2*x4*x8*X2*X5 + 2*x4*x8*X2*X6 + 2*x4*x8*X2*X7 + 2*x4*x8*X3*X5 + 2*x4*x8*X3*X6 + 2*x4*x8*X3*X7 + 2*x5*x7*X2*X6 + 2*x5*x7*X3*X6 + 2*x5*x7*X4*X6-2*x5*x7*X6*X8 + 2*x5*x8*X2*X6 + 2*x5*x8*X2*X7 + 2*x5*x8*X3*X6 + 2*x5*x8*X3*X7 + 2*x5*x8*X4*X6 + 2*x5*x8*X4*X7 + 2*x6*x8*X2*X7 + 2*x6*x8*X3*X7 + 2*x6*x8*X4*X7 + 2*x6*x8*X5*X7");

     RingElem y_2 = ReadExpr(P10, "-h*x2*x4*X3-h*x2*x5*X3-h*x2*x5*X4-h*x2*x6*X3-h*x2*x6*X4-h*x2*x6*X5-h*x2*x7*X3-h*x2*x7*X4-h*x2*x7*X5-h*x2*x7*X6-h*x2*x8*X3-h*x2*x8*X4-h*x2*x8*X5-h*x2*x8*X6-h*x2*x8*X7-h*x3*x5*X4-h*x3*x6*X4-h*x3*x6*X5-h*x3*x7*X4-h*x3*x7*X5-h*x3*x7*X6-h*x3*x8*X4-h*x3*x8*X5-h*x3*x8*X6-h*x3*x8*X7 + h*x3*X2*X4 + h*x3*X2*X5 + h*x3*X2*X6 + h*x3*X2*X7 + h*x3*X2*X8-h*x4*x6*X5-h*x4*x7*X5-h*x4*x7*X6-h*x4*x8*X5-h*x4*x8*X6-h*x4*x8*X7 + h*x4*X2*X5 + h*x4*X2*X6 + h*x4*X2*X7 + h*x4*X2*X8 + h*x4*X3*X5 + h*x4*X3*X6 + h*x4*X3*X7 + h*x4*X3*X8-h*x5*x7*X6-h*x5*x8*X6-h*x5*x8*X7 + h*x5*X2*X6 + h*x5*X2*X7 + h*x5*X2*X8 + h*x5*X3*X6 + h*x5*X3*X7 + h*x5*X3*X8 + h*x5*X4*X6 + h*x5*X4*X7 + h*x5*X4*X8-h*x6*x8*X7 + h*x6*X2*X7 + h*x6*X2*X8 + h*x6*X3*X7 + h*x6*X3*X8 + h*x6*X4*X7 + h*x6*X4*X8 + h*x6*X5*X7 + h*x6*X5*X8 + h*x7*X2*X8 + h*x7*X3*X8 + h*x7*X4*X8 + h*x7*X5*X8 + h*x7*X6*X8 + 2*x2*x4*X3*X5 + 2*x2*x4*X3*X6 + 2*x2*x4*X3*X7 + 2*x2*x4*X3*X8 + 2*x2*x5*X3*X6 + 2*x2*x5*X3*X7 + 2*x2*x5*X3*X8 + 2*x2*x5*X4*X6 + 2*x2*x5*X4*X7 + 2*x2*x5*X4*X8 + 2*x2*x6*X3*X7 + 2*x2*x6*X3*X8 + 2*x2*x6*X4*X7 + 2*x2*x6*X4*X8 + 2*x2*x6*X5*X7 + 2*x2*x6*X5*X8 + 2*x2*x7*X3*X8 + 2*x2*x7*X4*X8 + 2*x2*x7*X5*X8 + 2*x2*x7*X6*X8-2*x3*x5*X2*X4 + 2*x3*x5*X4*X6 + 2*x3*x5*X4*X7 + 2*x3*x5*X4*X8-2*x3*x6*X2*X4-2*x3*x6*X2*X5 + 2*x3*x6*X4*X7 + 2*x3*x6*X4*X8 + 2*x3*x6*X5*X7 + 2*x3*x6*X5*X8-2*x3*x7*X2*X4-2*x3*x7*X2*X5-2*x3*x7*X2*X6 + 2*x3*x7*X4*X8 + 2*x3*x7*X5*X8 + 2*x3*x7*X6*X8-2*x3*x8*X2*X4-2*x3*x8*X2*X5-2*x3*x8*X2*X6-2*x3*x8*X2*X7-2*x4*x6*X2*X5-2*x4*x6*X3*X5 + 2*x4*x6*X5*X7 + 2*x4*x6*X5*X8-2*x4*x7*X2*X5-2*x4*x7*X2*X6-2*x4*x7*X3*X5-2*x4*x7*X3*X6 + 2*x4*x7*X5*X8 + 2*x4*x7*X6*X8-2*x4*x8*X2*X5-2*x4*x8*X2*X6-2*x4*x8*X2*X7-2*x4*x8*X3*X5-2*x4*x8*X3*X6-2*x4*x8*X3*X7-2*x5*x7*X2*X6-2*x5*x7*X3*X6-2*x5*x7*X4*X6 + 2*x5*x7*X6*X8-2*x5*x8*X2*X6-2*x5*x8*X2*X7-2*x5*x8*X3*X6-2*x5*x8*X3*X7-2*x5*x8*X4*X6-2*x5*x8*X4*X7-2*x6*x8*X2*X7-2*x6*x8*X3*X7-2*x6*x8*X4*X7-2*x6*x8*X5*X7");

     RingElem y_3 = ReadExpr(P10, "15*h^3 + 26*h^2*x2 + 26*h^2*x3 + 26*h^2*x4 + 26*h^2*x5 + 26*h^2*x6 + 8*h*x3*X2 + 8*h*x4*X2 + 8*h*x4*X3 + 8*h*x5*X2 + 8*h*x5*X3 + 8*h*x5*X4 + 8*h*x6*X2 + 8*h*x6*X3 + 8*h*x6*X4 + 8*h*x6*X5 + 16*x2*x4*X3 + 16*x2*x5*X3 + 16*x2*x5*X4 + 16*x2*x6*X3 + 16*x2*x6*X4 + 16*x2*x6*X5 + 16*x3*x5*X4 + 16*x3*x6*X4 + 16*x3*x6*X5 + 16*x4*x6*X5");

     RingElem y_4 = ReadExpr(P10, "15*h^3 + 26*h^2*X2 + 26*h^2*X3 + 26*h^2*X4 + 26*h^2*X5 + 26*h^2*X6 + 8*h*x2*X3 + 8*h*x2*X4 + 8*h*x2*X5 + 8*h*x2*X6 + 8*h*x3*X4 + 8*h*x3*X5 + 8*h*x3*X6 + 8*h*x4*X5 + 8*h*x4*X6 + 8*h*x5*X6 + 16*x3*X2*X4 + 16*x3*X2*X5 + 16*x3*X2*X6 + 16*x4*X2*X5 + 16*x4*X2*X6 + 16*x4*X3*X5 + 16*x4*X3*X6 + 16*x5*X2*X6 + 16*x5*X3*X6 + 16*x5*X4*X6");

     RingElem y_5 = ReadExpr(P10, "h^2 + 2*h*x2 + 2*h*x3 + 2*h*x4 + 2*h*x5 + 2*h*x6 + 2*h*x7 + 2*h*x8-4*x2*X3-4*x2*X4-4*x2*X5-4*x2*X6-4*x2*X7-4*x2*X8 + 4*x3*X2-4*x3*X4-4*x3*X5-4*x3*X6-4*x3*X7-4*x3*X8 + 4*x4*X2 + 4*x4*X3-4*x4*X5-4*x4*X6-4*x4*X7-4*x4*X8 + 4*x5*X2 + 4*x5*X3 + 4*x5*X4-4*x5*X6-4*x5*X7-4*x5*X8 + 4*x6*X2 + 4*x6*X3 + 4*x6*X4 + 4*x6*X5-4*x6*X7-4*x6*X8 + 4*x7*X2 + 4*x7*X3 + 4*x7*X4 + 4*x7*X5 + 4*x7*X6-4*x7*X8 + 4*x8*X2 + 4*x8*X3 + 4*x8*X4 + 4*x8*X5 + 4*x8*X6 + 4*x8*X7");

     RingElem y_6 = ReadExpr(P10, "h^2 + 2*h*X2 + 2*h*X3 + 2*h*X4 + 2*h*X5 + 2*h*X6 + 2*h*X7 + 2*h*X8 + 4*x2*X3 + 4*x2*X4 + 4*x2*X5 + 4*x2*X6 + 4*x2*X7 + 4*x2*X8-4*x3*X2 + 4*x3*X4 + 4*x3*X5 + 4*x3*X6 + 4*x3*X7 + 4*x3*X8-4*x4*X2-4*x4*X3 + 4*x4*X5 + 4*x4*X6 + 4*x4*X7 + 4*x4*X8-4*x5*X2-4*x5*X3-4*x5*X4 + 4*x5*X6 + 4*x5*X7 + 4*x5*X8-4*x6*X2-4*x6*X3-4*x6*X4-4*x6*X5 + 4*x6*X7 + 4*x6*X8-4*x7*X2-4*x7*X3-4*x7*X4-4*x7*X5-4*x7*X6 + 4*x7*X8-4*x8*X2-4*x8*X3-4*x8*X4-4*x8*X5-4*x8*X6-4*x8*X7");

     RingElem y_7 = ReadExpr(P10, "-h^2 + x8*X8");     
     RingElem y_8 = ReadExpr(P10, "-h^2 + x7*X7");
     RingElem y_9 = ReadExpr(P10, "-h^2 + x6*X6");
     RingElem y_10 = ReadExpr(P10, "-h^2 + x5*X5");
     RingElem y_11 = ReadExpr(P10, "-h^2 + x4*X4");
     RingElem y_12 = ReadExpr(P10, "-h^2 + x3*X3");
     RingElem y_13 = ReadExpr(P10, "-h^2 + x2*X2");
     RingElem y_14 = ReadExpr(P10, "h + 2*X2 + 2*X3 + 2*X4 + 2*X5 + 2*X6 + 2*X7 + 2*X8");
     RingElem y_15 = ReadExpr(P10, "h + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 + 2*x8");

     std::vector<RingElem> yons;
     yons.push_back(y_1);
     yons.push_back(y_2);
     yons.push_back(y_3);
     yons.push_back(y_4);
     yons.push_back(y_5);
     yons.push_back(y_6);
     yons.push_back(y_7);
     yons.push_back(y_8);
     yons.push_back(y_9);
     yons.push_back(y_10);
     yons.push_back(y_11);
     yons.push_back(y_12);
     yons.push_back(y_13);
     yons.push_back(y_14);
     yons.push_back(y_15);
     
     ideal Y = ideal(yons);

     examp.push_back(Y); // case 23


     //Examples from paper: 
     //positive characteristic:
     //ring ZZMod101 = NewZZmod(101);
     ring ZZMod23 = NewZZmod(23);
     SparsePolyRing P8 = NewPolyRing(ZZMod101, symbols("x, y, z"));

     ideal L6 = ideal(RingElems(P8, "x,y,z"));

     //CoCoA_ASSERT_ALWAYS(radical(L6) == L6);
     examp.push_back(L6); // case 24

     
     //Characteristic is 0:

     ideal H2 = ideal(ReadExpr(P4, "x^4 +83*x^3 +73*y^2 -85*z^2 -437*t"), ReadExpr(P4, "y^3-x"), ReadExpr(P4, "z^3 +z -t"), ReadExpr(P4, "t^3-324*z^2 +94*y^2 +76*x"));

     // std::vector<RingElem> resH2;
     // resH2.push_back(ReadExpr(P4, "t^3 +94*y^2 -324*z^2 +76*x"));
     // resH2.push_back(ReadExpr(P4, "z^3 +z -t"));
     // resH2.push_back(ReadExpr(P4, "y^3 -x"));
     // resH2.push_back(ReadExpr(P4, "x^4 +83*x^3 +73*y^2 -85*z^2 -437*t"));
     // resH2.push_back(ReadExpr(P4, "139623900223990351168240621*x^3*y^2*z^2*t^2 +11588783718591199146963971543*x^2*y^2*z^2*t^2 +139623900223990351168240621*x^3*y^2*t^2 -45238143672572873778509961204*x^3*y^2*z +11588783718591199146963971543*x^2*y^2*t^2 -3754765924823548523616326779932*x^2*y^2*z +10192544716351295635281565333*y*z^2*t^2 +4637188974239167542999607504652*y^2*z^2 +901970395446977668546834411660*y^2*z +5735470573401075645288988229438*y*z^2 +10192544716351295635281565333*y*t^2 +4637188974239167542999607504652*y^2 -2186789525308136879996984606102*y*z +5735470573401075645288988229438*y"));
    
     //CoCoA_ASSERT_ALWAYS(radical(H2) == ideal(resH2));
     examp.push_back(H2); // case 25

     
     //examples John

     ideal W1 = ideal(RingElems(P3, "x^2 +x*y +z^2,  y^2 -y*z -x,  y^2 +z^2 -x"));
     ideal W2 = ideal(RingElems(P3, "y^2 +x*z +y*z,  x^2 -y^2 +x*z, z^2 +x +y"));

     // std::vector<RingElem> resW1;
     // resW1.push_back(ReadExpr(P3, "y*z +z^2"));
     // resW1.push_back(ReadExpr(P3, "x*z +(-1/2)*x +(-1/2)*y"));
     // resW1.push_back(ReadExpr(P3, "y^2 +z^2 -x"));
     // resW1.push_back(ReadExpr(P3, "x*y -2*z^2 +(3/2)*x +(1/2)*y"));
     // resW1.push_back(ReadExpr(P3, "x^2 +3*z^2 +(-3/2)*x +(-1/2)*y"));
     // resW1.push_back(ReadExpr(P3, "z^3 +(-1/4)*x +(-1/4)*y"));
     // resW1.push_back(ReadExpr(P3, "(-1/2)*z^2 +(1/4)*x +(1/4)*y +(1/4)*z"));
   
     //CoCoA_ASSERT_ALWAYS(radical(W1) == ideal(resW1));
     examp.push_back(W1); // case 26

     // std::vector<RingElem> resW2;
     // resW2.push_back(ReadExpr(P3, "z^2 +x +y"));
     // resW2.push_back(ReadExpr(P3, "x*z -y*z -x"));
     // resW2.push_back(ReadExpr(P3, "y^2 +2*y*z +x"));
     // resW2.push_back(ReadExpr(P3, "x*y +3*y*z -y"));
     // resW2.push_back(ReadExpr(P3, "x^2 +3*y*z +2*x"));
     // resW2.push_back(ReadExpr(P3, "y*z -2*x -3*y -z"));
   
     //CoCoA_ASSERT_ALWAYS(radical(W2) == ideal(resW2));
     examp.push_back(W2); // case 27

     
     /*
              |0dim?| runs?
     case 0:  | yes | yes
     case 1:  | yes | yes
     case 2:  | yes | yes
     case 3:  | yes | yes
     case 4:  | yes | yes
     case 5:  | yes | yes
     case 6:  | yes | yes
     case 7:  | yes | yes
     case 8:  | yes | yes
     case 9:  | yes | yes
     case 10: | yes | yes
     case 11: | no  | yes
     case 12: | no  | yes
     case 13: | no  | yes
     case 14: | no  | yes
     case 15: | yes | yes
     case 16: | no  | yes
     case 17: | no  | yes
     case 18: | no  | yes
     case 19: | no  | yes
     case 20: | no  | yes
     case 21: | no  | yes
     case 22: | no  | yes
     case 23: | no  | no  (it does a lot of calculations, but didn't stop after an hour)
     case 24: | yes | yes
     case 25: | yes | yes
     case 26: | yes | yes
     case 27: | yes | yes

     Comments:
     Cases that don't terminate successfully yet: 23
     IsRadEqToRadOfRad always returns true (though in case 23 we can't know yet of course)
     */

  }

  
  void CallExample(const std::vector<ideal>& examp, const int EX)
  {
    VerboseLog VERBOSE("CallExample");
    ring QQ = RingQQ();
    
    std::string L6_str = "x,y,z";
     
    std::string resH2_str =
       "t^3 +94*y^2 -324*z^2 +76*x,"
       "z^3 +z -t,"
       "y^3 -x,"
       "x^4 +83*x^3 +73*y^2 -85*z^2 -437*t,"
       "139623900223990351168240621*x^3*y^2*z^2*t^2 +11588783718591199146963971543*x^2*y^2*z^2*t^2 +139623900223990351168240621*x^3*y^2*t^2 -45238143672572873778509961204*x^3*y^2*z +11588783718591199146963971543*x^2*y^2*t^2 -3754765924823548523616326779932*x^2*y^2*z +10192544716351295635281565333*y*z^2*t^2 +4637188974239167542999607504652*y^2*z^2 +901970395446977668546834411660*y^2*z +5735470573401075645288988229438*y*z^2 +10192544716351295635281565333*y*t^2 +4637188974239167542999607504652*y^2 -2186789525308136879996984606102*y*z +5735470573401075645288988229438*y";

     std::string resW1_str =
      "y*z +z^2,"
      "x*z +(-1/2)*x +(-1/2)*y,"
      "y^2 +z^2 -x,"
      "x*y -2*z^2 +(3/2)*x +(1/2)*y,"
      "x^2 +3*z^2 +(-3/2)*x +(-1/2)*y,"
      "z^3 +(-1/4)*x +(-1/4)*y,"
      "(-1/2)*z^2 +(1/4)*x +(1/4)*y +(1/4)*z";

     std::string resW2_str =
      "z^2 +x +y,"
      "x*z -y*z -x,"
      "y^2 +2*y*z +x,"
      "x*y +3*y*z -y,"
      "x^2 +3*y*z +2*x,"
      "y*z -2*x -3*y -z";

     ideal RadI(zero(RingQQ()));
               
     switch(EX)
     {
     case  0: //terminates
     case  1: //terminates
     case  2: //rad(I)=rad(rad(I)), but computation slow
     case  3: //terminates
     case  4: //terminates
     case  5: //terminates
     case  6: //terminates
     case  7: //terminates
     case  8: //terminates
     case  9: //terminates
     case 10: //terminates
     case 11: //terminates
     case 12: //terminates
     case 13: //terminates
     case 14: //terminates
     case 15: //terminates
     case 16: //terminates
     case 17: //terminates
     case 18: //terminates
     case 19: //terminates
     case 20: //terminates
     case 21: //terminates
     case 22: //terminates
     case 23: //IsZeroDim too slow
     case 24: //terminates
       TestRadical(RadI, examp[EX]); break;
     case 25: TestRadical(RadI, examp[EX]); //terminates
       VERBOSE(20) << "Calculation correct? "
            << (radical(examp[EX]) == ideal(RingElems(RingOf(examp[EX]),resH2_str)))
            << endl;
       break;
     case 26: TestRadical(RadI, examp[EX]); //terminates
       VERBOSE(20) << "Calculation correct? "
            << (radical(examp[EX]) == ideal(RingElems(RingOf(examp[EX]),resW1_str)))
            << endl;
       break;
     case 27: TestRadical(RadI, examp[EX]); //terminates
       VERBOSE(20) << "Calculation correct? "
            << (radical(examp[EX]) == ideal(RingElems(RingOf(examp[EX]),resW2_str)))
            << endl;
       break;
     default: VERBOSE(20) << "No example with this index" << endl; break;
     }
  }


  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    cout << std::boolalpha; // so that bools print out as true/false
    std::vector<ideal> examples;
    FillExampleList(examples); 

    //    SetVerbosityLevel(20);
    SetVerbosityLevel(0);
    for (long i=1; i<27; ++i)
    {
      cout << "----- " << i << " -----------" << endl;
      if (i==23)  cout << "** skip **  (IsZeroDim too slow)" << endl;
      else
      {
        if (!IsZeroDim(examples[i]))  CallExample(examples, i);
      }
    }
  }
  
     
     /*
     ideal Yrad = radical(Y);
     cout << "Yrad = " << Yrad << endl;
*/
     
     /*
     
     //Examples from paper: 
     //positive characteristic:
     ring ZZMod101 = NewZZmod(101);
     ring ZZMod23 = NewZZmod(23);
     SparsePolyRing P8 = NewPolyRing(ZZMod101, symbols("x, y, z"));

     ideal L6 = ideal(ReadExpr(P8, "x"), ReadExpr(P8, "y"), ReadExpr(P8, "z"));

     CoCoA_ASSERT_ALWAYS(radical(L6) == L6);
    
     //Characteristic is 0:

     ideal H2 = ideal(ReadExpr(P4, "x^4 +83*x^3 +73*y^2 -85*z^2 -437*t"), ReadExpr(P4, "y^3-x"), ReadExpr(P4, "z^3 +z -t"), ReadExpr(P4, "t^3-324*z^2 +94*y^2 +76*x"));

     std::vector<RingElem> resH2;
     resH2.push_back(ReadExpr(P4, "t^3 +94*y^2 -324*z^2 +76*x"));
     resH2.push_back(ReadExpr(P4, "z^3 +z -t"));
     resH2.push_back(ReadExpr(P4, "y^3 -x"));
     resH2.push_back(ReadExpr(P4, "x^4 +83*x^3 +73*y^2 -85*z^2 -437*t"));
     resH2.push_back(ReadExpr(P4, "139623900223990351168240621*x^3*y^2*z^2*t^2 +11588783718591199146963971543*x^2*y^2*z^2*t^2 +139623900223990351168240621*x^3*y^2*t^2 -45238143672572873778509961204*x^3*y^2*z +11588783718591199146963971543*x^2*y^2*t^2 -3754765924823548523616326779932*x^2*y^2*z +10192544716351295635281565333*y*z^2*t^2 +4637188974239167542999607504652*y^2*z^2 +901970395446977668546834411660*y^2*z +5735470573401075645288988229438*y*z^2 +10192544716351295635281565333*y*t^2 +4637188974239167542999607504652*y^2 -2186789525308136879996984606102*y*z +5735470573401075645288988229438*y"));
    
     CoCoA_ASSERT_ALWAYS(radical(H2) == ideal(resH2));
    
     //examples John

     ideal W1 = ideal(ReadExpr(P3,"x^2 +x*y +z^2"), ReadExpr(P3, "y^2 -y*z -x"), ReadExpr(P3, "y^2 +z^2 -x"));
     ideal W2 = ideal(ReadExpr(P3,"y^2 +x*z +y*z"), ReadExpr(P3, "x^2 -y^2 +x*z"), ReadExpr(P3, "z^2 +x +y"));

     std::vector<RingElem> resW1;
     resW1.push_back(ReadExpr(P3, "y*z +z^2"));
     resW1.push_back(ReadExpr(P3, "x*z +(-1/2)*x +(-1/2)*y"));
     resW1.push_back(ReadExpr(P3, "y^2 +z^2 -x"));
     resW1.push_back(ReadExpr(P3, "x*y -2*z^2 +(3/2)*x +(1/2)*y"));
     resW1.push_back(ReadExpr(P3, "x^2 +3*z^2 +(-3/2)*x +(-1/2)*y"));
     resW1.push_back(ReadExpr(P3, "z^3 +(-1/4)*x +(-1/4)*y"));
     resW1.push_back(ReadExpr(P3, "(-1/2)*z^2 +(1/4)*x +(1/4)*y +(1/4)*z"));
   
     CoCoA_ASSERT_ALWAYS(radical(W1) == ideal(resW1));

     std::vector<RingElem> resW2;
     resW2.push_back(ReadExpr(P3, "z^2 +x +y"));
     resW2.push_back(ReadExpr(P3, "x*z -y*z -x"));
     resW2.push_back(ReadExpr(P3, "y^2 +2*y*z +x"));
     resW2.push_back(ReadExpr(P3, "x*y +3*y*z -y"));
     resW2.push_back(ReadExpr(P3, "x^2 +3*y*z +2*x"));
     resW2.push_back(ReadExpr(P3, "y*z -2*x -3*y -z"));
   
     CoCoA_ASSERT_ALWAYS(radical(W2) == ideal(resW2));
     */

  

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
