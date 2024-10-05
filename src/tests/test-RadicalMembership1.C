//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/QuotientRing.H"
#include "CoCoA/SparsePolyOps-ideal-RadicalMembership.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/time.H"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::flush;

#include <vector>
using std::vector;

//----------------------------------------------------------------------
// This test checks that IsInRadical and MinPowerInIdeal work in some
// simple test cases.
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.


  bool IsInRadicalVerbose(ConstRefRingElem f, const ideal& I)
  {
    //    cout << "IsInRadical: f is " << f << endl;
    //    double t = CpuTime();
    bool b = IsInRadical(f,I);
    //    cout << "  is " << b << "\ttime " << CpuTime()-t << endl;
    return b;
  }


  long MinPowerInIdealVerbose(ConstRefRingElem f, const ideal& I)
  {
    //    cout << "MinPowerInIdeal: f is " << f << endl;
    //    double t = CpuTime();
    long m = MinPowerInIdeal(f,I);
    //    cout << "  is " << m << "\ttime " << CpuTime()-t << endl;
    return m;
  }


  void test_MonomialIdeal(const PolyRing& P)
  {
    //cout << "  --====Monomial-ideal====--" << endl;
    ideal I = ideal(RingElems(P, "x^5, y^8, z^13"));
    
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(RingElem(P,"x"), I));
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(RingElem(P,"x^2"), I));
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(RingElem(P,"x^99"), I));
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(RingElem(P,"y"), I));
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(RingElem(P,"y^2"), I));
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(RingElem(P,"y^99"), I));
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(RingElem(P,"z"), I));
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(RingElem(P,"z^2"), I));
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(RingElem(P,"z^99"), I));
    
    RingElem f1 = RingElem(P,"x+y");
    RingElem f2 = RingElem(P,"(x+y+z)^2");
    
    CoCoA_ASSERT_ALWAYS(!IsInRadicalVerbose(f1+1,I));
    CoCoA_ASSERT_ALWAYS(!IsInRadicalVerbose(f2+2,I));
    
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(f1,I));
    CoCoA_ASSERT_ALWAYS(MinPowerInIdealVerbose(f1,I) == 12);
    
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(f2,I));
    CoCoA_ASSERT_ALWAYS(MinPowerInIdealVerbose(f2,I) == 12);
    
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(f1-f2,I));
    CoCoA_ASSERT_ALWAYS(MinPowerInIdealVerbose(f1-f2,I) == 18);
    
    CoCoA_ASSERT_ALWAYS(!IsInRadicalVerbose(f1+1,I));
    CoCoA_ASSERT_ALWAYS(MinPowerInIdealVerbose(f1+1,I) == -1);
    
    CoCoA_ASSERT_ALWAYS(!IsInRadicalVerbose(f2-1,I));
    CoCoA_ASSERT_ALWAYS(MinPowerInIdealVerbose(f2-1,I) == -1);
  }
  

  void test_HomogIdeal(const PolyRing& P)
  {
    //cout << "  --====Homogeneous-ideal====--" << endl;
    RingElem g1 = RingElem(P, "2*x^2 +3*y*z -x*z");
    RingElem g2 = RingElem(P, "3*x*y*z^2 -5*x*z^3 +2*y^4");
    ideal I = ideal(power(g1,4) + power(g2,2),
                    power(g1,5));
    
    RingElem f1 = g1;
    CoCoA_ASSERT_ALWAYS(!IsInRadicalVerbose(f1 + RingElem(P,"x^2"),I));
    // CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(f1,I));  called by MinPower
    CoCoA_ASSERT_ALWAYS(MinPowerInIdealVerbose(f1,I) == 5);
    
    RingElem f2 = g2*g2 + g1 - g2;
    CoCoA_ASSERT_ALWAYS(!IsInRadicalVerbose(f2 + RingElem(P,"y^2"),I));
    // CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(f2,I));  called by MinPower
    CoCoA_ASSERT_ALWAYS(MinPowerInIdealVerbose(f2,I) == 6);
  }

  
  void test_NonHomogIdeal(const PolyRing& P)
  {
    //cout << "  --====Non-homogeneous-ideal====--" << endl;
    RingElem g1 = RingElem(P,"2*x^2 +3*y*z -x");
    RingElem g2 = RingElem(P,"3*x*y*z -x*z +y");
    ideal I = ideal(power(g1,4),  power(g2,3));
    
    RingElem f1 = g1 + g2;
    CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(f1,I));

    // (default: off) slow & non-essential tests for further checks:
    
    // CoCoA_ASSERT_ALWAYS(!IsInRadicalVerbose(f1 + RingElem(P,"x"),I));
    // CoCoA_ASSERT_ALWAYS(IsElem(power(f1,MinPowerInIdealVerbose(f1,I)), I));
    
    // RingElem f2 = g2*g2 + g1 - g2;
    // CoCoA_ASSERT_ALWAYS(IsInRadicalVerbose(f2,I));
    // CoCoA_ASSERT_ALWAYS(!IsInRadicalVerbose(f2 + RingElem(P,"y^2"),I));
    // CoCoA_ASSERT_ALWAYS(IsElem(power(f2,MinPowerInIdealVerbose(f2,I)), I));
  }

  void test_IJ()
  {
    ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    // First test case: monomial, principal ideals
    ideal I(indet(P,0));
    ideal J(indet(P,1));
    CoCoA_ASSERT_ALWAYS(!IsInRadical(I,J));
    CoCoA_ASSERT_ALWAYS(!IsInRadical(J,I));
    // Second test case: non-monomial, non-principal ideals
    ideal II(RingElemVec(P,"[x+1,y+1]"));
    ideal JJ(RingElemVec(P,"[y+1,z+1]"));
    CoCoA_ASSERT_ALWAYS(!IsInRadical(II,JJ));
    CoCoA_ASSERT_ALWAYS(!IsInRadical(JJ,II));
  }
  
  void program()
  {
    GlobalManager CoCoAFoundations;
    
    cout << std::boolalpha; // so that bools print out as true/false

    // Use ZZ/(32003) for coeffs as QQ is noticeably slower
    PolyRing P = NewPolyRing(NewZZmod(32003), symbols("x, y, z"));

    test_MonomialIdeal(P);
    test_HomogIdeal(P);
    test_NonHomogIdeal(P);
    test_IJ(); // 2024-03-03 spotted bug in impl -- test it here!
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
