// Copyright (c) 2005  John Abbott,  Anna M Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program showing `AlexanderDual' and `PrimaryDecomposition' on \n"
  "monomial ideals.  If Frobby is available, it is used too.     \n";

const string LongDescription =
  "This example shows how to compute the `AlexanderDual' and the       \n"
  "`PrimaryDecomposition' of monomial ideals.  These operations are    \n"
  "offered by both CoCoALib and the external library Frobby.  This     \n"
  "program shows how to check whether Frobby is available, and if so,  \n"
  "how to call its functions with CoCoALib data.                       \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring Fp = NewZZmod(32003);          // coefficient ring
    SparsePolyRing Fpx = NewPolyRing(Fp, SymbolRange("x",0,7)); // Fp[x[0..7]]
    SparsePolyRing P = Fpx;
    double t0;  // for CpuTime
  
    const vector<RingElem>& x = indets(P);
    vector<RingElem> g;
    g = RingElems(P, 
                  "x[2]^2 * x[5]^4,"
                  "x[1]^3 * x[4]^4,"
                  "x[1]^3 * x[5]^4,"
                  "x[3]^3 * x[6]^4,"
                  "x[3]^4 * x[6]^3");

    ideal J1(g);
    ideal J2(x[1]*x[2], x[2]*x[3]*x[5], x[0]*x[1]*x[3], x[0]*x[3]*x[5]);
    cout << "J1  = " << J1 << endl; 
    cout << "J2  = " << J2 << endl << endl;

    t0 = CpuTime();
    ideal I = intersect(J1, J2);
    cout << "Cpu Time = " << CpuTime()-t0 << endl;
    cout << "intersect(J1, J2) = " << I << endl;
    cout << endl;

    cout << "Only for squarefree monomial ideals:" << endl;

    t0 = CpuTime();
    ideal AD = AlexanderDual(J2);
    cout << "Cpu Time = " << CpuTime()-t0 << endl;
    cout << "AlexanderDual(J2) = " << AD << endl;

    t0 = CpuTime();
    vector<ideal> PrimDec = PrimaryDecomposition(J2);
    cout << "Cpu Time = " << CpuTime()-t0 << endl;
    cout << "PrimaryDecomposition(J2) = " << PrimDec << endl;

    cout << endl;

#ifndef CoCoA_WITH_FROBBY
    cout << "External library Frobby is not available, so we skip the Frobby examples." << endl << endl;
#else
    cout << "Frobby can work on any monomial ideal:" << endl;
    t0 = CpuTime();
    AD = FrbAlexanderDual(J2);
    cout << "Cpu Time = " << CpuTime()-t0 << endl;
    cout << "FrbAlexanderDual(J2) = " << AD << endl;

    PrimDec.empty();
    t0 = CpuTime();
    FrbPrimaryDecomposition(PrimDec, J2);
    cout << "Cpu Time = " << CpuTime()-t0 << endl;
    cout << "FrbPrimaryDecomposition(PrimDec, J2) => " << PrimDec << endl;

    t0 = CpuTime();
    AD = FrbAlexanderDual(J1);
    cout << "Cpu Time = " << CpuTime()-t0 << endl;
    cout << "FrbAlexanderDual(J1) = " << AD << endl;

    t0 = CpuTime();
    PPMonoidElem pp = power(product(indets(PPM(P))), 6);
    AD = FrbAlexanderDual(J1, pp);
    cout << "Cpu Time = " << CpuTime()-t0 << endl;
    cout << "pp = " << pp << endl;
    cout << "FrbAlexanderDual(J1, pp) = " << AD << endl;

    t0 = CpuTime();
    FrbPrimaryDecomposition(PrimDec, J1);
    cout << "Cpu Time = " << CpuTime()-t0 << endl;
    cout << "FrbPrimaryDecomposition(PrimDec, J1) => " << PrimDec << endl;
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
