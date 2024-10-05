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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-AlexanderDual.C,v 1.19 2022/02/13 09:56:55 abbott Exp $
// $Log: ex-AlexanderDual.C,v $
// Revision 1.19  2022/02/13 09:56:55  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.18  2019/10/18 14:12:41  bigatti
// -- updated after having "RingElems"
//
// Revision 1.17  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.16  2017/04/27 15:24:52  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.15  2015/06/25 16:40:27  abbott
// Summary: Moved code into namespace CoCoA (redmine 739)
// Author: JAA
//
// Revision 1.14  2014/03/21 16:44:32  bigatti
// -- improved input using ReadExpr
//
// Revision 1.13  2013/06/28 12:06:28  abbott
// Updated Frobby fn names.
//
// Revision 1.12  2012/09/21 13:37:44  abbott
// Improved descriptions; prints out clearer message when Frobby is not present.
//
// Revision 1.11  2012/02/08 17:40:23  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.10  2011/09/06 13:35:35  abbott
// Minor layout change (in a string).
//
// Revision 1.9  2011/07/05 16:25:14  bigatti
// -- reorganized since now AlexanderDual is part of cocoalib
//
// Revision 1.8  2011/07/05 15:36:57  bigatti
// -- modified example: now AlexanderDual is in CoCoALib
//
// Revision 1.7  2011/05/24 14:58:25  abbott
// Consequential changes from removal of several ctors for principal ideals.
//
// Revision 1.6  2011/03/10 17:15:50  abbott
// Replaced size_t by long.
//
// Revision 1.5  2010/12/26 13:02:03  abbott
// Finished changing "GlobalXXXput()" into the corresponding standard C++ stream.
//
// Revision 1.4  2010/12/17 16:07:55  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2010/11/26 15:21:58  bigatti
// -- added "Example" to function names for collision with Frobby functions
//
// Revision 1.2  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/02/28 13:51:59  bigatti
// -- added function IsMonomial
//
// Revision 1.6  2007/02/26 15:49:08  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.5  2007/02/12 15:29:07  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.4  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.3  2006/11/21 10:57:35  cocoa
// -- cleaned up ex-AlexanderDual.C and added to Makefile
//
// Revision 1.2  2006/11/17 18:14:55  cocoa
// -- added const & to some arguments (embarrassing!)
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/04/28 11:26:06  cocoa
// -- first import
//
