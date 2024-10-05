// Copyright (c) 2005  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "The example in this file shows how to create and use some \n"
  "homomorphisms between rings.                              \n";

const string LongDescription =
  "We compute these polynomials (with parameters) in some rings: \n"
  "f = (2*a/3-1)*x[0] + 1/a;  g = x[0]-a;                        \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  // A simple function to test embedding ring homs in CoCoALib.
  // We know that Kx is of the form R(a)[x[0..N]]
  void SimpleTest(const PolyRing& Kx) 
  {
    cout << "Kx is " << Kx << endl;
    // Give a handy name to the indeterminate x[0]
    RingElem x0 = indet(Kx, 0);
    cout << "indet(Kx,0) is " << indet(Kx,0) << endl;

    // the following 6 lines show how to work with several rings and RingHoms

    ring K  = CoeffRing(Kx); // this is R(a)
    ring Ra = BaseRing(K);  // this is R[a]
    RingHom RaToK = CanonicalHom(Ra, K);     // K = FractionField(R[a])
    RingHom KToKx = CanonicalHom(K, Kx);
    RingElem a_R = indet(Ra, 0); // a as RingElem of R[a]
    RingElem a   = KToKx(RaToK(a_R));        // a as RingElem of Kx

    // for this particular example this would have been much simpler:
    //   RingElem a(Kx, symbol("a"));

    RingElem f = (2*a/3-1)*x0 + 1/a;         // This is a RingElem of Kx.
    RingElem g = x0-a;                       // This is another one.

    cout << "f = " << f << endl;
    cout << "f*g = " << f*g << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // this changes "cout" so that bools are printed true/false

    ring QQ = RingQQ();
    ring Fp = NewZZmod(32003);

    //---------------------------------------------------------------------------
    {
      cout << "  --- Coeffs in " << Fp << endl;
      PolyRing Fpa = NewPolyRing(Fp, symbols("a"));              // Fpa is Fp[a]  where Fp = ZZ/(32003)
      ring K = NewFractionField(Fpa);                            // K   is Fp(a)
      PolyRing Kxyzt = NewPolyRing(K, symbols("x,y,z,t")); // Kxyzt is Fp(a)[x,y,z,t]
      SimpleTest(Kxyzt);
      cout << endl;
    }

    //---------------------------------------------------------------------------
    {
      cout << "Coeffs in " << RingZZ() << endl;
      PolyRing ZZa = NewPolyRing(RingZZ(), symbols("a"));    // ZZ[a]
      PolyRing Qax4 = NewPolyRing(NewFractionField(ZZa), NewSymbols(4)); // QQ(a)[x[0..3]]
      SimpleTest(Qax4);
      cout << endl << endl;
    }

    //---------------------------------------------------------------------------
    {
      cout << "  --- Coeffs in " << QQ << endl;
      PolyRing Qa = NewPolyRing(QQ, symbols("a"));           // QQ[a]
      PolyRing Qax4 = NewPolyRing(NewFractionField(Qa), NewSymbols(4));  // QQ(a)[x[0..3]]
      SimpleTest(Qax4);
      cout << endl;
    }

    cout << "  --- Now we supply an unsuitable input" << endl;
    cout << "      (it is unsuitable because CoeffRing is not K(a))" << endl;
    try
    {
      SimpleTest(NewPolyRing(Fp, symbols("a")));
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      if (err != ERR::ReqPolyRing) throw;
      cout << "\nOK!  Our unsuitable input produced the expected exception." << endl;
    }
    cout << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingHom5.C,v 1.17 2024/07/11 15:39:31 bigatti Exp $
// $Log: ex-RingHom5.C,v $
// Revision 1.17  2024/07/11 15:39:31  bigatti
// Summary: error codes: NotElemPolyRing to ReqElemPolyRing
//
// Revision 1.16  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.15  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.14  2017/07/03 12:27:26  abbott
// Summary: Minor correction
//
// Revision 1.13  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.12  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.11  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.10  2014/07/22 11:14:48  bigatti
// -- changed message for unsuitable input
//
// Revision 1.9  2014/07/08 15:47:40  abbott
// Summary: Modified "unsuitable" example (not sure I've done the right thing though)
// Author: JAA
//
// Revision 1.8  2014/07/08 12:48:18  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.7  2014/07/08 12:47:26  abbott
// Summary: Removed AsPolyRing, AsSparsePolyRing, AsQuotientRing
// Author: JAA
//
// Revision 1.6  2012/04/04 08:32:02  bigatti
// -- added comment for boolalpha
// -- uncommended example about ZZ(a)  -- now it works
// -- cleaning Q --> QQ
//
// Revision 1.5  2012/02/08 17:45:51  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2010/10/01 15:47:57  bigatti
// -- added a comment about using RingElem(Kx, symbol("a"))
//
// Revision 1.2  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.1  2007/03/07 14:31:01  bigatti
// -- was ex-RingHom1.C
//
// Revision 1.6  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/02/26 17:39:46  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.4  2007/02/12 16:11:12  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.3  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.2  2007/01/18 14:34:18  cocoa
// -- example of alternative syntax
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.3  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.2  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.2  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.3  2004/12/09 15:08:42  cocoa
// -- added log info
//
// Revision 1.2  2004/11/22 17:10:04  cocoa
// -- changed  "PolyRing" --> "AsPolyRing",  "FractionField" --> "AsFractionField"
//
