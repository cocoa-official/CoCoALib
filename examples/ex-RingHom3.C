// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how we define a ring homomorphism       \n"
  "to perform a change of coordinates in a polynomial ring.   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Create coeff ring, then poly ring.
    ring Fp = NewZZmod(32003);
    PolyRing Fpx = NewPolyRing(Fp, symbols("a,b,c,d")); // Fpx is Z/(32003)[a,b,c,d]

    cout << "We shall work in the PolyRing " << Fpx << endl;

    // Line below gives easy access to the indets: they are just x[0], x[1],...
    const vector<RingElem>& x = indets(Fpx);

    // Careful with rationals in C++ because C++ treats 2/3 as an integer division
    const RingElem TwoThirds = RingElem(Fpx,2)/3;

    RingElem f = 2*x[0] + 4*x[1]*x[1] - TwoThirds*x[0]*x[2];
    cout << "Original poly is f = " << f << endl;

    // Now create the poly algebra homomorphism; start with vector of images:
    vector<RingElem> images;
    images.push_back(zero(Fpx));       // x[0] |--> 0
    images.push_back(x[0] + 2*x[1]);   // x[1] |--> x[0] + 2*x[1]
    images.push_back(one(Fpx));        // x[2] |--> 1
    images.push_back(TwoThirds*x[3]);  // x[3] |--> (2/3)*x[3]

    RingHom phi = PolyAlgebraHom(Fpx, Fpx, images);

    cout << "We have built the following homomorphism"
         << " which we have called phi:" << endl
         << phi << endl
         << "Here are some values calculated using phi:" << endl
         << endl;

    cout << x[0] << " |--> " << phi(x[0]) << endl;
    cout << x[1] << " |--> " << phi(x[1]) << endl;
    cout << x[2] << " |--> " << phi(x[2]) << endl;
    cout << x[3] << " |--> " << phi(x[3]) << endl;
    cout << "f    |--> " << phi(f) << endl;
    cout << endl;

    ideal I(x[0]+x[1], x[2]-x[3]);

    const vector<RingElem>& h=gens(I);
    long NGens = len(h);
    vector<RingElem> v;
    v.reserve(NGens); // not necessary but good style: reserves memory (C++ STL)
    for (long i=0; i < NGens; ++i)  v.push_back(phi(h[i])); 
    cout << "v = " << v << endl;
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
    cerr << "***ERROR***  UNCAUGHT CoCoAError";
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingHom3.C,v 1.15 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingHom3.C,v $
// Revision 1.15  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.14  2019/10/29 16:02:08  bigatti
// -- added brackets in a comment
//
// Revision 1.13  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.12  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.11  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.10  2012/07/04 12:18:02  abbott
// Improved printing: actually use the indet names instead of hardwiring the names in a string.
//
// Revision 1.9  2012/05/04 09:28:45  bigatti
// -- slightly modified comment for reserve
//
// Revision 1.8  2012/03/30 16:19:00  bigatti
// -- some cleaning and a helfun comment for "reserve"
//
// Revision 1.7  2012/02/08 17:45:51  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.6  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.5  2010/02/04 11:19:46  bigatti
// -- added example of homomorphism on an ideal
//
// Revision 1.4  2009/06/13 00:35:21  abbott
// Changed name of "coordinate change hom" into phi.
//
// Revision 1.3  2009/06/11 14:09:39  abbott
// Minor reformatting to avoid short lines (because this file is
// used as an example with projectors, and the window is narrow).
//
// Revision 1.2  2007/06/21 21:31:05  abbott
// Cleaned up ex-RingHom3 for public presentation at CoCoA School 2007 (at RISC).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.12  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.11  2007/03/08 10:25:41  bigatti
// -- minor: removed space
//
// Revision 1.10  2007/03/07 14:52:56  bigatti
// -- minor: text and comments
//
// Revision 1.9  2007/03/07 14:32:15  bigatti
// -- now uses CanonicalHom and PolyAlgebraHom
// -- minor cleaning
//
// Revision 1.8  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.7  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.6  2007/02/26 17:39:46  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.5  2007/02/12 16:11:12  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.4  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.3  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.2  2006/06/21 12:10:33  cocoa
// -- simplified thanks to function "indets"
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.5  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.4  2005/09/30 13:02:48  cocoa
// -- minor changes
//
// Revision 1.3  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.2  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.1  2005/05/04 16:35:06  cocoa
// -- new examples
//
