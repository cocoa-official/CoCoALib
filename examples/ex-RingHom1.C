// Copyright (c) 2007  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "The example in this file shows how to create and use some \n"
  "homomorphisms between rings.                              \n";

const string LongDescription =
  "CanonicalHom is an easy way to make these homomorphisms:  \n"
  "R --> R/I,    R --> R[x],   R --> FractionField(R),       \n"
  "R --> R,      QQ --> R,     ZZ --> R,                     \n"
  "PolyAlgebraHom makes the R-algebra homomorphisms:         \n"
  "R[x] --> R,   R[x] --> R[y]                               \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Create some rings
    ring QQ = RingQQ();
    ring Fp = NewZZmod(101);
    PolyRing P = NewPolyRing(Fp, symbols("a,b"));   // P is Fp[a,b]

    cout << "-- CanonicalHom into P = " << P << endl;

    RingHom FpToP = CanonicalHom(Fp, P);  // same as CoeffEmbeddingHom(P)
    RingHom QQToP = CanonicalHom(QQ, P);  // same as QQEmbeddingHom(P)
    // NB!! QQToP is a partial homomorphism:
    // e.g. we cannot compute QQToP(one(QQ)/101)

    RingElem a = RingElem(Fp, 13);
    cout << "FpToP(" << a << ") = " << FpToP(a) << endl;

    RingElem q = RingElem(QQ, 5)/2;
    cout << "QQToP(" << q << ") = " << QQToP(q) << endl;
    cout << "  same as (RingElem calls CanonicalHom)" << endl;
    cout << "  RingElem(P,q) = " << RingElem(P,q) << endl;
    cout << endl;
 
    vector<RingElem> IndetImages;
    IndetImages.push_back(RingElem(Fp,2));
    IndetImages.push_back(RingElem(Fp,5));
    cout << "-- PolyAlgebraHom:    "
         << indet(P,0) << " |--> " << IndetImages[0] << "   "
         << indet(P,1) << " |--> " << IndetImages[1] << endl;

    RingHom PToFp = PolyAlgebraHom(P, Fp, IndetImages);
  
    RingElem f = 100*indet(P,0) + indet(P,1);
    cout << "PToFp(" << f << ") = " << PToFp(f) << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingHom1.C,v 1.9 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingHom1.C,v $
// Revision 1.9  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.8  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.7  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.6  2015/06/29 15:47:58  bigatti
// -- code in namespace CoCoA
//
// Revision 1.5  2013/03/15 15:19:26  bigatti
// -- added example for new RingElem ctor calling CanonicalHom
//
// Revision 1.4  2012/02/08 17:45:51  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2009/01/23 16:08:41  abbott
// Clarified a comment.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/08 10:24:16  bigatti
// -- simplest example
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
