// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how we define a ring homomorphism       \n"
  "to evaluate polynomials.                                   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  RingElem EvalAtZero(RingElem f) //--> the constant term
  {
    if (!IsSparsePolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "ImplicitDirectLPP");
    for (SparsePolyIter it= BeginIter(f); !IsEnded(it); ++it)
      if (IsOne(PP(it))) return coeff(it);
    return zero(CoeffRing(owner(f)));
  }
  

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Create some coefficient rings
    ring Fp = NewZZmod(32003);

    PolyRing Fpx = NewPolyRing(Fp, SymbolRange("x", 0,3)); // Fpx is Z/(32003)[x[0..3]]
    cout << "  --- PolyRing is " << Fpx << endl;

    const vector<RingElem>& x = indets(Fpx);
    RingElem f = 4*x[1]*x[1] + x[0]*x[2] - 2*x[0] +5;
    cout << "  --- f = " << f << endl;

    cout << "  --- EvalAtZero(f) = " << EvalAtZero(f) << endl;

    cout << "  --- Definition of PolyAlgebraHom  eval" << endl;
    RingHom eval = PolyAlgebraHom(Fpx, Fp, "1, 0, 2/3, 100");

    for (long i=0; i<4; ++i)
      cout << x[i] << " |-> " << eval(x[i]) << endl;
    cout << "f    |-> " << eval(f) << endl;

    cout << " --- Explicit construction of a C++ vector  images" << endl;
    vector<RingElem> images;
    images.push_back(one(Fp));
    images.push_back(zero(Fp));
    images.push_back(RingElem(Fp,"2/3")); // x[2] |-> 2/3
    // using string syntax, because C++ reads 2/3 as an integer division!
    images.push_back(RingElem(Fp,100)); // x[3] |-> 100

    eval = PolyAlgebraHom(Fpx, Fp, images);
    // same as PolyRingHom(Fpx, Fp, IdentityHom(Fp), images);
    cout << "f    |-> " << eval(f) << endl;

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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingHom4.C,v 1.14 2024/07/11 15:39:31 bigatti Exp $
// $Log: ex-RingHom4.C,v $
// Revision 1.14  2024/07/11 15:39:31  bigatti
// Summary: error codes: NotElemPolyRing to ReqElemPolyRing
//
// Revision 1.13  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.12  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.11  2020/01/13 16:38:53  bigatti
// -- improved style, better examples ... I think ;-)
//
// Revision 1.10  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.9  2017/04/27 15:28:00  abbott
// Summary: Minor cleaning
//
// Revision 1.8  2017/04/27 15:24:52  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.7  2017/04/21 13:57:03  bigatti
// -- added example for EvalAtZero
//
// Revision 1.6  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.5  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.4  2012/07/04 12:18:02  abbott
// Improved printing: actually use the indet names instead of hardwiring the names in a string.
//
// Revision 1.3  2012/02/08 17:45:51  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.8  2007/03/07 14:52:56  bigatti
// -- minor: text and comments
//
// Revision 1.7  2007/03/07 14:32:53  bigatti
// -- now uses PolyAlgebraHom
// -- minor cleaning
//
// Revision 1.6  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/26 17:39:46  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.3  2007/02/12 16:15:37  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.2  2006/04/27 16:12:07  cocoa
// -- changed NewIdentityRingHom --> NewIdentityHom
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.2  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.1  2005/09/30 13:01:34  cocoa
// -- first import
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
